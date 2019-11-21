# -*- coding: utf-8 -*-
from __future__ import division
import os
import numpy as np
import json
import scipy as sp
import xmltodict
import statsmodels.api as sm
from statsmodels.formula.api import ols
import requests, requests.exceptions
import shutil
import multiprocessing
import nibabel as nib
import logging
import util

"""
Find a set of differentially expressed genes between two user defined volumes of interest based on JuBrain maps. The tool downloads expression values of user specified sets of genes from Allen Brain API. Then, it uses zscores to find which genes are expressed differentially between the user specified regions of interests. This tool is available as a Python package.
Example:
    $ python driver.py
"""

def get_specimen_data(specimen_metadata):
    """
    For each specimen, extract the name and alignment matrix and stores a dict
    Args:
          specimen_metadata (dict): Contains metadata for specimens used in Allen Brain. They can be downloaded through Allen Brain API.
    Returns:
          dict: specimen dict with two keys -
                a) name: name of the specimen
                b) alignment3d: transformation matrix to convert from MRI to MNI152 space.
    """
    specimen = dict()
    specimen['name'] = specimen_metadata['name']
    x = specimen_metadata['alignment3d']
    specimen['alignment3d'] = np.array([
    [x['tvr_00'], x['tvr_01'], x['tvr_02'], x['tvr_09']],
    [x['tvr_03'], x['tvr_04'], x['tvr_05'], x['tvr_10']],
    [x['tvr_06'], x['tvr_07'], x['tvr_08'], x['tvr_11']],
    [0, 0, 0, 1]])
    return specimen

def unwrap_self_do_anova_with_permutation_rep(*args, **kwargs):
    """
    Helper function to enable usage of the multiprocessing module inside a class.
    Args:
          *args: Variable length argument list.
          **kwargs: Arbitrary keyword arguments.
    Returns:
              do_anova_with_permutation_rep()
    """
    return Analysis.do_anova_with_permutation_rep(*args, **kwargs)


class Analysis:

    def __init__(self, gene_cache_dir, filter_threshold=0.2, single_probe_mode=False, verbose=False, n_rep=1000):
        """
        Initialize the Analysis class with various internal variables -
        Args:
            gene_cache_dir (str): Disk location where the gene and specimen data, downloaded from Allen Brain API, will be cached for future reuse.
            verbose (bool): True for verbose output during execution of the code.
        Attributes:
            probe_ids (list) : list of probe ids associated with the given list of genesymbols which are not present in gene_cache, if it exists.
            gene_list (list) : list of genesymbols to perform differential analysis with, provided by the user.
            gene_list_to_download (list) : list of genesymbols whose information is not in the cache yet, needs to be downloaded.
            gene_symbols (list) :  same length as probe_ids, each probe is represented with the corresponding genesymbol, used by get_mean_zscores().
            gene_cache (dict) :  dictionary to indicate which genes are present in the cache.
            donor_ids (list) :  six donor ids of Allen Brain API.
            allen_brain_api_data (dict) : dictionary to store Allen Brain API data, with two keys - samples_and_zscores contain samples and zscores from allen brain  api for all the given probes and specimen_info contains age, race, sex of the six donors represented by donor_ids.
            rois (list) : list of two nii volumes for each region of interest used in differential analysis.
            filtered_coords_and_zscores (list) : internal variable for storing MNI152 coordinates and zscores corresponsing to each region of interest.
            filter_threshold (float) : internal variable used at filter_coordinates_and_zscores() to select or reject a sample.
            n_rep (int) : number of iterations of FWE correction.
            cache (str) : disk location where data from Allen Brain API has been downloaded and stored.
            anova_factors (dict) : internal dictionary used by the ANOVA module - contains five keys - 'Age', 'Race', 'Specimen', 'Area', 'Zscores'.
            genesymbol_and_mean_zscores (dict) : dictionary with two keys - uniqueid, combinedzscores where each gene and the winsorzed mean zscores over all probes associated with that gene is stored.
            gene_id_and_pvalues (dict) : dict for storing gene ids and associated p values.
        """
        self.probe_keys = []
        self.probe_ids = []
        self.gene_list = []
        self.gene_list_to_download = []
        self.gene_symbols = []
        self.gene_cache = {}
        self.donor_ids = ['15496', '14380', '15697', '9861', '12876', '10021'] #HARDCODING donor_ids

        # depends on donor id
        self.samples_zscores_and_specimen_dict = dict.fromkeys(['samples_and_zscores', 'specimen_info'])
        self.specimen_factors = dict.fromkeys(['id', 'name', 'race', 'gender', 'age'])
        self.samples_zscores_and_specimen_dict['specimen_info'] = []
        self.samples_zscores_and_specimen_dict['samples_and_zscores'] = []
        self.rois = []
        self.affine = []
        self.filtered_coords_and_zscores = []
        self.filter_threshold = float(filter_threshold)
        self.n_rep = n_rep
        self.cache_dir = gene_cache_dir
        self.verbose = verbose        
        self.single_probe_mode = single_probe_mode
        self.anova_factors = dict.fromkeys(['Age', 'Race', 'Specimen', 'Area', 'Zscores'])
        self.genesymbol_and_mean_zscores = dict.fromkeys(['uniqueId', 'combined_zscores'])
        self.result_for_web_version = []
        logging.basicConfig(level=logging.INFO)

        '''
        Removes any folder that has not been written to properly, most likely due to force quit
        '''
        probe_path = os.path.join(self.cache_dir, '{}/probes.txt'.format(self.donor_ids[0]))
        if os.path.exists(self.cache_dir) and not os.path.exists(probe_path):
            shutil.rmtree(self.cache_dir, ignore_errors = False)

        # wait, if cache path doesn't exist, should we create the cache path?
        # also, it is a good idea to download required gene probes etc in background
        '''
        Creates a gene cache to indicate which genes are present in the cache
        '''
        if os.path.exists(self.cache_dir):
            self.create_gene_cache()
            logging.getLogger(__name__).info('{} genes exist in {}'.format(len(self.gene_cache), self.cache_dir))
        else:
            logging.getLogger(__name__).info('{} does not exist. It will take some time '.format(self.cache_dir))

    def DifferentialAnalysis(self, gene_list, roi1, roi2):
        """
        Driver routine
        Args:
              gene_list (list): list of gene symbols to perform differential analysis with, provided by the user.
              roi1 (nibabel.nifti1.Nifti1Image): Nifti1Image representing the probability map of the first region of interest, chosen by the user.
              roi2 (nibabel.nifti1.Nifti1Image): Nifti1Image representing the probability map of the second region of interest, chosen by the user.
        Returns:
             dict: A dictionary representing the gene symbols and their corresponding p values
        """
        if not gene_list:
            raise ValueError('Atleast one gene is needed for the analysis')
        if not (isinstance(roi1['data'], nib.nifti1.Nifti1Image) and isinstance(roi2['data'], nib.nifti1.Nifti1Image)):
            raise ValueError('Atleast two valid regions of interest are needed')
        self.set_candidate_genes(gene_list)
        self.set_roi_MNI152(roi1, 0)
        self.set_roi_MNI152(roi2, 1)
        self.affine = roi1['data'].affine
        logging.getLogger(__name__).info('Starting the analysis. This may take some time.....')
        self.anova()
        #return self.gene_id_and_pvalues
        return json.dumps(self.result_for_web_version)

    # only creates gene cache for the first donor_id
    # and populates gene_cache with "MAOA": None
    def create_gene_cache(self):
        """
        Create a dictionary with an entry for each gene whose api information has been downloaded. We scan self.cache_dir/15496/probes.txt for the gene symbols.
        """
        with open(os.path.join(self.cache_dir, '{}/probes.txt'.format(self.donor_ids[0])), 'r') as f:
            probes = json.load(f)
        for probe in probes:
            self.gene_cache.update({probe['gene-symbol'] : None})
        if self.verbose:
            logging.getLogger(__name__).info('gene_cache: {}'.format(self.gene_cache))

    def retrieve_probe_ids(self):
        """
        Retrieve probe ids for the given gene lists, update self.probe_ids which will be used by download_and_save_zscores_samples() or download_and_save_zscores_samples_partial() to
        form the url and update self.gene_symbols to be used by get_mean_zscores()
        """

        for gene in self.gene_list:
            
            data = util.from_brainmap_retrieve_gene(gene=gene, verbose=self.verbose)

            if int(data['Response']['@num_rows']) <= 0:
                raise ValueError('Please check the spelling of {}. No such gene exists in Allen Brain API.'.format(gene))
            self.probe_keys = self.probe_keys + [donor['id'] for donor in data['Response']['probes']['probe']]
            self.probe_ids = self.probe_ids + [donor['id'] for donor in data['Response']['probes']['probe'] if gene in self.gene_list_to_download]
            self.gene_symbols = self.gene_symbols + [gene for donor in data['Response']['probes']['probe']]
        if self.verbose:
            logging.getLogger(__name__).info('probe_ids: {}'.format(self.probe_ids))
            logging.getLogger(__name__).info('gene_symbols: {}'.format(self.gene_symbols))


    def read_cached_zscores_samples_and_specimen_data(self):
        """
        Read cached Allen Brain Api data from disk location and update self.samples_zscores_and_specimen_dict['specimen_info'] and self.samples_zscores_and_specimen_dict['samples_and_zscores']
        """
        probe_keys_dict = {}
        for p in self.probe_keys:
            probe_keys_dict.update({int(p):[]})
        for donor in self.donor_ids:
            donor_path = os.path.join(self.cache_dir, donor)
            specimen = dict.fromkeys(['name', 'alignment3d'])
            with open(os.path.join(donor_path, 'specimenMat.txt'), 'r') as f:
                specimen['alignment3d'] = np.loadtxt(f)
            with open(os.path.join(donor_path, 'specimenName.txt'), 'r') as f:
                specimen['name'] = f.read()
            self.samples_zscores_and_specimen_dict['specimen_info'] = self.samples_zscores_and_specimen_dict['specimen_info'] + [specimen]
            with open(os.path.join(donor_path, 'samples.txt'), 'r') as f:
                samples = json.load(f)
            with open(os.path.join(donor_path, 'probes.txt'), 'r') as f:
                probes = json.load(f)
            newprobes = [p for p in probes if p['id'] in probe_keys_dict.keys()]
            zscores = np.array([[float(p['z-score'][j]) for p in newprobes] for j in range(len(samples))])
            self.samples_zscores_and_specimen_dict['samples_and_zscores'] = self.samples_zscores_and_specimen_dict['samples_and_zscores']  + [{'samples' : samples, 'zscores' : zscores}]
            if self.verbose:
                logging.getLogger(__name__).info('inside readcachedata {}, {}'.format(len(self.samples_zscores_and_specimen_dict['samples_and_zscores'][-1]['samples']), self.samples_zscores_and_specimen_dict['samples_and_zscores'][-1]['zscores'].shape))

    def set_roi_MNI152(self, roi, index):
        """
        For each specimen and for the given roi populate self.filtered_coords_and_zscores with zscores and valid coordinates
        Args:
              roi (nib.nifti1.Nifti1Image): probability map of a region of interest.
              index (int): 0 or 1, representing which region of interest it is.
        """
        if index < 0 or index > 1:
            raise ValueError('only 0 and 1 are valid choices')
        for i in range(len(self.samples_zscores_and_specimen_dict['specimen_info'])):
            self.filtered_coords_and_zscores.append(
                util.filter_coordinates_and_zscores(
                    roi_name=roi['name'],
                    roi_nii=roi['data'],
                    index_to_samples_zscores_and_specimen_dict=self.samples_zscores_and_specimen_dict['samples_and_zscores'][i],
                    specimen=self.samples_zscores_and_specimen_dict['specimen_info'][i],
                    index=index,
                    filter_threshold=self.filter_threshold
                )
            )

    def __download_and_save_zscores_and_samples(self, donor_id):
        text = util.from_brainmap_retrieve_microarray_filterby_donorids_probeids(probe_ids=self.probe_ids, donor_id=donor_id)
        data = text['msg']
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        donor_path = os.path.join(self.cache_dir, donor_id)
        if not os.path.exists(donor_path):
            os.makedirs(donor_path)

        zscores = np.array([[float(data['probes'][i]['z-score'][j]) for i in range(len(data['probes']))] for j in range(len(data['samples']))])

        with open(os.path.join(donor_path, 'zscores.txt'), 'wb') as f:
            np.savetxt(f, zscores, fmt = '%.5f')
        with open(os.path.join(donor_path, 'samples.txt'), 'w') as outfile:
            json.dump(data['samples'], outfile)
        with open(os.path.join(donor_path, 'probes.txt'), 'w') as outfile:
            json.dump(data['probes'], outfile)

        if self.verbose:
            logging.getLogger(__name__).info('For {} samples_length: {}  probes_length: {} zscores_shape: {} '.format(donor_id,len(data['samples']),len(data['probes']), zscores.shape))
        return {'samples' : data['samples'], 'probes' : data['probes'], 'zscores' : zscores}

    def __download_and_save_zscores_and_samples_partial(self, donor_id):
        text = util.from_brainmap_retrieve_microarray_filterby_donorids_probeids(probe_ids=self.probe_ids, donor_id=donor_id)
        data = text['msg']
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        donor_path = os.path.join(self.cache_dir, donor_id)
        if not os.path.exists(donor_path):
            os.makedirs(donor_path)
        probes = data['probes']
        zscores = np.array([[float(data['probes'][i]['z-score'][j]) for i in range(len(data['probes']))] for j in range(len(data['samples']))])
        with open(os.path.join(donor_path, 'probes.txt'), 'r') as f:
            probes_cached = json.load(f)
        with open(os.path.join(donor_path, 'zscores.txt'), 'r') as f:
            zscores_cached = np.loadtxt(f)
        with open(os.path.join(donor_path, 'samples.txt'), 'w') as outfile:
            json.dump(data['samples'], outfile)
        probes = probes_cached + probes
        with open(os.path.join(donor_path, 'probes.txt'), 'w') as outfile:
            json.dump(probes, outfile)
        zscores = np.append(zscores_cached, zscores, axis=1)
        np.savetxt(os.path.join(donor_path, 'zscores.txt'), zscores)

    def download_and_save_specimens(self):
        """
        Download names and transformation matrix for each specimen/donor from Allen Brain Api and save them on disk as specimenName.txt
        and specimenMat.txt respectively, load.
        """
        specimens  = ['H0351.1015', 'H0351.1012', 'H0351.1016', 'H0351.2001', 'H0351.1009', 'H0351.2002']
        self.samples_zscores_and_specimen_dict['specimen_info'] = []
        for specimen_id in specimens:
            text = util.from_brainmap_retrieve_specimen(specimen_id, verbose=self.verbose)
            self.samples_zscores_and_specimen_dict['specimen_info'] = self.samples_zscores_and_specimen_dict['specimen_info'] + [get_specimen_data(text['msg'][0])]
        if self.verbose:
            logging.getLogger(__name__).info('{}'.format(self.samples_zscores_and_specimen_dict['specimen_info']))

        for donor, specimen in zip(self.donor_ids, self.samples_zscores_and_specimen_dict['specimen_info']):
            with open(os.path.join(self.cache_dir, '{}/specimenName.txt'.format(donor)), 'w') as outfile:
                outfile.write(specimen['name'])
            np.savetxt(os.path.join(self.cache_dir, '{}/specimenMat.txt'.format(donor)), specimen['alignment3d'])


    def download_and_save_zscores_and_samples(self):
        """
        Loop through the donors and call download_and_save_zscores_samples() and populate samples_zscores_and_specimen_dict
        """
        self.samples_zscores_and_specimen_dict['samples_and_zscores'] = [self.__download_and_save_zscores_and_samples(donor) for donor in self.donor_ids]

    def download_and_save_zscores_samples_and_specimen_data(self):
        """
        Download data from Allen Brain Api for the given set of genes and specimen
        """
        self.download_and_save_zscores_and_samples()
        self.download_and_save_specimens()

    def set_candidate_genes(self, gene_list):
        """
        Set list of genesymbols and prepare to read/download data for them.
        Args:
              gene_list (list) : list of genesymbols to perform differential analysis with, provided by the user.
        """
        self.gene_list = gene_list
        donorprobe = os.path.join(self.cache_dir, '{}/probes.txt'.format(self.donor_ids[0]))
        '''
        If the cache doesnt exist then get all the probes associated with the genes and download and save the api and specimen information
        and populate samples_zscores_and_specimen_dict['samples_and_zscores'] and allen_brain_api_data['specimen_info'].
        '''
        if not os.path.exists(self.cache_dir):
            self.gene_list_to_download = self.gene_list[:]
            self.retrieve_probe_ids()
            self.download_and_save_zscores_samples_and_specimen_data()
        else:
            '''
            If the cache exists there are two possibilities.
            a) All the requested genes are present in the cache. In that case, read the downloaded api and specimen data from disk and
            populate samples_zscores_and_specimen_dict['samples_and_zscores'] and samples_zscores_and_specimen_dict['specimen_info']
            b) A few of the requested genes are present in the cache. In that case, read the downloaded api and specimen data from disk and
            populate samples_zscores_and_specimen_dict['samples_and_zscores'] and samples_zscores_and_specimen_dict['specimen_info'] for genes already present in the cache. For the rest, get all the probes
            associated with the genes, download and save the api and specimen information and populate samples_zscores_and_specimen_dict['samples_and_zscores'] and
            samples_zscores_and_specimen_dict['specimen_info'].
            '''
            self.gene_list_to_download = [key for key in self.gene_list if key not in self.gene_cache.keys()]
            if self.gene_list_to_download:
                logging.getLogger(__name__).info('Microarray expression values of {} gene(s) need(s) to be downloaded'.format(len(self.gene_list_to_download)))
            if self.verbose:
                logging.getLogger(__name__).info('genes to be downloaded:{} '.format(self.gene_list_to_download))
            self.retrieve_probe_ids()
            if self.gene_list_to_download:
                for donor in self.donor_ids:
                    self.__download_and_save_zscores_and_samples_partial(donor)            
            self.read_cached_zscores_samples_and_specimen_data()


    def get_mean_zscores(self, combined_zscores):
        """
        Compute Winsorzed mean of zscores over all probes associated with a given gene. combined_zscores have zscores for all the probes and all the valid coordinates.
        As a gene_id_and_pvalues you get a numpy array of size len(self.filtered_coords_and_zscores)xlen(self.gene_list). self.genesymbol_and_mean_zscores['combined_zscores'][i][j] returns the     winsorzed mean of jth gene taken over all the probes corresponding to the ith valid sample.
        Args:
             combined_zscores (list): lists of zscores corresponding to each region of interest, populated from filtered_coords_and_zscores
        """
        unique_gene_symbols = np.unique(self.gene_symbols)
        '''
        A = [a,a,a,b,b,b,c,c]
        B = [a,b,c]
        Following line of code will give  indices = [[0,1,2],[3,4,5],[6,7]]
        '''
        indices = [np.where(np.in1d(self.gene_symbols, genesymbol))[0] for genesymbol in unique_gene_symbols]
        '''
        for i in range (len(unique_gene_symbols)):
                for j in range(len(combined_zscores)):
                    for k in range(len(indices[i])):
                        tmp[j] = combined_zscores[j][indices[i][k]][:]
                winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1))
        '''
        winsorzed_mean_zscores =  np.array([[np.mean(sp.stats.mstats.winsorize([combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))], limits=0.1)) for i in range (len(unique_gene_symbols))] for j in range(len(combined_zscores))])
        self.genesymbol_and_mean_zscores['uniqueId'] = unique_gene_symbols
        self.genesymbol_and_mean_zscores['combined_zscores'] = winsorzed_mean_zscores

    def accumulate_roicoords_and_name(self):
        """
        Populate roi names, coordinates and sample well and polygon id for display
        """
        areainfo = {}
        for roi_coord_zscore in self.filtered_coords_and_zscores:
            key = roi_coord_zscore['realname']
            if key not in areainfo:
                areainfo[key] = []
            i = 0
            for c in roi_coord_zscore['coords']:
                if self.single_probe_mode:
                    areainfo[key].append({'xyz' : np.transpose(np.matmul(self.affine,np.transpose(np.append(c,1))))[0:3].tolist(), 'winsorzed_mean' : self.combined_zscores[i].tolist()})
                else:
                    areainfo[key].append({'xyz' : np.transpose(np.matmul(self.affine,np.transpose(np.append(c,1))))[0:3].tolist(), 'winsorzed_mean' : self.genesymbol_and_mean_zscores['combined_zscores'][i].tolist()})
                #areainfo[key].append({'xyz' : c.tolist(), 'winsorzed_mean' : self.genesymbol_and_mean_zscores['combined_zscores'][i].tolist()})
                i = i+1

        self.result_for_web_version.append(areainfo)

    def initialize_anova_factors(self):
        """
        Prepare self.anova_factors. Populate Age, Race, Area, Specimen, Zcores keys of self.anova_factors
        """
        #Populates self.genesymbol_and_mean_zscores (uniqueid and zscores)
        if self.single_probe_mode:
            self.combined_zscores = np.array([roi_coord_and_zscore['zscores'][i] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))])
            self.n_genes = len(self.combined_zscores[0])
        else:
            combined_zscores = [roi_coord_and_zscore['zscores'][i] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
            self.get_mean_zscores(combined_zscores)
            self.n_genes = len(self.genesymbol_and_mean_zscores['combined_zscores'][0])
        #combined_zscores = [roi_coord_and_zscore['zscores'][i] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
        #Populates self.specimenFactors (id, race, gender, name, age)
        self.read_specimen_factors(self.cache_dir)
        if self.verbose:
            logging.getLogger(__name__).info("number of specimens: {} name: {}".format(len(self.specimen_factors), len(self.specimen_factors['name'])))

        self.anova_factors['Area'] = [roi_coord_and_zscore['name'] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
        self.anova_factors['Specimen'] = [roi_coord_and_zscore['specimen'] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
        '''
        Both Age and Race should have len(self.filtered_coords_and_zscores) entries. The following three lines are used to get the correct values from specimenFactors['Age'] and specimenFactors['Race'] using specimenFactors['name'] and repeat them the required number of times as given by self.anova_factors['Specimen']
        '''
        st = set(self.specimen_factors['name'])
        self.anova_factors['Age'] = [self.specimen_factors['age'][self.specimen_factors['name'].index(specimen_name)] for ind, specimen_name in enumerate(self.anova_factors['Specimen'])]
        self.anova_factors['Race'] = [self.specimen_factors['race'][self.specimen_factors['name'].index(specimen_name)] for ind, specimen_name in enumerate(self.anova_factors['Specimen'])]
        self.accumulate_roicoords_and_name()

    def first_iteration(self):
        """
        Perform one iteration of ANOVA. Use output of this to populate F_vec_ref_anovan which becomes initial estimate of n_rep passes of FWE.
        """

        self.F_vec_ref_anovan = np.zeros(self.n_genes)
        for i in range(self.n_genes):
            if self.single_probe_mode:
                self.anova_factors['Zscores'] = self.combined_zscores[:,i]
            else:
                self.anova_factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,i]
            #self.anova_factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,i]
            mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.anova_factors).fit()
            aov_table = sm.stats.anova_lm(mod, typ=1)
            if self.verbose:
                logging.getLogger(__name__).info('aov table: {}'.format(aov_table))
            #F_vec_ref_anovan is used as an initial condition to F_mat_perm_anovan in fwe_correction
            self.F_vec_ref_anovan[i] = aov_table['F'][0]

    def do_anova_with_permutation_gene(self, index_to_gene_list):
        """
        Perform one repetition of anova for each gene
        Args:
              index_to_gene_list (int) : Index into the genesymbol_and_mean_zscores['combined_zscores'] array, representing mean zscore of a gene.
        Returns:
                 float: F value extracted from the anova table
        """
        self.anova_factors['Area'] = np.random.permutation(self.anova_factors['Area'])
        if self.single_probe_mode:
            self.anova_factors['Zscores'] = self.combined_zscores[:,index_to_gene_list]
        else:
            self.anova_factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,index_to_gene_list]
            # self.anova_factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,index_to_gene_list]
        mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.anova_factors).fit()
        aov_table = sm.stats.anova_lm(mod, typ=1)
        return aov_table['F'][0]

    def do_anova_with_permutation_rep(self):
        """
        Perform one repetition of anova for all genes
        Returns:
                 list: a list of F_values, one for each gene.
        """
        return list(map(self.do_anova_with_permutation_gene, range(0,self.n_genes)))


    def fwe_correction(self):
        """
        Perform n_rep passes of FWE using gene_id_and_pvalues of first_iteration() as an initial guess
        """
        invn_rep = 1/self.n_rep
        initial_guess_F_vec = self.F_vec_ref_anovan
        pool = multiprocessing.Pool()
        #self.F_mat_perm_anovan = np.array(pool.map(unwrap_self_do_anova_with_permutation_rep, zip([self]*self.n_rep, range(1,self.n_rep)) #for parameter
        self.F_mat_perm_anovan = np.array(pool.map(unwrap_self_do_anova_with_permutation_rep, [self]*(self.n_rep-1)))
        self.F_mat_perm_anovan = np.insert(self.F_mat_perm_anovan, 0, initial_guess_F_vec, axis=0)
        self.accumulate_gene_id_and_pvalues()

    def div_func(self, arr):
        """
        Helper function to compute average F-value, averaged over self.n_rep for each gene
        Args:
              arr (numpy.ndarray) : F-values for all the repetition for a single gene
        Returns:
                float: average F-value for that gene
        """
        return np.count_nonzero(arr)/self.n_rep

    def accumulate_gene_id_and_pvalues(self):
        """
        Populate pvalues and geneids for the gene_id_and_pvalues dict after n_rep passes of FWE
        """
        #ref represenets maximum p value for each gene across n_rep repetitions
        #compute family wise error corrected p value
        self.FWE_corrected_p = np.apply_along_axis(self.div_func, 0, self.F_mat_perm_anovan.max(1)[:, np.newaxis] >= np.array(self.F_vec_ref_anovan))
        #self.FWE_corrected_p =  [np.count_nonzero(max_F_mat_perm_anovan >= f)/self.n_rep for f in self.F_vec_ref_anovan]
        if self.single_probe_mode:
            self.gene_id_and_pvalues = dict(zip(self.probe_keys, self.FWE_corrected_p))
        else:
            self.gene_id_and_pvalues = dict(zip(self.genesymbol_and_mean_zscores['uniqueId'], self.FWE_corrected_p))
        #self.gene_id_and_pvalues = dict(zip(self.genesymbol_and_mean_zscores['uniqueId'], self.FWE_corrected_p))
        self.result_for_web_version.append(self.gene_id_and_pvalues)
        if self.single_probe_mode:
            self.result_for_web_version.append(self.probe_keys)
        if self.verbose:
            logging.getLogger(__name__).info('gene_id_and_pvalues: {}'.format(self.gene_id_and_pvalues))

    def anova(self):
        """
        Perform one way anova on zscores as the dependent variable and specimen factors such as age, race, name and area
        as independent variables
        """

        self.initialize_anova_factors()
        self.first_iteration()
        self.fwe_correction()

    def build_specimen_factors(self, cache):
        """
        Download various factors such as age, name, race, gender of the six specimens from Allen Brain Api, save them at cache/specimenFactors.txt and create a dict.
        Args:
             cache (str): Location where the specimen_factors dict will be stored.
        """
        url_build_specimen_factors = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"
        try:
            text = requests.get(url_build_specimen_factors).json()
        except requests.exceptions.RequestException as e:
            logging.getLogger(__name__).error(e)
            raise
        factorPath = os.path.join(cache, 'specimenFactors.txt')
        with open(factorPath, 'w') as outfile:
            json.dump(text, outfile)
        res = text['msg']

        self.specimen_factors['id'] = [r['id'] for r in res]
        self.specimen_factors['name'] = [r['name'] for r in res]
        self.specimen_factors['race'] = [r['race_only'] for r in res]
        self.specimen_factors['gender'] = [r['sex'] for r in res]
        self.specimen_factors['age'] = [r['age']['days']/365 for r in res]


    def read_specimen_factors(self, cache):
        """
        Read various factors such as age, name, race, gender of the six specimens from cache/specimenFactors.txt.
        Args:
             cache (str): Location where the specimen_factors dict is stored.
        """
        fileName = os.path.join(cache, 'specimenFactors.txt')
        if not os.path.exists(fileName):
            self.build_specimen_factors(cache)
        f = open(fileName, "r")
        content = json.load(f)
        f.close()
        res = content['msg']
        self.specimen_factors = dict()
        self.specimen_factors['id'] = [r['id'] for r in res]
        self.specimen_factors['name'] = [r['name'] for r in res]
        self.specimen_factors['race'] = [r['race_only'] for r in res]
        self.specimen_factors['gender'] = [r['sex'] for r in res]
        self.specimen_factors['age'] = [r['age']['days']/365 for r in res]
