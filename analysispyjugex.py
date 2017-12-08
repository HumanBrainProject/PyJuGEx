# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
import os, glob
import numpy as np
from numpy import *
import json
from numpy.linalg import inv
import scipy as sp
import scipy.stats.mstats
import xmltodict
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats
import sys
import requests, requests.exceptions
import pandas as pd
import shutil
import multiprocessing
import nibabel as nib
import logging

def get_specimen_data(specimen_metadata):
    """
    For each specimen, extract the name and alignment matrix and put into a dict object
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

def transform_samples_MRI_to_MNI52(samples, transformation_mat):
    """
    Convert the MRI coordinates of samples to MNI152 space
    """
    np_T = np.array(transformation_mat[0:3, 0:4])
    mri = np.vstack(s['sample']['mri'] for s in samples)
    add = np.ones((len(mri), 1), dtype=np.int)
    mri = np.append(mri, add, axis=1)
    mri = np.transpose(mri)
    coords = np.matmul(np_T, mri)
    coords = coords.transpose()
    return coords

class Analysis:

    def __init__(self, gene_cache, verbose=False):
        """
        initialize_anova_factors the Analysis class with various internal variables -

        probe_ids = list of probe ids associated with the give list of genes which are not present in gene_cache, if it exists.
        gene_list = given list of genes
        gene_list_to_download = list of genes whose information is not in the cache yet, needs to be downloaded
        gene_symbols = same length as probe_ids, each probe is represented with the corresponding genesymbol, used by get_mean_zscores()
        gene_cache = dictionary to indicate which keys are present in the cache
        donor_ids = six donor ids of Allen Brain API
        allen_brain_api_data = dictionary to store Allen Brain API data, with two keys - samples_and_zscores contain samples and zscores from allen brain  api for all the given probes and specimen_info contains age, race, sex of the six donors represented by donor_ids
        rois = list of two nii volumes for each region of interest
        filtered_coords_and_zscores = Internal variable for storing MNI52 coordinates and zscores corresponsing to each region of interest.
        filter_threshold = Internal variable used at filter_coordinates_and_zscores() to select or reject a sample
        n_rep = number of iterations of FWE correction
        cache = Disk location where data from Allen Brain API has been downloaded and stored.
        anova_factors = Internal variable used by the ANOVA module - contains five keys - 'Age', 'Race', 'Specimen', 'Area', 'Zscores'
        genesymbol_and_mean_zscores = dictionary with two keys - uniqueid, combinedzscores where each gene and the winsorzed mean zscores over all probes associated with that gene is stored
        gene_id_and_pvalues = dict for storing gene ids and associated p values.
        """
        self.probe_ids = []
        self.gene_list = []
        self.gene_list_to_download = []
        self.gene_symbols = []
        self.gene_cache = {}
        self.donor_ids = ['15496', '14380', '15697', '9861', '12876', '10021'] #HARDCODING donor_ids
        self.allen_brain_api_dict = dict.fromkeys(['samples_and_zscores', 'specimen_info'])
        self.specimen_factors = dict.fromkeys(['id', 'name', 'race', 'gender', 'age'])
        self.allen_brain_api_dict['specimen_info'] = []
        self.allen_brain_api_dict['samples_and_zscores'] = []
        self.rois = []
        self.filtered_coords_and_zscores = []
        self.filter_threshold = 0.2
        self.n_rep = 1000
        self.cache_dir = gene_cache
        self.verbose = verbose
        self.anova_factors = dict.fromkeys(['Age', 'Race', 'Specimen', 'Area', 'Zscores'])
        self.genesymbol_and_mean_zscores = dict.fromkeys(['uniqueId', 'combined_zscores'])
        logging.basicConfig(level=logging.INFO)
        '''
        Removes any folder that has not been written to properly, most likely due to force quit
        '''
        probe_path = os.path.join(self.cache_dir, '{}/probes.txt'.format(self.donor_ids[0]))
        if os.path.exists(self.cache_dir) and not os.path.exists(probe_path):
            shutil.rmtree(self.cache_dir, ignore_errors = False)
        '''
        Creates a gene cache to indicate which genes are present in the cache
        '''
        if not os.path.exists(self.cache_dir):
            logging.getLogger(__name__).info('{} does not exist. It will take some time '.format(self.cache_dir))
        else:
            self.create_gene_cache()
            logging.getLogger(__name__).info('{} genes exist in {}'.format(len(self.gene_cache), self.cache_dir))

    def DifferentialAnalysis(self, gene_list, roi1, roi2):
        """
        Driver routine
        """
        if not gene_list:
            raise ValueError('Atleast one gene is needed for the analysis')
        if not (isinstance(roi1, nib.nifti1.Nifti1Image) or isinstance(roi2, nib.nifti1.Nifti1Image)):
            raise ValueError('Atleast two valid regions of interest are needed')
        self.set_candidate_genes(gene_list)
        self.set_roi_MNI152(roi1, 0)
        self.set_roi_MNI152(roi2, 1)
        logging.getLogger(__name__).info('Starting the analysis. This may take some time.....')
        self.anova()
        return self.gene_id_and_pvalues

    def create_gene_cache(self):
        """
        Create a dictionary with an entry for each gene whose api information has been downloaded.
        We scan self.cache_dir/15496/probes.txt for the gene symbols
        """
        with open(os.path.join(self.cache_dir, '{}/probes.txt'.format(self.donor_ids[0])), 'r') as f:
            probes = json.load(f)
        for probe in probes:
            self.gene_cache.update({probe['gene-symbol'] : None})
        if self.verbose:
            logging.getLogger(__name__).info('gene_cache: {}'.format(self.gene_cache))

    def retrieve_probe_ids(self):
        """
        Retrieve probe ids for the given gene lists, update self.probe_ids which will be used by query_allen_brain_api() or query_allen_brain_api_partial() to
        form the url and update self.gene_symbols to be used by get_mean_zscores()
        """
        base_retrieve_probe_ids = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq"
        end_retrieve_probe_ids = "],rma::options[only$eq'probes.id']"

        for gene in self.gene_list:
            url = '{base}{gene_symbol}{end}'.format(base = base_retrieve_probe_ids, gene_symbol = gene, end = end_retrieve_probe_ids)
            if self.verbose:
                logging.getLogger(__name__).info('url: {}'.format(url))
                #print(url)
            try:
                response = requests.get(url)
            except requests.exceptions.RequestException as e:
                logging.getLogger(__name__).error(e)
                raise
            data = xmltodict.parse(response.text)
            self.probe_ids = self.probe_ids + [donor['id'] for donor in data['Response']['probes']['probe'] if gene in self.gene_list_to_download]
            self.gene_symbols = self.gene_symbols + [gene for donor in data['Response']['probes']['probe']]

        if self.verbose:
            logging.getLogger(__name__).info('probe_ids: {}'.format(self.probe_ids))
            logging.getLogger(__name__).info('gene_symbols: {}'.format(self.gene_symbols))


    def read_cached_api_specimen_data(self):
        """
        Read cached Allen Brain Api data from disk location and update self.allen_brain_api_dict['specimen_info'] and self.allen_brain_api_dict['samples_and_zscores']
        """
        for donor in self.donor_ids:
            donorpath = os.path.join(self.cache_dir, donor)
            specimen = dict.fromkeys(['name', 'alignment3d'])
            with open(os.path.join(donorpath, 'specimenMat.txt'), 'r') as f:
                specimen['alignment3d'] = np.loadtxt(f)
            with open(os.path.join(donorpath, 'specimenName.txt'), 'r') as f:
                specimen['name'] = f.read()
            self.allen_brain_api_dict['specimen_info'] = self.allen_brain_api_dict['specimen_info'] + [specimen]
            with open(os.path.join(donorpath, 'samples.txt'), 'r') as f:
                samples = json.load(f)
            with open(os.path.join(donorpath, 'zscores.txt'), 'r') as f:
                zscores = np.loadtxt(f)
            self.allen_brain_api_dict['samples_and_zscores'] = self.allen_brain_api_dict['samples_and_zscores']  + [{'samples' : samples, 'zscores' : zscores}]
            if self.verbose:
                logging.getLogger(__name__).info('inside readcachedata {}, {}'.format(len(self.allen_brain_api_dict['samples_and_zscores'][-1]['samples']), self.allen_brain_api_dict['samples_and_zscores'][-1]['zscores'].shape))

    def set_roi_MNI152(self, voi, index):
        """
        For each specimen and for the given voi populate self.filtered_coords_and_zscores with zscores and valid coordinates
        """
        if index < 0 or index > 1:
            raise ValueError('only 0 and 1 are valid choices')
        for i in range(len(self.allen_brain_api_dict['specimen_info'])):
            self.filtered_coords_and_zscores.append(self.filter_coordinates_and_zscores(voi, self.allen_brain_api_dict['samples_and_zscores'][i], self.allen_brain_api_dict['specimen_info'][i], index))

    def query_allen_brain_api(self, donorId):
        """
        Query Allen Brain Api for given set of genes
        """
        base_query_api = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        probes = ''.join('{},'.format(probe) for probe in self.probe_ids)[:-1]
        end_query_api = "][donors$eq{}]".format(donorId)
        url = '{}{}{}'.format(base_query_api, probes, end_query_api)
        try:
            response = requests.get(url)
            text = requests.get(url).json()
        except requests.exceptions.RequestException as e:
            logging.getLogger(__name__).info(e)
            raise
        data = text['msg']
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        donorPath = os.path.join(self.cache_dir, donorId)
        if not os.path.exists(donorPath):
            os.makedirs(donorPath)

        zscores = np.array([[float(data['probes'][i]['z-score'][j]) for i in range(len(data['probes']))] for j in range(len(data['samples']))])

        with open(os.path.join(donorPath, 'zscores.txt'), 'wb') as f:
            np.savetxt(f, zscores, fmt = '%.5f')
        with open(os.path.join(donorPath, 'samples.txt'), 'w') as outfile:
            json.dump(data['samples'], outfile)
        with open(os.path.join(donorPath, 'probes.txt'), 'w') as outfile:
            json.dump(data['probes'], outfile)

        if self.verbose:
            logging.getLogger(__name__).info('For {} samples_length: {}  probes_length: {} zscores_shape: {} '.format(donorId,len(data['samples']),len(data['probes']), zscores.shape))
        return {'samples' : data['samples'], 'probes' : data['probes'], 'zscores' : zscores}

    def filter_coordinates_and_zscores(self, img, index_to_allen_brain_api_dict, specimen, index):
        """
        Populate self.filtered_coords_and_zscores with zscores and coords for samples which belong to a particular specimen and spatially represented in the given voi.
        """
        revised_allen_brain_api_dict = dict.fromkeys(['zscores', 'coords', 'specimen', 'name'])
        revised_allen_brain_api_dict['name'] = 'img{}'.format(str(index+1))
        img_arr = img.get_data()
        invimgMni = inv(img.affine)
        T = np.dot(invimgMni, specimen['alignment3d'])
        coords = transform_samples_MRI_to_MNI52(index_to_allen_brain_api_dict['samples'], T)
        coords = (np.rint(coords)).astype(int)
        #How to use numpy.where
        coords = [np.array([-1, -1, -1]) if (coord > 0).sum() != 3 or img_arr[coord[0],coord[1],coord[2]] <= self.filter_threshold or img_arr[coord[0],coord[1],coord[2]] == 0 else coord for coord in coords]
        revised_allen_brain_api_dict['coords'] = [coord for coord in coords if (coord > 0).sum() == 3]
        revised_allen_brain_api_dict['zscores'] = [zscore for (coord, zscore) in zip(coords, index_to_allen_brain_api_dict['zscores']) if (coord > 0).sum() == 3]
        revised_allen_brain_api_dict['specimen'] = specimen['name']
        return revised_allen_brain_api_dict


    def query_allen_brain_api_partial(self, donorId):
        """
        Query Allen Brain Api for the given set of genes if they have not already been downloaded to gene_cache and save them on disk
        """
        base_url_query_allen_brain_api_partial = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        probes = ''.join('{},'.format(probe) for probe in self.probe_ids)[:-1]
        end_url_query_allen_brain_api_partial = "][donors$eq{}]".format(donorId)
        url = '{}{}{}'.format(base_url_query_allen_brain_api_partial, probes, end_url_query_allen_brain_api_partial)
        try:
            text = requests.get(url).json()
        except requests.exceptions.RequestException as e:
            logging.getLogger(__name__).error(e)
            raise
        data = text['msg']
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        donorpath = os.path.join(self.cache_dir, donorId)
        if not os.path.exists(donorpath):
            os.makedirs(donorpath)
        probes = data['probes']
        zscores = np.array([[float(data['probes'][i]['z-score'][j]) for i in range(len(data['probes']))] for j in range(len(data['samples']))])
        #READ PROBES
        with open(os.path.join(donorpath, 'probes.txt'), 'r') as f:
            probesC = json.load(f)
        #READ ZSCORES
        with open(os.path.join(donorpath, 'zscores.txt'), 'r') as f:
            zscoresC = np.loadtxt(f)
        #WRITE SAMPLES
        with open(os.path.join(donorpath, 'samples.txt'), 'w') as outfile:
            json.dump(data['samples'], outfile)
        #READ AND WRITE PROBES
        probes = probesC + probes       
        with open(os.path.join(donorpath, 'probes.txt'), 'w') as outfile:
            json.dump(probes, outfile)
        #WRITE ZSCORES
        zscores = np.append(zscoresC, zscores, axis=1)      
        np.savetxt(os.path.join(donorpath, 'zscores.txt'), zscores)

    def download_specimens(self):
        """
        Download names and transformation matrix for each specimen/donor from Allen Brain Api and save them on disk as specimenName.txt
        and specimenMat.txt respectively, load.
        """
        base_url_download_specimens = "http://api.brain-map.org/api/v2/data/Specimen/query.json?criteria=[name$eq"+"'"
        end_url_download_specimens = "']&include=alignment3d"
        specimens  = ['H0351.1015', 'H0351.1012', 'H0351.1016', 'H0351.2001', 'H0351.1009', 'H0351.2002']
        self.allen_brain_api_dict['specimen_info'] = []
        for specimen_id in specimens:
            url = '{}{}{}'.format(base_url_download_specimens, specimen_id, end_url_download_specimens)
            try:
                text = requests.get(url).json()
            except requests.exceptions.RequestException as e:
                logging.getLogger(__name__).info(e)
                raise
            self.allen_brain_api_dict['specimen_info'] = self.allen_brain_api_dict['specimen_info'] + [get_specimen_data(text['msg'][0])]
        if self.verbose:
            logging.getLogger(__name__).info('{}'.format(self.allen_brain_api_dict['specimen_info']))

        for donor, specimen in zip(self.donor_ids, self.allen_brain_api_dict['specimen_info']):
            with open(os.path.join(self.cache_dir, '{}/specimenName.txt'.format(donor)), 'w') as outfile:
                outfile.write(specimen['name'])
            np.savetxt(os.path.join(self.cache_dir, '{}/specimenMat.txt'.format(donor)), specimen['alignment3d'])


    def download_allen_brain_api_data(self):
        """
        Loop through the donors and call query_allen_brain_api() and populate allen_brain_api_dict
        """
        self.allen_brain_api_dict['samples_and_zscores'] = [self.query_allen_brain_api(donor) for donor in self.donor_ids]

    def download_and_retrieve_gene_data(self):
        """
        Download data from Allen Brain Api for the given set of genes and specimen
        """
        self.download_allen_brain_api_data()
        self.download_specimens()

    def set_candidate_genes(self, gene_list):
        """
        Set list of genes and prepare to read/download data for them.
        """
        self.gene_list = gene_list
        donorprobe = os.path.join(self.cache_dir, '{}/probes.txt'.format(self.donor_ids[0]))
        """
        If the cache doesnt exist then get all the probes associated with the genes and download and save the api and specimen information
        and populate allen_brain_api_dict['samples_and_zscores'] and allen_brain_api_data['specimen_info'].
        """
        if not os.path.exists(self.cache_dir):
            self.gene_list_to_download = self.gene_list[:]
            self.retrieve_probe_ids()
            self.download_and_retrieve_gene_data()
        else:
            """
            If the cache exists there are two possibilities.
            a) All the requested genes are present in the cache. In that case, read the downloaded api and specimen data from disk and
            populate allen_brain_api_dict['samples_and_zscores'] and allen_brain_api_dict['specimen_info']
            b) A few of the requested genes are present in the cache. In that case, read the downloaded api and specimen data from disk and
            populate allen_brain_api_dict['samples_and_zscores'] and allen_brain_api_dict['specimen_info'] for genes already present in the cache. For the rest, get all the probes
            associated with the genes, download and save the api and specimen information and populate allen_brain_api_dict['samples_and_zscores'] and
            allen_brain_api_dict['specimen_info'].
            """
            self.gene_list_to_download = [key for key in self.gene_list if key not in self.gene_cache.keys()]
            if self.gene_list_to_download:
                print('Microarray expression values of',len(self.gene_list_to_download),'gene(s) need(s) to be downloaded')
            if self.verbose:
                logging.getLogger(__name__).info('genes to be downloaded:{} '.format(self.gene_list_to_download))
            self.retrieve_probe_ids()
            if self.gene_list_to_download:
                for donor in self.donor_ids:
                    self.query_allen_brain_api_partial(donor)
            self.read_cached_api_specimen_data()


    def get_mean_zscores(self, combined_zscores):
        """
        Compute Winsorzed mean of zscores over all probes associated with a given gene. combined_zscores have zscores for all the probes and all the valid coordinates.
        As a gene_id_and_pvalues you get a numpy array of size len(self.filtered_coords_and_zscores)xlen(self.gene_list). self.genesymbol_and_mean_zscores['combined_zscores'][i][j] returns the winsorzed mean of
        jth gene taken over all the probes corresponding to the ith valid sample.
        """
        unique_gene_symbols = np.unique(self.gene_symbols)
        """
        A = [a,a,a,b,b,b,c,c]
        B = [a,b,c]
        Following line of code will give  indices = [[0,1,2],[3,4,5],[6,7]]
        """
        indices = [np.where(np.in1d(self.gene_symbols, genesymbol))[0] for genesymbol in unique_gene_symbols]
        """
        for i in range (len(unique_gene_symbols)):
                for j in range(len(combined_zscores)):
                    for k in range(len(indices[i])):
                        tmp[j] = combined_zscores[j][indices[i][k]][:]
                winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1))
        """
        winsorzed_mean_zscores =  np.array([[np.mean(sp.stats.mstats.winsorize([combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))], limits=0.1)) for i in range (len(unique_gene_symbols))] for j in range(len(combined_zscores))])
        self.genesymbol_and_mean_zscores['uniqueId'] = unique_gene_symbols
        self.genesymbol_and_mean_zscores['combined_zscores'] = winsorzed_mean_zscores

    def initialize_anova_factors(self):
        """
        Prepare self.anova_factors. Populate Age, Race, Area, Specimen, Zcores keys of self.anova_factors
        """
        combined_zscores = [roi_coord_and_zscore['zscores'][i] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
        #Populates self.specimenFactors (id, race, gender, name, age)
        self.read_specimen_factors(self.cache_dir)
        if self.verbose:
            logging.getLogger(__name__).info("number of specimens: {} name: {}".format(len(self.specimen_factors), len(self.specimen_factors['name'])))
        #Populates self.genesymbol_and_mean_zscores (uniqueid and zscores)
        self.get_mean_zscores(combined_zscores)

        self.n_genes = len(self.genesymbol_and_mean_zscores['combined_zscores'][0])
        self.anova_factors['Area'] = [roi_coord_and_zscore['name'] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
        self.anova_factors['Specimen'] = [roi_coord_and_zscore['specimen'] for roi_coord_and_zscore in self.filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
        '''
        Both Age and Race should have len(self.filtered_coords_and_zscores) entries. The following three lines are used to get the correct values from specimenFactors['Age'] and specimenFactors['Race'] using specimenFactors['name'] and repeat them the required number of times as given by self.anova_factors['Specimen']
        '''
        st = set(self.specimen_factors['name'])
        self.anova_factors['Age'] = [self.specimen_factors['age'][self.specimen_factors['name'].index(specimen_name)] for ind, specimen_name in enumerate(self.anova_factors['Specimen'])]
        self.anova_factors['Race'] = [self.specimen_factors['race'][self.specimen_factors['name'].index(specimen_name)] for ind, specimen_name in enumerate(self.anova_factors['Specimen'])]

    def first_iteration(self):
        """
        Perform one iteration of ANOVA. Use output of this to populate F_vec_ref_anovan which becomes initial estimate of n_rep passes of FWE.
        """
        self.F_vec_ref_anovan = np.zeros(self.n_genes)
        for i in range(self.n_genes):
            self.anova_factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,i]
            mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.anova_factors).fit()
            aov_table = sm.stats.anova_lm(mod, typ=1)
            if self.verbose:
                logging.getLogger(__name__).info('aov table: {}'.format(aov_table))
            #F_vec_ref_anovan is used as an initial condition to F_mat_perm_anovan in fwe_correction
            self.F_vec_ref_anovan[i] = aov_table['F'][0]

    def do_anova_with_permutation_gene(self, index_to_uniqueId):
        self.anova_factors['Area'] = np.random.permutation(self.anova_factors['Area'])
        self.anova_factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,index_to_uniqueId]
        mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.anova_factors).fit()
        aov_table = sm.stats.anova_lm(mod, typ=1)
        return aov_table['F'][0]

    #Pool inside a  pool, is it a good idea?
    def do_anova_with_permutation_rep(self, rep):
        self.F_vec_perm_anovan = list(map(self.do_anova_with_permutation_gene, range(0,self.n_genes)))
        return self.F_vec_perm_anovan

    def fwe_correction(self):
        """
        Perform n_rep passes of FWE using gene_id_and_pvalues of first_iteration() as an initial guess
        """
        invn_rep = 1/self.n_rep
        initial_guess_F_vec = self.F_vec_ref_anovan
        pool = multiprocessing.Pool()
        self.F_mat_perm_anovan = np.array(list(pool.map(self.do_anova_with_permutation_rep, range(1,self.n_rep))))
        self.F_mat_perm_anovan = np.insert(self.F_mat_perm_anovan, 0, initial_guess_F_vec, axis=0)
        self.accumulate_gene_id_and_pvalues()

    def accumulate_gene_id_and_pvalues(self):
        """
        Populate pvalues and geneids for the gene_id_and_pvalues dict after n_rep passes of FWE
        """
        invn_rep = 1/self.n_rep
        #ref represenets maximum p value for each gene across n_rep repetitions
        ref = self.F_mat_perm_anovan.max(1)
        #compute family wise error corrected p value
        self.FWE_corrected_p =  [len([1 for a in ref if a >= f])/self.n_rep if sys.version_info[0] >= 3 else len([1 for a in ref if a >= f])*invn_rep for f in self.F_vec_ref_anovan]
        self.gene_id_and_pvalues = dict(zip(self.genesymbol_and_mean_zscores['uniqueId'], self.FWE_corrected_p))
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
        Download various factors such as age, name, race, gender of the six specimens from Allen Brain Api and create a dict.
        """
        url_build_specimen_factors = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"
        try:
            text = requests.get(url_build_specimen_factors).json()
        except requests.exceptions.RequestException as e:
            print('In build_specimen_factors')
            print(e)
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
        Read various factors such as age, name, race, gender of the six specimens from disk.
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
