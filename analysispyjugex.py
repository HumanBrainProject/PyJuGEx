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
def get_specimen_data(info):
    """
    For each specimen, extract the name and alignment matrix and put into a dict object
    """
    specimenD = dict()
    specimenD['name'] = info['name']
    x = info['alignment3d']
    specimenD['alignment3d'] = np.array([
    [x['tvr_00'], x['tvr_01'], x['tvr_02'], x['tvr_09']],
    [x['tvr_03'], x['tvr_04'], x['tvr_05'], x['tvr_10']],
    [x['tvr_06'], x['tvr_07'], x['tvr_08'], x['tvr_11']],
    [0, 0, 0, 1]])
    return specimenD

def transformSamples(samples, T):
    """
    Convert the MRI coordinates of samples to MNI152 space
    """
    np_T = np.array(T[0:3, 0:4])
    mri = np.vstack(s["sample"]["mri"] for s in samples)
    add = np.ones((len(mri), 1), dtype=np.int)
    mri = np.append(mri, add, axis=1)
    mri = np.transpose(mri)
    coords = np.matmul(np_T, mri)
    coords = coords.transpose()
    return coords

class Analysis:

    def __init__(self, gene_cache, verbose=False):
        """
        initialize_anova_data the Analysis class with various internal variables -

        probeids = list of probe ids associated with the give list of genes which are not present in genecache, if it exists.
        genelist = given list of genes
        downloadgenelist = list of genes whose information is not in the cache yet, needs to be downloaded
        genesymbols = same length as probeids, each probe is represented with the corresponding genesymbol, used by getmeanzscores()
        genecache = dictionary to indicate which keys are present in the cache
        donorids = six donor ids of Allen Brain API
        apidata = dictionary to store Allen Brain API data, with two keys - apiinfo contains zscores, samples etc and specimeninfo contains
        age, race, sex of the six donors represented by donorids
        vois = list of two nii volumes for each region of interest
        main_r = Internal variable for storing mni coordinates and zscores corresponsing to each region of interest.
        mapthreshold = Internal variable used at expressionSpmCorrelation() to select or reject a sample
        n_rep = number of iterations of FWE correction
        cache = Disk location where data from Allen Brain API has been downloaded and stored.
        anova_data = Internal variable used by the ANOVA module - contains five keys - 'Age', 'Race', 'Specimen', 'Area', 'Zscores'
        all_probe_data = dictionary with two keys - uniqueid, combinedzscores where each gene and the winsorzed mean zscores over all probes associated with that gene is stored
        result = dict for storing gene ids and associated p values.
        """
        self.probeids = []
        self.genelist = []
        self.downloadgenelist = []
        self.genesymbols = []
        self.genecache = {}
        self.donorids = ['15496', '14380', '15697', '9861', '12876', '10021'] #HARDCODING DONORIDS
        self.apidata = dict.fromkeys(['apiinfo', 'specimeninfo'])        
        self.specimenFactors = dict.fromkeys(['id', 'name', 'race', 'gender', 'age'])
        self.apidata['specimenInfo'] = []
        self.apidata['apiinfo'] = []
        self.vois = []
        self.main_r = []
        self.mapthreshold = 0.2
        self.n_rep = 1000
        self.cachedir = gene_cache
        self.verboseflag = verbose
        self.anova_data = dict.fromkeys(['Age', 'Race', 'Specimen', 'Area', 'Zscores'])
        self.all_probe_data = dict.fromkeys(['uniqueId', 'combined_zscores'])
        '''
        Removes any folder that has not been written to properly, most likely due to force quit
        '''
        probepath = os.path.join(self.cachedir, self.donorids[0]+'/probes.txt')
        if os.path.exists(self.cachedir) and not os.path.exists(probepath):
            shutil.rmtree(self.cachedir, ignore_errors = False)
        '''
        Creates a gene cache to indicate which genes are present in the cache
        '''
        if not os.path.exists(self.cachedir):
            print(self.cachedir,' does not exist. It will take some time ')
        else:
            self.creategenecache()
            print(len(self.genecache),' genes exist in ', self.cachedir)

    def DifferentialAnalysis(self, genelist, roi1, roi2):
        """
        Driver routine
        """
        if not genelist:
            raise ValueError('Atleast one gene is needed for the analysis')
        if not (isinstance(roi1, nib.nifti1.Nifti1Image) or isinstance(roi2, nib.nifti1.Nifti1Image)):
            raise ValueError('Atleast two valid regions of interest are needed')
        self.set_candidate_genes(genelist)
        self.set_ROI_MNI152(roi1, 0)
        self.set_ROI_MNI152(roi2, 1)
        print('Starting the analysis. This may take some time.....')
        self.perform_anova()
        return self.result

    def creategenecache(self):
        """
        Create a dictionary with an entry for each gene whose api information has been downloaded.
        We scan self.cachedir/15496/probes.txt for the gene symbols
        """
        with open(os.path.join(self.cachedir, self.donorids[0]+'/probes.txt'), 'r') as f:
            probes = json.load(f)
        for p in probes:
            self.genecache.update({p['gene-symbol'] : None})
        if self.verboseflag:
            print(self.genecache)

    def retrieveprobeids(self):
        """
        Retrieve probe ids for the given gene lists, update self.probeids which will be used by queryapi() or queryapipartial() to
        form the url and update self.genesymbols to be used by getmeanzscores()
        """
        base_retrieve_probeids = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq"
        end_retrieve_probeids = "],rma::options[only$eq'probes.id']"
        if self.verboseflag:
            print('genelist ',self.genelist)
        for gene in self.genelist:
            url = '{base}{gene_symbol}{end}'.format(base = base_retrieve_probeids, gene_symbol = gene, end = end_retrieve_probeids)
            if self.verboseflag:
                print(url)
            try:
                response = requests.get(url)
            except requests.exceptions.RequestException as e:
                print('In retreiveprobeids')
                print(e)
                raise
            data = xmltodict.parse(response.text)
            self.probeids = self.probeids + [d['id'] for d in data['Response']['probes']['probe'] if gene in self.downloadgenelist]
            self.genesymbols = self.genesymbols + [gene for d in data['Response']['probes']['probe']]

        if self.verboseflag:
            print('probeids: ',self.probeids)
            print('genesymbols: ',self.genesymbols)


    def readCachedApiSpecimenData(self):
        """
        Read cached Allen Brain Api data from disk location and update self.apidata['specimenInfo'] and self.apidata['apiInfo']
        """
        for d in self.donorids:
            donorpath = os.path.join(self.cachedir, d)
            specimen = dict.fromkeys(['name', 'alignment3d'])
            with open(os.path.join(donorpath, 'specimenMat.txt'), 'r') as f:
                specimen['alignment3d'] = np.loadtxt(f)
            with open(os.path.join(donorpath, 'specimenName.txt'), 'r') as f:
                specimen['name'] = f.read()
            self.apidata['specimenInfo'] = self.apidata['specimenInfo'] + [specimen]

            #LOAD SAMPLES
            with open(os.path.join(donorpath, 'samples.txt'), 'r') as f:
                samplesC = json.load(f)
            #LOAD PROBES
            with open(os.path.join(donorpath, 'probes.txt'), 'r') as f:
                probesC = json.load(f)
            #LOAD ZSCORES
            with open(os.path.join(donorpath, 'zscores.txt'), 'r') as f:
                zscoresC = np.loadtxt(f)

            self.apidata['apiinfo'] = self.apidata['apiinfo']  + [{'samples' : samplesC, 'zscores' : zscoresC}]
            if self.verboseflag:
                print('inside readcachedata ',len(self.apidata['apiinfo'][-1]['samples']), ' ', self.apidata['apiinfo'][-1]['zscores'].shape)

    def set_ROI_MNI152(self, voi, index):
        """
        For each specimen and for the given voi populate self.main_r with zscores and valid coordinates
        """
        if index < 0 or index > 1:
            raise ValueError('only 0 and 1 are valid choices')
        for i in range(0, len(self.apidata['specimenInfo'])):
            self.main_r.append(self.expressionSpmCorrelation(voi, self.apidata['apiinfo'][i], self.apidata['specimenInfo'][i], index))

    def queryapi(self, donorId):
        """
        Query Allen Brain Api for given set of genes
        """
        base_query_api = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        probes = ''.join(p+"," for p in self.probeids)[:-1]
        end_query_api = "][donors$eq"+donorId+"]"
        url = '{base}{probes}{end}'.format(base = base_query_api, probes = probes, end = end_query_api)
        '''
        url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in" + ''.join(p+"," for p in self.probeids)
        url = url[:-1]
        url += "][donors$eq"+donorId+"]"
        '''
        try:
            response = requests.get(url)
            text = requests.get(url).json()
        except requests.exceptions.RequestException as e:
            print('In queryapi ')
            print(e)
            raise
        data = text['msg']
        if not os.path.exists(self.cachedir):
            os.makedirs(self.cachedir)
        donorPath = os.path.join(self.cachedir, donorId)
        if not os.path.exists(donorPath):
            os.makedirs(donorPath)

        zscores = np.array([[float(data['probes'][i]['z-score'][j]) for i in range(len(data['probes']))] for j in range(len(data['samples']))])

        with open(os.path.join(donorPath, 'zscores.txt'), 'wb') as f:
            np.savetxt(f, zscores, fmt = '%.5f')
        with open(os.path.join(donorPath, 'samples.txt'), 'w') as outfile:
            json.dump(data['samples'], outfile)
        with open(os.path.join(donorPath, 'probes.txt'), 'w') as outfile:
            json.dump(data['probes'], outfile)

        if self.verboseflag:
            print('For ',donorId,' samples_length: ',len(data['samples']),' probes_length: ',len(data['probes']),' zscores_shape: ',zscores.shape)
        return {'samples' : data['samples'], 'probes' : data['probes'], 'zscores' : zscores}

    def expressionSpmCorrelation(self, img, apidataind, specimen, index):
        """
        Populate self.main_r with zscores and coords for samples which belong to a particular specimen and spatially represented in the given voi.
        """
        revisedApiData = dict.fromkeys(['zscores', 'coords', 'specimen', 'name'])
        revisedApiData['name'] = 'img'+str(index+1)
        dataImg = img.get_data()
        invimgMni = inv(img.affine)
        T = np.dot(invimgMni, specimen['alignment3d'])
        coords = transformSamples(apidataind['samples'], T)
        coords = (np.rint(coords)).astype(int)
        #How to use numpy.where
        coords = [np.array([-1, -1, -1]) if (coord > 0).sum() != 3 or dataImg[coord[0],coord[1],coord[2]] <= self.mapthreshold or dataImg[coord[0],coord[1],coord[2]] == 0 else coord for coord in coords]
        revisedApiData['coords'] = [c for c in coords if (c > 0).sum() == 3]
        revisedApiData['zscores'] = [z for (c, z) in zip(coords, apidataind['zscores']) if (c > 0).sum() == 3]
        revisedApiData['specimen'] = specimen['name']
        return revisedApiData


    def queryapipartial(self, donorId):
        """
        Query Allen Brain Api for the given set of genes if they have not already been downloaded to genecache and save them on disk
        """

        base_url_queryapipartial = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        probes = ''.join(p+"," for p in self.probeids)[:-1]
        end_url_queryapipartial = "][donors$eq"+donorId+"]"
        url = '{base}{probes}{end}'.format(base = base_url_queryapipartial, probes = probes, end = end_url_queryapipartial)

        if self.verboseflag:
            print(url)
        try:
            text = requests.get(url).json()
        except requests.exceptions.RequestException as e:
            print('In queryapipartial ')
            print(e)
            raise
        data = text['msg']
        if not os.path.exists(self.cachedir):
            os.makedirs(self.cachedir)
        donorpath = os.path.join(self.cachedir, donorId)
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

    def downloadspecimens(self):
        """
        Download names and transformation matrix for each specimen/donor from Allen Brain Api and save them on disk as specimenName.txt
        and specimenMat.txt respectively, load.
        """
        base_url_downloadspecimens = "http://api.brain-map.org/api/v2/data/Specimen/query.json?criteria=[name$eq"+"'"
        end_url_downloadspecimens = "']&include=alignment3d"
        specimens  = ['H0351.1015', 'H0351.1012', 'H0351.1016', 'H0351.2001', 'H0351.1009', 'H0351.2002']
        self.apidata['specimenInfo'] = []
        for specimen_id in specimens:
            url = '{base}{specimen_id}{end}'.format(base = base_url_downloadspecimens, specimen_id = specimen_id, end = end_url_downloadspecimens)
            if self.verboseflag:
                print(url)
            try:
                text = requests.get(url).json()
            except requests.exceptions.RequestException as e:
                print('In downloadspecimens ')
                print(e)
                raise
            self.apidata['specimenInfo'] = self.apidata['specimenInfo'] + [get_specimen_data(text['msg'][0])]
        if self.verboseflag:
            print(self.apidata['specimenInfo'])

        for d, s in zip(self.donorids, self.apidata['specimenInfo']):
            with open(os.path.join(self.cachedir, d+'/specimenName.txt'), 'w') as outfile:
                outfile.write(s['name'])
            np.savetxt(os.path.join(self.cachedir, d+'/specimenMat.txt'), s['alignment3d'])


    def getapidata(self):
        """
        Loop through the donors and call queryapi() and populate apidata
        """
        self.apidata['apiinfo'] = [self.queryapi(d) for d in self.donorids]

    def download_and_retrieve_gene_data(self):
        """
        Download data from Allen Brain Api for the given set of genes and specimen
        """
        self.getapidata()
        self.downloadspecimens()

    def set_candidate_genes(self, genelist):
        """
        Set list of genes and prepare to read/download data for them.
        """
        self.genelist = genelist
        donorprobe = os.path.join(self.cachedir, self.donorids[0]+'/probes.txt')
        """
        If the cache doesnt exist then get all the probes associated with the genes and download and save the api and specimen information
        and populate apidata['apiinfo'] and apidata['specimeninfo'].
        """
        if not os.path.exists(self.cachedir):
            self.downloadgenelist = self.genelist[:]
            if self.verboseflag:
                print('self.cachedir does not exist')
            self.retrieveprobeids()
            self.download_and_retrieve_gene_data()
        else:
            """
            If the cache exists there are two possibilities.
            a) All the requested genes are present in the cache. In that case, read the downloaded api and specimen data from disk and
            populate apidata['apiinfo'] and apidata['specimeninfo']
            b) A few of the requested genes are present in the cache. In that case, read the downloaded api and specimen data from disk and
            populate apidata['apiinfo'] and apidata['specimeninfo'] for genes already present in the cache. For the rest, get all the probes
            associated with the genes, download and save the api and specimen information and populate apidata['apiinfo'] and
            apidata['specimeninfo'].
            """
            self.downloadgenelist = [k for k in self.genelist if k not in self.genecache.keys()]
            if self.downloadgenelist:
                print('Microarray expression values of',len(self.downloadgenelist),'gene(s) need(s) to be downloaded')
            if self.verboseflag:
                print(self.downloadgenelist)
            self.retrieveprobeids()
            if self.downloadgenelist:
                for d in self.donorids:
                    self.queryapipartial(d)
            self.readCachedApiSpecimenData()


    def getmeanzscores(self, combined_zscores):
        """
        Compute Winsorzed mean of zscores over all probes associated with a given gene. combined_zscores have zscores for all the probes and all the valid coordinates.
        As a result you get a numpy array of size len(self.main_r)xlen(self.genelist). self.all_probe_data['combined_zscores'][i][j] returns the winsorzed mean of
        jth gene taken over all the probes corresponding to the ith valid sample.
        """
        unique_gene_symbols = np.unique(self.genesymbols)
        """
        A = [a,a,a,b,b,b,c,c]
        B = [a,b,c]
        Following line of code will give  indices = [[0,1,2],[3,4,5],[6,7]]
        """
        indices = [np.where(np.in1d(self.genesymbols, x))[0] for x in unique_gene_symbols]
        """
        for i in range (len(unique_gene_symbols)):
                for j in range(len(combined_zscores)):
                    for k in range(len(indices[i])):
                        tmp[j] = combined_zscores[j][indices[i][k]][:]
                winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1))
        """
        winsorzed_mean_zscores =  np.array([[np.mean(sp.stats.mstats.winsorize([combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))], limits=0.1)) for i in range (len(unique_gene_symbols))] for j in range(len(combined_zscores))])
        self.all_probe_data['uniqueId'] = unique_gene_symbols
        self.all_probe_data['combined_zscores'] = winsorzed_mean_zscores

    def initialize_anova_data(self):
        """
        Prepare self.anova_data. Populate Age, Race, Area, Specimen, Zcores keys of self.anova_data
        """
        combined_zscores = [r['zscores'][i] for r in self.main_r for i in range(len(r['zscores']))]
        #Populates self.specimenFactors (id, race, gender, name, age)
        self.read_specimen_factors(self.cachedir)
        if self.verboseflag:
            print("number of specimens ", len(self.specimenFactors), " name: ", len(self.specimenFactors['name']))  
        #Populates self.all_probe_data (uniqueid and zscores)
        self.getmeanzscores(combined_zscores)

        self.n_genes = len(self.all_probe_data['combined_zscores'][0])
        self.anova_data['Area'] = [r['name'] for r in self.main_r for i in range(len(r['zscores']))]
        self.anova_data['Specimen'] = [r['specimen'] for r in self.main_r for i in range(len(r['zscores']))]
        #Both Age and Race should have len(self.main_r) entries. The following three lines are used to get the correct values from specimenFactors['Age'] and specimenFactors['Race'] using specimenFactors['name'] and repeat them the required number of times as given by self.anova_data['Specimen']
        st = set(self.specimenFactors['name'])
        self.anova_data['Age'] = [self.specimenFactors['age'][self.specimenFactors['name'].index(a)] for ind, a in enumerate(self.anova_data['Specimen'])]
        self.anova_data['Race'] = [self.specimenFactors['race'][self.specimenFactors['name'].index(a)] for ind, a in enumerate(self.anova_data['Specimen'])]
        if self.verboseflag:
            print('race')
            print(self.anova_data['Race'])
            print('age')
            print(self.anova_data['Age'])
            print(len(self.genelist))

    def first_iteration(self):
        """
        Perform one iteration of ANOVA. Use output of this to populate F_vec_ref_anovan which becomes initial estimate of n_rep passes of FWE.
        """
        self.F_vec_ref_anovan = np.zeros(self.n_genes)
        for i in range(0, self.n_genes):
            self.anova_data['Zscores'] = self.all_probe_data['combined_zscores'][:,i]
            mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.anova_data).fit()
            aov_table = sm.stats.anova_lm(mod, typ=1)
            if self.verboseflag:
                print(aov_table)
            #F_vec_ref_anovan is used as an initial condition to F_mat_perm_anovan in fwe_correction
            self.F_vec_ref_anovan[i] = aov_table['F'][0]

    def do_anova_with_permutation_gene(self, index_to_uniqueId):
        self.anova_data['Area'] = np.random.permutation(self.anova_data['Area'])
        self.anova_data['Zscores'] = self.all_probe_data['combined_zscores'][:,index_to_uniqueId]
        mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.anova_data).fit()
        aov_table = sm.stats.anova_lm(mod, typ=1)
        return aov_table['F'][0]

    #Pool inside a  pool, is it a good idea?
    def do_anova_with_permutation_rep(self, rep):
        self.F_vec_perm_anovan = list(map(self.do_anova_with_permutation_gene, range(0,self.n_genes)))
        return self.F_vec_perm_anovan

    def fwe_correction(self):
        """
        Perform n_rep passes of FWE using result of first_iteration() as an initial guess
        """
        invn_rep = 1/self.n_rep
        initial_guess_F_vec = self.F_vec_ref_anovan
        pool = multiprocessing.Pool()
        self.F_mat_perm_anovan = np.array(list(pool.map(self.do_anova_with_permutation_rep, range(1,self.n_rep))))
        self.F_mat_perm_anovan = np.insert(self.F_mat_perm_anovan, 0, initial_guess_F_vec, axis=0)
        self.accumulate_result()

    def accumulate_result(self):
        """
        Populate pvalues and geneids for the result dict after n_rep passes of FWE
        """
        invn_rep = 1/self.n_rep
        #ref represenets maximum p value for each gene across n_rep repetitions
        ref = self.F_mat_perm_anovan.max(1)
        #compute family wise error corrected p value
        self.FWE_corrected_p =  [len([1 for a in ref if a >= f])/self.n_rep if sys.version_info[0] >= 3 else len([1 for a in ref if a >= f])*invn_rep for f in self.F_vec_ref_anovan]
        self.result = dict(zip(self.all_probe_data['uniqueId'], self.FWE_corrected_p))
        if self.verboseflag:
            print(self.result)

    def perform_anova(self):
        """
        Perform one way anova on zscores as the dependent variable and specimen factors such as age, race, name and area
        as independent variables
        """
        self.initialize_anova_data()
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

        self.specimenFactors['id'] = [r['id'] for r in res]
        self.specimenFactors['name'] = [r['name'] for r in res]
        self.specimenFactors['race'] = [r['race_only'] for r in res]
        self.specimenFactors['gender'] = [r['sex'] for r in res]
        self.specimenFactors['age'] = [r['age']['days']/365 for r in res]


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
        self.specimenFactors = dict()
        self.specimenFactors['id'] = [r['id'] for r in res]
        self.specimenFactors['name'] = [r['name'] for r in res]
        self.specimenFactors['race'] = [r['race_only'] for r in res]
        self.specimenFactors['gender'] = [r['sex'] for r in res]
        self.specimenFactors['age'] = [r['age']['days']/365 for r in res]
