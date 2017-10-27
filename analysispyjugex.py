# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
import os
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
if sys.version_info[0] >=3 :
    import urllib

def getmeanzscores(gene_symbols, combined_zscores, area1len, area2len):
    """
    Compute Winsorzed mean of zscores over all genes.
    """
    unique_gene_symbols = np.unique(gene_symbols)
    indices = [np.where(np.in1d(gene_symbols, x))[0] for x in unique_gene_symbols]
    winsorzed_mean_zscores = np.zeros((len(combined_zscores), len(unique_gene_symbols)))
    tmp = [[] for _ in range(len(combined_zscores))]
    for i in range (0, len(unique_gene_symbols)):
        for j in range(0, len(combined_zscores)):
            tmp[j] = [combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))]
            winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1)) #maybe add a special case for one element in tmp
    res = dict.fromkeys(['uniqueId', 'combined_zscores', 'area1_zscores', 'area2_zscores'])
    res['uniqueId'] = unique_gene_symbols
    res['combined_zscores'] = winsorzed_mean_zscores
    res['area1_zscores'] = winsorzed_mean_zscores[0:area1len]
    res['area2_zscores'] = winsorzed_mean_zscores[area1len:]
    return res

def getSpecimenData(info):
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

def buildSpecimenFactors(cache):
    """
    Download various factors such as age, name, race, gender of the six specimens from Allen Brain Api
    and create a dict.
    """
    url = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"
    specimenFactors = dict()
    specimenFactors['id'] = []
    specimenFactors['name'] = []
    specimenFactors['race'] = []
    specimenFactors['gender'] = []
    specimenFactors['age'] = []
    if sys.version_info[0] >= 3:
        response = urllib.request.urlopen(url).read().decode('utf8')
        text = json.loads(response)
    else:
        #response = requests.get(url)
        #text = json.loads(response.text)
        text = requests.get(url).json()
    factorPath = os.path.join(cache, 'specimenFactors.txt')
    with open(factorPath, 'w') as outfile:
        json.dump(text, outfile)
    res = text['msg']
    for i in range(0, len(res)):
        specimenFactors['id'].append(res[i]['id'])
        specimenFactors['name'].append(res[i]['name'])
        specimenFactors['race'].append(res[i]['race_only'])
        specimenFactors['gender'].append(res[i]['sex'])
        specimenFactors['age'].append(res[i]['age']['days']/365)
    return specimenFactors

def readSpecimenFactors(cache):
    """
    Read various factors such as age, name, race, gender of the six specimens from disk.
    """
    specimenFactors = dict.fromkeys(['id', 'name', 'race', 'gender', 'age'])
    specimenFactors['id'] = []
    specimenFactors['name'] = []
    specimenFactors['race'] = []
    specimenFactors['gender'] = []
    specimenFactors['age'] = []
    fileName = os.path.join(cache, 'specimenFactors.txt')
    if not os.path.exists(fileName):
        specimenFactors = buildSpecimenFactors(cache)
    f = open(fileName, "r")
    content = json.load(f)
    f.close()
    res = content['msg']
    for i in range(0, len(res)):
        specimenFactors['id'].append(res[i]['id'])
        specimenFactors['name'].append(res[i]['name'])
        specimenFactors['race'].append(res[i]['race_only'])
        specimenFactors['gender'].append(res[i]['sex'])
        specimenFactors['age'].append(res[i]['age']['days']/365)
    return specimenFactors

class Analysis:

    def __init__(self, gene_cache, verbose=False):
        """
        Initialize the Analysis class with various internal variables -
        refreshcache = True if no cache location is given, false otherwise.
        cache = Disk location where data from Allen Brain API has been downloaded and stored.
        probeids = list of probe ids associated with the give list of genes.
        genelist = given list of genes
        genesymbols =
        donorids = size donor ids of Allen Brain API
        vois = list fo two nii volumes for each region of interest
        main_r = Internal variable for storing mni coordinates and zscores corresponsing to each region of interest.
        mapthreshold = Internal variable to select or reject a sample
        result = dict for storing gene ids and associated p values.
        """
        self.probeids = []
        self.genelist = []
        self.downloadgenelist = []
        self.genesymbols = []
        self.genecache = {}
        self.donorids = ['15496', '14380', '15697', '9861', '12876', '10021'] #HARDCODING DONORIDS
        self.apidata = dict.fromkeys(['apiinfo', 'specimeninfo'])
        self.apidata['specimenInfo'] = []
        self.apidata['apiinfo'] = []
        self.vois = []
        self.main_r = []
        self.mapthreshold = 0.2
        self.result = None
        self.cache = gene_cache
        self.verboseflag = verbose
        if not os.path.exists(gene_cache):
            print(gene_cache,' does not exist. It will take some time ')
        else:
            self.creategenecache()
            print(len(self.genecache),' genes exist in ', self.cache)

    def DifferentialAnalysis(self, genelist, roi1, roi2):
        self.set_candidate_genes(genelist)
        self.set_ROI_MNI152(roi1, 0)
        self.set_ROI_MNI152(roi2, 1)

        print('Starting the analysis. This may take some time.....')
        self.run()
        result = self.pvalues()
        return result

    def creategenecache(self):
        donorpath = os.path.join(self.cache, self.donorids[0])
        filename = os.path.join(donorpath, 'probes.txt')
        f = open(filename, "r")
        probes = json.load(f)
        for p in probes:
            self.genecache.update({p['gene-symbol'] : None})
        if(self.verboseflag):
            print(self.genecache)
    '''
    def readgenetoprobeidscache(self):
        f = open('genetoprobecache.txt', 'r')
        genecachedup = json.load(f)
        for g in self.genelist:
            if g in self.downloadgenelist:
                for c in range(0, len(genecachedup[g])):
                    self.probeids = self.probeids + [genecachedup[g][c]['id']]
            for c in range(0, len(genecachedup[g])):
                self.genesymbols = self.genesymbols + [g]
        if self.verboseflag:
            print(self.genesymbols)
            print(self.probeids)
    '''
    def retrieveprobeids(self):
        """
        Retrieve probe ids for the given gene lists
        """
        connection = False
        if(self.verboseflag):
            print('genelist ',self.genelist)
        for g in self.genelist:
            url = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq"
            url += g
            url += "],rma::options[only$eq'probes.id']"
            if(self.verboseflag):
                print(url)
            if sys.version_info[0] >= 3:
                try:
                    response = urllib.request.urlopen(url).read()
                except URLError as e:
                    if hasattr(e, 'reason'):
                        print('We failed to reach a server.')
                        print('Reason: ', e.reason)
                    elif hasattr(e, 'code'):
                        print('The server couldn\'t fulfill the request.')
                        print('Error code: ', e.code)
            else:
                try:
                    response = requests.get(url)
                except requests.exceptions.RequestException as e:
                    print(e)
                    connection = True
            '''
            if connection is True:
                self.readgenetoprobeidscache()
            else:
            '''
            if sys.version_info[0] < 3:
                data = xmltodict.parse(response.text)
            else:
                data = xmltodict.parse(response)
            for d in data['Response']['probes']['probe']:
                if g in self.downloadgenelist:
                    self.probeids = self.probeids + [d['id']]
                self.genesymbols = self.genesymbols + [g]
            if(self.verboseflag):
                print('probeids: ',self.probeids)
                print('genesymbols: ',self.genesymbols)


    def readCachedApiSpecimenData(self):
        """
        Read cached Allen Brain Api data from disk location
        """
        #self.downloadspecimens()
        for d in self.donorids:
            donorpath = os.path.join(self.cache, d)
            fileNameM = os.path.join(donorpath, 'specimenMat.txt')
            mat = np.loadtxt(fileNameM)
            fileNameN = os.path.join(donorpath, 'specimenName.txt')
            f = open(fileNameN, 'r')
            name = f.read()
            f.close()
            specimen = dict.fromkeys(['name', 'alignment3d'])
            specimen['name'] = name
            specimen['alignment3d'] = mat
            self.apidata['specimenInfo'].append(specimen)
            #LOAD SAMPLES
            fileName = os.path.join(donorpath, 'samples.txt')
            f = open(fileName, "r")
            samplesC = json.load(f)
            f.close()
            #LOAD PROBES
            fileName = os.path.join(donorpath, 'probes.txt')
            f = open(fileName, "r")
            probesC = json.load(f)
            f.close()
            #LOAD ZSCORES
            fileName = os.path.join(donorpath, 'zscores.txt')
            zscoresC = np.loadtxt(fileName)
            apiDataC = dict()
            apiDataC['samples'] = samplesC
            apiDataC['probes'] = probesC
            apiDataC['zscores'] = zscoresC
            self.apidata['apiinfo'].append(apiDataC)
            if(self.verboseflag):
                print('inside readcachedata ',len(apiDataC['samples']), ' ', apiDataC['zscores'].shape, ' ', len(apiDataC['probes']))

    def set_ROI_MNI152(self, voi, index):
        """
        Set the region of interest from the downloaded nii files
        """
        for i in range(0, len(self.apidata['specimenInfo'])):
            revisedApiDataCombo = dict()
            revisedApiData = self.expressionSpmCorrelation(voi, self.apidata['apiinfo'][i], self.apidata['specimenInfo'][i]) #Maybe an index will work instead of expressionspmcorrelation
            revisedApiDataCombo['zscores'] = revisedApiData['zscores'][:]
            revisedApiDataCombo['coords'] = revisedApiData['coords'][:]
            revisedApiDataCombo['samples'] = revisedApiData['samples'][:]
            revisedApiDataCombo['probes'] = revisedApiData['probes'][:]
            revisedApiDataCombo['specimen'] = revisedApiData['specimen']
            if index == 0:
                revisedApiDataCombo['name'] = 'img1'
            elif index == 1:
                revisedApiDataCombo['name'] = 'img2'
            else:
                print('only 0 and 1 are valid choices')
                exit()
            if(self.verboseflag):
                print('extractexplevel img1: ',revisedApiDataCombo['specimen'],' ',len(revisedApiDataCombo['coords']))
            self.main_r.append(revisedApiDataCombo)


    def queryapi(self, donorId):
        url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        for p in self.probeids:
            url += p
            url += ","
        url = url[:-1]
        url += "][donors$eq"
        url += donorId
        url += "]"
        if sys.version_info[0] >= 3:
            response = urllib.request.urlopen(url).read().decode('utf8')
            text = json.loads(response)
        else:
            response = requests.get(url)
            text = requests.get(url).json()
        data = text['msg']
        samples = []
        probes = []
        if not os.path.exists(self.cache):
            os.makedirs(self.cache)
        donorPath = os.path.join(self.cache, donorId)
        if not os.path.exists(donorPath):
            os.makedirs(donorPath)
        nsamples = len(data['samples'])
        nprobes = len(data['probes'])
        samples = data['samples']
        probes = data['probes']

        zscores = np.zeros((nsamples, nprobes))
        for i in range(0, nprobes):
            for j in range(0, nsamples):
                zscores[j][i] = probes[i]['z-score'][j]


        fileName = os.path.join(donorPath, 'samples.txt')
        with open(fileName, 'w') as outfile:
            json.dump(data['samples'], outfile)

        fileName = os.path.join(donorPath, 'probes.txt')
        with open(fileName, 'w') as outfile:
            json.dump(data['probes'], outfile)

        fileName = os.path.join(donorPath, 'zscores.txt')
        np.savetxt(fileName, zscores)
    
        apiData = dict()
        apiData['samples'] = samples
        apiData['probes'] = probes
        apiData['zscores'] = zscores
        if(self.verboseflag):
            print('For ',donorId,' samples_length: ',len(apiData['samples']),' probes_length: ',len(apiData['probes']),' zscores_shape: ',apiData['zscores'].shape)
        return apiData



    def expressionSpmCorrelation(self, img, apidataind, specimen):
        """
        Create internal data structures with valid coordinates in MNI152 space corresponding to the regions of interest
        """
        revisedApiData = dict.fromkeys(['zscores', 'coords', 'samples', 'probes', 'specimen'])
        revisedApiData['zscores'] = []
        revisedApiData['coords'] = []
        revisedApiData['samples'] = []
        revisedApiData['probes'] = []
        revisedApiData['specimen'] = []
        dataImg = img.get_data()
        imgMni = img.affine
        invimgMni = inv(imgMni)
        Mni = specimen['alignment3d']
        T = np.dot(invimgMni, Mni)
        coords = transformSamples(apidataind['samples'], T)
        coords = np.rint(coords)       
        for i in range(0, len(coords)):
            coord = coords[i]
            sum = (coord > 0).sum()
            if sum != 3 or dataImg[int(coord[0]), int(coord[1]), int(coord[2])] <= self.mapthreshold or dataImg[int(coord[0]),int(coord[1]),int(coord[2])] == 0:
                coords[i] = [-1, -1, -1]
        for i in range(0, len(coords)):
            coord = coords[i]
            sum = (coord > 0).sum()
            if sum == 3:
                revisedApiData['zscores'].append(apidataind['zscores'][i])
                revisedApiData['coords'].append(coord)
        revisedApiData['samples'] = apidataind['samples'][:]
        revisedApiData['probes'] = apidataind['probes'][:]
        revisedApiData['specimen'] = specimen['name']
        return revisedApiData



    def queryapipartial(self, donorId):
        """
        Query Allen Brain Api for the given set of genes
        """
        url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        for p in self.probeids:
            url += p
            url += ","
        url = url[:-1]
        url += "][donors$eq"
        url += donorId
        url += "]"
        if(self.verboseflag):
            print(url)
        response = urllib.request.urlopen(url).read().decode('utf8')
        text = json.loads(response)
        data = text['msg']
        samples = []
        probes = []

        if not os.path.exists(self.cache):
            os.makedirs(self.cache)
        donorpath = os.path.join(self.cache, donorId)
        if not os.path.exists(donorpath):
            os.makedirs(donorpath)
        nsamples = len(data['samples'])
        nprobes = len(data['probes'])
        samples = data['samples']
        probes = data['probes']

        zscores = np.zeros((nsamples, nprobes))
        for i in range(0, nprobes):
            for j in range(0, nsamples):
                zscores[j][i] = probes[i]['z-score'][j]

        #LOAD PROBES
        fileName = os.path.join(donorpath, 'probes.txt')
        f = open(fileName, "r")
        probesC = json.load(f)
        f.close()
        #LOAD ZSCORES
        fileName = os.path.join(donorpath, 'zscores.txt')
        zscoresC = np.loadtxt(fileName)
        '''
        mri1 = np.array([s['sample']['mri'] for s in samples])
        mri2 = np.array([s['sample']['mri'] for s in samplesC])
        if not np.equal(mri1.all(), mri2.all()):
            print('mri1 and mri2 are different')
        np.savetxt('mri1.txt', mri1)
        np.savetxt('mri2.txt', mri2)
        exit()
        '''
        '''
        apiDataC = dict()
        apiDataC['samples'] = samples
        apiDataC['probes'] = probesC + probes
        apiDataC['zscores'] = np.append(zscoresC, zscores, axis=1)
        '''
        fileName = os.path.join(donorpath, 'samples.txt')
        with open(fileName, 'w') as outfile:
            json.dump(samples, outfile)
        probes = probesC + probes
        fileName = os.path.join(donorpath, 'probes.txt')
        with open(fileName, 'w') as outfile:
            json.dump(probes, outfile)
        zscores = np.append(zscoresC, zscores, axis=1)
        filename = os.path.join(donorpath, 'zscores.txt')
        np.savetxt(filename, zscores)

        '''
        if  not os.path.exists(filename):
            np.savetxt(filename, zscores)
        else:
            zscoresE = np.loadtxt(filename)
            print(zscoresE.shape)
            zscoresE.append(zscores)
            np.savetxt(filename, zscores)


        apiData = dict()
        apiData['samples'] = samples
        apiData['probes'] = probes
        apiData['zscores'] = zscores
        print('For ',donorId,' samples_length: ',len(apiData['samples']),' probes_length: ',len(apiData['probes']),' zscores_shape: ',apiData['zscores'].shape)
        return apiData
        '''

    def downloadspecimens(self):
        """
        Downlaod names and alignment matrix for each specimen/donor from Allen Brain Api and save them on disk as specimenName.txt
        and specimenMat.txt respectively, load.
        """
        specimens  = ['H0351.1015', 'H0351.1012', 'H0351.1016', 'H0351.2001', 'H0351.1009', 'H0351.2002']
        self.apidata['specimenInfo'] = []
        for i in range(0, len(specimens)):
            url = "http://api.brain-map.org/api/v2/data/Specimen/query.json?criteria=[name$eq"
            url+= "'"
            url += specimens[i]
            url += "']&include=alignment3d"
            if(self.verboseflag):
                print(url)
            if sys.version_info[0] >= 3:
                response = urllib.request.urlopen(url).read().decode('utf8')
                text = json.loads(response)
            else:
                text = requests.get(url).json()
            data = text['msg'][0]
            res = getSpecimenData(data)
            self.apidata['specimenInfo'].append(res)
        if(self.verboseflag):
            print(self.apidata['specimenInfo'])
        for i in range(0, len(self.donorids)):
            factorPath = os.path.join(self.cache, self.donorids[i]+'/specimenName.txt')
            with open(factorPath, 'w') as outfile:
                outfile.write(self.apidata['specimenInfo'][i]['name'])
            factorPath = os.path.join(self.cache, self.donorids[i]+'/specimenMat.txt')
            np.savetxt(factorPath, self.apidata['specimenInfo'][i]['alignment3d'])

    def getapidata(self):
        """
        Loop through the donors and call queryapi() and populate apidata
        """
        self.apidata['apiinfo'] = []
        for i in range(0, len(self.donorids)):
            self.apidata['apiinfo'].append(self.queryapi(self.donorids[i]))

    def download_and_retrieve_gene_data(self):
        """
        Download data from Allen Brain Api for the given set of genes and specimen information
        """
        self.getapidata()
        self.downloadspecimens()

    def set_candidate_genes(self, genelist):
        """
        Set list of genes and prepare to read/download data for them.
        """
        self.genelist = genelist
        self.downloadgenelist = self.genelist[:]
        if not os.path.exists(self.cache):
            if(self.verboseflag):
                print('self.cache does not exist')
            self.retrieveprobeids()
            self.download_and_retrieve_gene_data()
        else:
            self.creategenecache()
            for k, v in self.genecache.items():
                if k in self.genelist:
                    self.downloadgenelist.remove(k)            
            if self.downloadgenelist:
                print('Microarray expression values of',len(self.downloadgenelist),'gene(s) need(s) to be downloaded')
            if(self.verboseflag):
                print(self.downloadgenelist)
            self.retrieveprobeids()
            if self.downloadgenelist:
                for i in range(0, len(self.donorids)):
                    self.queryapipartial(self.donorids[i])
            self.readCachedApiSpecimenData()

    def run(self):
        self.performAnova()

    def performAnova(self):
        """
        Perform one way anova on zscores as the dependent variable and specimen factors such as age, race, name and area
        as independent variables
        """
        area1_zscores = []
        area2_zscores = []
        area1_specimen = []
        area2_specimen = []
        area1_area = []
        area2_area = []
        combined_zscores = []
        if(self.verboseflag):
            print(" ",len(self.main_r)," ",self.main_r[0]['name']," ",self.main_r[1]['name'])
        for i in range(0, len(self.main_r)):
            if self.main_r[i]['name'] == 'img1':
                area1_zscores = area1_zscores + self.main_r[i]['zscores'][:]
                area1_specimen = area1_specimen + [self.main_r[i]['specimen']]*len(self.main_r[i]['zscores'])
                area1_area = area1_area + [self.main_r[i]['name']]*len(self.main_r[i]['zscores'])
                combined_zscores = combined_zscores + self.main_r[i]['zscores'][:]
            elif self.main_r[i]['name'] == 'img2':
                area2_zscores = area2_zscores + self.main_r[i]['zscores'][:]
                area2_specimen = area2_specimen + [self.main_r[i]['specimen']]*len(self.main_r[i]['zscores'])
                area2_area = area2_area + [self.main_r[i]['name']]*len(self.main_r[i]['zscores'])
                combined_zscores = combined_zscores + self.main_r[i]['zscores'][:]
        factor_age = []
        factor_race = []
        factor_gender = []
        factor_age_numeric = []
        factor_specimen = area1_specimen + area2_specimen
        factor_area = area1_area + area2_area
        if(self.verboseflag):
            print(factor_specimen)
            print(factor_area)

        n_samples = len(combined_zscores)
        n_samples_area1 = len(area1_area)
        n_samples_area2 = len(area2_area)
        if(self.verboseflag):
            print("some variables ",n_samples," , ",n_samples_area1," , ",n_samples_area2, " , ", len(factor_specimen))
        specimenFactors = readSpecimenFactors(self.cache)
        if(self.verboseflag):
            print("number of specimens ", len(specimenFactors), " name: ", len(specimenFactors['name']))

        st = set(specimenFactors['name'])
        for ind, a in enumerate(factor_specimen):
            info_index = 0
            if a in st:
                info_index = specimenFactors['name'].index(a)
            factor_age_numeric = factor_age_numeric + [specimenFactors['age'][info_index]]
            factor_race = factor_race + [specimenFactors['race'][info_index]]

        if(self.verboseflag):
            print('race')
            print(factor_race)
            print('age')
            print(factor_age_numeric)
            print(len(self.genelist))

        allProbeData = getmeanzscores(self.genesymbols, combined_zscores, len(area1_zscores), len(area2_zscores))
        combined_zscores = np.zeros((len(allProbeData['combined_zscores']), len(allProbeData['combined_zscores'])))
        combined_zscores = np.copy(allProbeData['combined_zscores'])
        area1_zscores = np.zeros((len(allProbeData['area1_zscores']), len(allProbeData['area1_zscores'])))
        area1_zscores = np.copy(allProbeData['area1_zscores'])
        area2_zscores = np.zeros((len(allProbeData['area2_zscores']), len(allProbeData['area2_zscores'])))
        area2_zscores = np.copy(allProbeData['area2_zscores'])
        if(self.verboseflag):
            print('combined_zscores shape ',combined_zscores.shape,' ',area1_zscores.shape,' ',area2_zscores.shape)
        uniqueId = np.copy(allProbeData['uniqueId'])

        geneIds = []
        st = set(self.genesymbols)
        for ind, a in enumerate(uniqueId):
            index = 0
            if a in st:
                index = self.genesymbols.index(a)
            geneIds = geneIds + [self.genesymbols[index]]

        n_genes = len(combined_zscores[0]) #SHOULD NOT THIS BE 285???
        Reference_Anovan_p = np.zeros(n_genes)
        Reference_Anovan_eta2 = np.zeros(n_genes)
        Reference_Anovan_CI_l = np.zeros(n_genes)
        Reference_Anovan_CI_h = np.zeros(n_genes)
        Reference_Anovan_diff_mean = np.zeros(n_genes)
        F_vec_ref_anovan = np.zeros(n_genes)
        data = {}
        data['Area'] = factor_area
        data['Specimen'] = factor_specimen
        data['Age'] = factor_age_numeric
        data['Race'] = factor_race

        for i in range(0, n_genes):
            data['Zscores'] = combined_zscores[:,i]
            mod = ols('Zscores ~ Area + Specimen + Age + Race', data=data).fit()
            aov_table = sm.stats.anova_lm(mod, typ=1)
            if(self.verboseflag):
                print(aov_table)

            F_vec_ref_anovan[i] = aov_table['F'][0]
            ss_total = aov_table['sum_sq'][0]+aov_table['sum_sq'][1]+aov_table['sum_sq'][4]
            ss_between_group_area =  aov_table['sum_sq'][0]
            Reference_Anovan_eta2[i] = ss_between_group_area/ss_total

            var1 = []
            var2 = []
            row = combined_zscores[:,i]
            var1 = var1 + [combined_zscores[j][i] for j in range(0, len(row)) if factor_area[j] == 'img1']
            var2 = var2 + [combined_zscores[j][i] for j in range(0, len(row)) if factor_area[j] == 'img2']
            mse = (np.var(var1, ddof=1) + np.var(var2, ddof=1))*0.5
            sm1m2 = 2.011*sqrt((2*mse)/n_genes)
            mean1 = np.mean(var1)
            mean2 = np.mean(var2)
            v = mean1 - mean2
            Reference_Anovan_CI_l[i] = v - sm1m2
            Reference_Anovan_CI_h[i] = v + sm1m2
            Reference_Anovan_diff_mean[i] = v
            Reference_Anovan_p[i] = aov_table['PR(>F)'][0] #p(1)

        n_rep = 1000
        invn_rep = 1/n_rep
        FWE_corrected_p = np.zeros(n_genes)
        F_mat_perm_anovan = np.zeros((n_rep, n_genes))
        p_mat_perm_anovan = np.zeros((n_rep, n_genes))
        F_mat_perm_anovan[0] = F_vec_ref_anovan
        for rep in range(1, n_rep):
            F_vec_perm_anovan = np.zeros(n_genes)
            p_vec_perm_anovan = np.zeros(n_genes)
            for j in range(0, n_genes):
                shuffle = np.random.permutation(factor_area)
                data['Area'] = shuffle
                data['Zscores'] = combined_zscores[:,j]
                mod = ols('Zscores ~ Area + Specimen + Age + Race', data=data).fit()
                aov_table = sm.stats.anova_lm(mod, typ=1)
                F_vec_perm_anovan[j] = aov_table['F'][0]
                p_vec_perm_anovan[j] = aov_table['PR(>F)'][0]
            F_mat_perm_anovan[rep] = F_vec_perm_anovan
            p_mat_perm_anovan[rep] = p_vec_perm_anovan

        ref = F_mat_perm_anovan.max(1)
        Uncorrected_permuted_p = np.zeros(n_genes)

        for j in range(0, n_genes):
            val = F_vec_ref_anovan[j]
            list1 = F_mat_perm_anovan[:,j]
            sum = len([1 for a in list1 if a >= val])
            Uncorrected_permuted_p[j] = sum*invn_rep
            sum = len([1 for a in ref if a >= val])
            FWE_corrected_p[j] = sum*invn_rep
        self.result = dict(zip(geneIds, FWE_corrected_p))
        if(self.verboseflag):
            print(self.result)


    def pvalues(self):
        return self.result;
