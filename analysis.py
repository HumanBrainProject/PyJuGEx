# -*- coding: utf-8 -*-
import os
import numpy as np
from numpy import *
import json
from numpy.linalg import inv
import rpy2.robjects as robjects
import rpy2.rlike.container as rlc
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import urllib.request
import scipy as sp
import scipy.stats.mstats
import xmltodict

def switch2gensymbol(entrez_id, combined_zscores, area1len, area2len):
    unique_entrez_id = np.unique(entrez_id)
    print(entrez_id)
    print(unique_entrez_id)
    indices = [np.where(np.in1d(entrez_id, x))[0] for x in unique_entrez_id]
    print(indices)
    winsorzed_mean_zscores = np.zeros((len(combined_zscores), len(unique_entrez_id)))
    print(len(combined_zscores),' ',len(combined_zscores[0]),' ',len(indices))
    tmp = [[] for _ in range(len(combined_zscores))]
    for i in range (0, len(unique_entrez_id)):
        for j in range(0, len(combined_zscores)):
            tmp[j] = [combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))]
            winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1)) #maybe add a special case for one element in tmp
    print(len(winsorzed_mean_zscores),' ',len(winsorzed_mean_zscores[0]))
    res = dict.fromkeys(['uniqueId', 'combined_zscores', 'area1_zscores', 'area2_zscores'])
    res['uniqueId'] = unique_entrez_id
    res['combined_zscores'] = winsorzed_mean_zscores
    res['area1_zscores'] = winsorzed_mean_zscores[0:area1len]
    res['area2_zscores'] = winsorzed_mean_zscores[area1len:]
    return res

def switch2genesymbol(gene_symbols, combined_zscores, area1len, area2len):
    unique_gene_symbols = np.unique(gene_symbols)
    indices = [np.where(np.in1d(gene_symbols, x))[0] for x in unique_gene_symbols]
    print(indices)
    winsorzed_mean_zscores = np.zeros((len(combined_zscores), len(unique_gene_symbols)))
    print(len(combined_zscores),' ',len(combined_zscores[0]),' ',len(indices))
    tmp = [[] for _ in range(len(combined_zscores))]
    for i in range (0, len(unique_gene_symbols)):
        for j in range(0, len(combined_zscores)):
            tmp[j] = [combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))]
            winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1)) #maybe add a special case for one element in tmp
    print(len(winsorzed_mean_zscores),' ',len(winsorzed_mean_zscores[0]))
    res = dict.fromkeys(['uniqueId', 'combined_zscores', 'area1_zscores', 'area2_zscores'])
    res['uniqueId'] = unique_gene_symbols
    res['combined_zscores'] = winsorzed_mean_zscores
    res['area1_zscores'] = winsorzed_mean_zscores[0:area1len]
    res['area2_zscores'] = winsorzed_mean_zscores[area1len:]
    return res

def getSpecimenData(info):
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
    np_T = np.array(T[0:3, 0:4])
    mri = np.vstack(s["sample"]["mri"] for s in samples)
    add = np.ones((len(mri), 1), dtype=np.int)
    mri = np.append(mri, add, axis=1)
    mri = np.transpose(mri)
    coords = np.matmul(np_T, mri)
    coords = coords.transpose()
    return coords

#NO NEED TO DOWNLOAD AT ALL, READ FROM DISK

def buildSpecimenFactors():
    url = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"
    specimenFactors = dict()
    specimenFactors['id'] = []
    specimenFactors['name'] = []
    specimenFactors['race'] = []
    specimenFactors['gender'] = []
    specimenFactors['age'] = []
    response = urllib.request.urlopen(url).read().decode('utf8')
    text = json.loads(response)
    rootDir = os.path.dirname('./')
    factorPath = os.path.join(rootDir, 'specimenFactors.txt')
    with open(factorPath, 'w') as outfile:
        json.dump(text, outfile)
    res = text['msg']
    for i in range(0, len(res)):
        specimenFactors['id'].append(res[i]['id'])
        specimenFactors['name'].append(res[i]['name'])
        specimenFactors['race'].append(res[i]['race_only'])
        specimenFactors['gender'].append(res[i]['sex'])
        specimenFactors['age'].append(res[i]['age']['days']/365)
    return specimenFactors;

def readSpecimenFactors():
    specimenFactors = dict.fromkeys(['id', 'name', 'race', 'gender', 'age'])
    specimenFactors['id'] = []
    specimenFactors['name'] = []
    specimenFactors['race'] = []
    specimenFactors['gender'] = []
    specimenFactors['age'] = []
    rootDir = os.path.dirname('./')
    fileName = os.path.join(rootDir, 'specimenFactors.txt')
    if not os.path.exists(fileName):
        specimenFactors = buildSpecimenFactors()
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
    return specimenFactors;

class Analysis:
    def __init__(self, gene_cache):
        self.refreshcache = False
        if not gene_cache:
            self.refreshcache = True
        self.cache = gene_cache
        self.probeids = []
        self.probeidsD = dict.fromkeys(['gene_symbols', 'probe_ids'])
        self.probeidsD['gene_symbols'] = []
        self.probeidsD['probe_ids'] = []
        self.donorids = ['15496','14380','15697','9861','12876','10021'] #HARDCODING DONORIDS
        self.genelist = None #Generate the empty dicts here"
        self.apidata = dict.fromkeys(['apiinfo', 'specimeninfo'])
        self.vois = []
        self.main_r = []
        self.mapthreshold = 2
        self.result = None

    def retrieveprobeids(self):
        genesymbols = np.unique(np.array(self.genelist['gene_symbols']))
        print(genesymbols)
        for g in genesymbols:
            url = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq"
            url += g
            url += "],rma::options[only$eq'probes.id']"
            print(url)
            response = urllib.request.urlopen(url).read()
            data = xmltodict.parse(response)
            for d in data['Response']['probes']['probe']:
                self.probeidsD['gene_symbols'] = self.probeidsD['gene_symbols'] + [g]
                self.probeidsD['probe_ids'] = self.probeidsD['probe_ids'] + [d['id']]
                self.probeids = self.probeids + [d['id']]

    def readCachedApiSpecimenData(self, rootdir):
        self.apidata['specimenInfo'] = []
        self.apidata['apiinfo'] = []
        for d in self.donorids:
            donorpath = os.path.join(rootdir, d)
            fileNameM = os.path.join(donorpath, 'specimenInfoMat.txt')
            mat = np.loadtxt(fileNameM)
            fileNameN = os.path.join(donorpath, 'specimenInfoName.txt')
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
            print('inside readcachedata ',len(apiDataC['samples']), ' ', apiDataC['zscores'].shape, ' ', len(apiDataC['probes']))
        for s in self.apidata['specimenInfo']:
            print(s['alignment3d'])
            print(s['name'])

    def set_ROI_MNI152(self, voi, index):
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
            else:
                revisedApiDataCombo['name'] = 'img2'
            print('extractexplevel img1: ',revisedApiDataCombo['specimen'],' ',len(revisedApiDataCombo['coords']))
            self.main_r.append(revisedApiDataCombo)

    def expressionSpmCorrelation(self, img, apidataind, specimen):
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
            if sum != 3 or dataImg[int(coord[0]), int(coord[1]), int(coord[2])] <= self.mapthreshold/10 or dataImg[int(coord[0]),int(coord[1]),int(coord[2])] == 0:
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
        for i in range(0, len(revisedApiData['coords'])):
            print(revisedApiData['coords'][i])
            print(revisedApiData['zscores'][i].shape)
        return revisedApiData

    def retrieve_gene_data(self, genelist):
        print('In retrieve_gene_data')
        self.genelist = genelist
        print(len(genelist))
        self.retrieveprobeids()

    def queryapi(self, donorId):
        url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        for p in self.probeids:
            url += p
            url += ","
        url = url[:-1]
        url += "][donors$eq"
        url += donorId
        url += "]"
        response = urllib.request.urlopen(url).read().decode('utf8')
        text = json.loads(response)
        data = text['msg']
        samples = []
        probes = []
        rootDir = os.path.dirname('AllenBrainApi/')
        if not os.path.exists(rootDir):
            os.makedirs(rootDir)
        donorPath = os.path.join(rootDir, donorId)
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
        f = open(fileName, "r")
        probesC = json.load(f)
        f.close()

        fileName = os.path.join(donorPath, 'zscores.txt')
        np.savetxt(fileName, zscores)
        zscoresC = np.loadtxt(fileName)

        fileName = os.path.join(donorPath, 'samples.txt')
        f = open(fileName, "r")
        samplesC = json.load(f)
        f.close()

        fileName = os.path.join(donorPath, 'probes.txt')
        f = open(fileName, "r")
        probesC = json.load(f)
        f.close()

        fileName = os.path.join(donorPath, 'zscores.txt')
        zscoresC = np.loadtxt(fileName)

        apiData = dict()
        apiData['samples'] = samples
        apiData['probes'] = probes
        apiData['zscores'] = zscores
        return apiData

    def downloadspecimens(self):
        specimens  = ['H0351.1015', 'H0351.1012', 'H0351.1016', 'H0351.2001', 'H0351.1009', 'H0351.2002']
        self.apidata['specimenInfo'] = []
        for i in range(0, len(specimens)):
            url = "http://api.brain-map.org/api/v2/data/Specimen/query.json?criteria=[name$eq"
            url+= "'"
            url += specimens[i]
            url += "']&include=alignment3d"
            print(url)
            response = urllib.request.urlopen(url).read().decode('utf8')
            text = json.loads(response)
            data = text['msg'][0]
            res = getSpecimenData(data)
            self.apidata['specimenInfo'].append(res)
        print(len(self.apidata['specimenInfo']))

    def getapidata(self):
        self.apidata['apiinfo'] = []
        for i in range(0, len(self.donorids)):
            self.apidata['apiinfo'].append(self.queryapi(self.donorids[i]))

    def download_and_retrieve_gene_data(self):
       self.downloadspecimens()
       self.apidata['apiinfo'] = []
       print('In downlaod_and_retrieve_gene_data')
       rootDir = os.path.dirname('AllenBrainApi/')
       self.getapidata()

    def set_candidate_genes(self, genelist):
        self.genelist = genelist
        print(len(genelist))
        self.retrieveprobeids()
        if self.refreshcache is False:
            self.readCachedApiSpecimenData(os.path.dirname(self.cache))
        else:
            self.download_and_retrieve_gene_data()

    def run(self):
        self.performAnova()

    def performAnova(self):
        r = robjects.r
        area1_zscores = []
        area2_zscores = []
        area1_specimen = []
        area2_specimen = []
        area1_area = []
        area2_area = []
        combined_zscores = []
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
        print(factor_specimen)
        print(factor_area)

        n_samples = len(combined_zscores)
        n_samples_area1 = len(area1_area)
        n_samples_area2 = len(area2_area)

        print("some variables ",n_samples," , ",n_samples_area1," , ",n_samples_area2, " , ", len(factor_specimen))
        specimenFactors = readSpecimenFactors()
        print("number of specimens ", len(specimenFactors), " name: ", len(specimenFactors['name']))

        st = set(specimenFactors['name'])
        for ind, a in enumerate(factor_specimen):
            info_index = 0
            if a in st:
                info_index = specimenFactors['name'].index(a)
            factor_age_numeric = factor_age_numeric + [specimenFactors['age'][info_index]]
            factor_race = factor_race + [specimenFactors['race'][info_index]]
        print('race')
        print(factor_race)
        print('age')
        print(factor_age_numeric)
        print(len(self.genelist))
        allProbeData = switch2genesymbol(self.genelist['gene_symbols'], combined_zscores, len(area1_zscores), len(area2_zscores))
        '''
        allProbeData = switch2gensymbol(self.genelist['entrez_id'], combined_zscores, len(area1_zscores), len(area2_zscores))
        allProbeDataD = switch2genesymbol(self.genelist['gene_symbols'], combined_zscores, len(area1_zscores), len(area2_zscores))
        if np.equal(np.array(allProbeData['combined_zscores']).all(), np.array(allProbeDataD['combined_zscores']).all()):
            print('samecz')
        if np.equal(np.array(allProbeData['area1_zscores']).all(), np.array(allProbeDataD['area1_zscores']).all()):
            print('samearea1')
        if np.equal(np.array(allProbeData['area2_zscores']).all(), np.array(allProbeDataD['area2_zscores']).all()):
            print('samearea2')
        exit()
        '''
        combined_zscores = np.zeros((len(allProbeData['combined_zscores']), len(allProbeData['combined_zscores'])))
        combined_zscores = np.copy(allProbeData['combined_zscores'])
        area1_zscores = np.zeros((len(allProbeData['area1_zscores']), len(allProbeData['area1_zscores'])))
        area1_zscores = np.copy(allProbeData['area1_zscores'])
        area2_zscores = np.zeros((len(allProbeData['area2_zscores']), len(allProbeData['area2_zscores'])))
        area2_zscores = np.copy(allProbeData['area2_zscores'])
        print('combined_zscores shape ',combined_zscores.shape,' ',area1_zscores.shape,' ',area2_zscores.shape)
        uniqueId = np.copy(allProbeData['uniqueId'])
        geneIds = []
        st = set(self.genelist['gene_symbols'])
        for ind, a in enumerate(uniqueId):
            index = 0
            if a in st:
                index = self.genelist['gene_symbols'].index(a)
            geneIds = geneIds + [self.genelist['gene_symbols'][index]]

        n_genes = len(combined_zscores[0]) #SHOULD NOT THIS BE 285???
        print(n_genes)
        Reference_Anovan_p = np.zeros(n_genes)
        Reference_Anovan_eta2 = np.zeros(n_genes)
        Reference_Anovan_CI_l = np.zeros(n_genes)
        Reference_Anovan_CI_h = np.zeros(n_genes)
        Reference_Anovan_diff_mean = np.zeros(n_genes)
        F_vec_ref_anovan = np.zeros(n_genes)

        for i in range(0, n_genes):
            od = rlc.OrdDict([('Area', robjects.StrVector(factor_area)),
                              ('Specimen', robjects.StrVector(factor_specimen)),
                              ('Age', robjects.IntVector(factor_age_numeric)),
                              ('Race', robjects.StrVector(factor_race)),
                              ('Zscores', robjects.FloatVector(combined_zscores[:,i]))])
            dataf = robjects.DataFrame(od)
            f = robjects.Formula('Zscores~Area+Specimen+Age+Race')
            a = r['aov'](f, data = dataf)
            summary = r['summary'](a)
            print(summary)
            F_vec_ref_anovan[i] = summary[0][3][0] #tab{2,6}
            ss_total = summary[0][1][0]+summary[0][1][1]+summary[0][1][2] #tab{7,2}
            ss_between_group_area = summary[0][1][0] #tab{2,2}
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
            Reference_Anovan_p[i] = summary[0][4][0] #p(1)

        n_rep = 1000
        FWE_corrected_p = np.zeros(n_genes)
        F_mat_perm_anovan = np.zeros((n_rep, n_genes))
        p_mat_perm_anovan = np.zeros((n_rep, n_genes))
        F_mat_perm_anovan[0] = F_vec_ref_anovan
        for rep in range(1, n_rep):
            F_vec_perm_anovan = np.zeros(n_genes)
            p_vec_perm_anovan = np.zeros(n_genes)
            for j in range(0, n_genes):
                shuffle = np.random.permutation(factor_area)
                f = robjects.StrVector(shuffle)
                od = rlc.OrdDict([('Area', f),
                                  ('Specimen', robjects.StrVector(factor_specimen)),
                                  ('Age', robjects.IntVector(factor_age_numeric)),
                                  ('Race', robjects.StrVector(factor_race)),
                                  ('Zscores', robjects.FloatVector(combined_zscores[:,j]))])
                dataf = robjects.DataFrame(od)
                f = robjects.Formula('Zscores~Area+Specimen+Age+Race')
                a = r['aov'](f, data = dataf)
                summary = r['summary'](a)
                F_vec_perm_anovan[j] = summary[0][3][0]
                p_vec_perm_anovan[j] = summary[0][4][0]                
            F_mat_perm_anovan[rep] = F_vec_perm_anovan
            p_mat_perm_anovan[rep] = p_vec_perm_anovan

        ref = F_mat_perm_anovan.max(1)
        Uncorrected_permuted_p = np.zeros(n_genes)

        for j in range(0, n_genes):
            val = F_vec_ref_anovan[j]
            list1 = F_mat_perm_anovan[:,j]
            sum = len([1 for a in list1 if a >= val])
            Uncorrected_permuted_p[j] = sum/n_rep
            sum = len([1 for a in ref if a >= val])
            FWE_corrected_p[j] = sum/n_rep
        self.result = dict(zip(geneIds, FWE_corrected_p))
        print(self.result)

    def pvalues(self):
        return self.result;
