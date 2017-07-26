# -*- coding: utf-8 -*-
import backendpackage
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

def switch2gensymbol(entrez_id, combined_zscores, area1Len, area2Len):
    unique_entrez_id = np.unique(entrez_id)
    print('length of entrez id is ', len(unique_entrez_id))
    print(unique_entrez_id)
    inds =  [[] for _ in range(len(unique_entrez_id))]
    for i in range(0, len(unique_entrez_id)):
        for j in range (0, len(entrez_id)):
            if unique_entrez_id[i] == entrez_id[j]:
                inds[i].append(j);
    winsorzed_mean_zscores = np.zeros((len(combined_zscores), len(unique_entrez_id)))
    for i in range (0, len(unique_entrez_id)):
        tmp = [[] for _ in range(len(combined_zscores))]
        for j in range(0, len(combined_zscores)):
            for k in range (0, len(inds[i])):
                tmp[j].append(combined_zscores[j][inds[i][k]]);
            ncols = len(tmp[j]);
        for j in range(0, len(combined_zscores)):
            winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1))#maybe add a special case for one element in tmp
    print(len(winsorzed_mean_zscores),' ',len(winsorzed_mean_zscores[0]))
    res = dict()
    res['uniqueId'] = []
    for i in range(0, len(unique_entrez_id)):
        res['uniqueId'].append(unique_entrez_id[i])
    res['combined_zscores'] = []
    for i in range(0, len(winsorzed_mean_zscores)):
        res['combined_zscores'].append(winsorzed_mean_zscores[i])
    print(len(res['combined_zscores']),' ',len(res['combined_zscores'][0]))
    res['area1_zscores'] = []
    res['area2_zscores'] = []
    for i in range(0, area1Len):
           res['area1_zscores'].append(winsorzed_mean_zscores[i])
    for i in range(0, area2Len):
           res['area2_zscores'].append(winsorzed_mean_zscores[i+area1Len])
    return res


def transformSamples(samples, T):
    nsamples = len(samples)
    coords = []
    T00 = T[0][0]
    T01 = T[0][1]
    T02 = T[0][2]
    T03 = T[0][3]
    T10 = T[1][0]
    T11 = T[1][1]
    T12 = T[1][2]
    T13 = T[1][3]
    T20 = T[2][0]
    T21 = T[2][1]
    T22 = T[2][2]
    T23 = T[2][3]
    for i in range(0, nsamples):
        mri = samples[i]['sample']['mri']
        x = mri[0]
        y = mri[1]
        z = mri[2]
        cx = T00*x + T01*y + T02*z + T03
        cy = T10*x + T11*y + T12*z + T13
        cz = T20*x + T21*y + T22*z + T23
        coord = [cx, cy, cz]
        coords.append(coord)
    print(len(coords))
    return coords

#NO NEED TO DOWNLOAD AT ALL, READ FROM DISK

def buildSpecimenFactors():
    url = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"
    specimenFactors = dict()
    response = urllib.request.urlopen(url).read().decode('utf8')
    text = json.loads(response)
    rootDir = os.path.dirname('./')
    if not os.path.exists(rootDir):
        os.makedirs(rootDir)
    factorPath = os.path.join(rootDir, 'specimenFactors.txt')
    with open(factorPath, 'w') as outfile:
        json.dump(text, outfile)
    res = text['msg']
    specimenFactors['id'] = []
    specimenFactors['name'] = []
    specimenFactors['race'] = []
    specimenFactors['gender'] = []
    specimenFactors['age'] = []
    for i in range(0, len(res)):
        specimenFactors['id'].append(res[i]['id'])
        specimenFactors['name'].append(res[i]['name'])
        specimenFactors['race'].append(res[i]['race_only'])
        specimenFactors['gender'].append(res[i]['sex'])
        specimenFactors['age'].append(res[i]['age']['days']/365)
    return specimenFactors;

def readSpecimenFactors():
    specimenFactors = dict()
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
    def __init__(self):
        self.genelist = None #Generate the empty dicts here"
        self.apispecimendata = None #Generate the empty dicts here"
        self.vois = []
        self.main_r = []
        self.mapthreshold = 2
        self.result = None

    def readCachedApiSpecimenData(self, donorIds, probeIds, rootdir):
        self.apidata = dict.fromkeys(['apiinfo', 'specimeninfo'])
        self.apidata['specimenInfo'] = []
        self.apidata['apiinfo'] = []
        for d in donorIds:
            donorPath = os.path.join(rootdir, d)
            fileNameM = os.path.join(donorPath, 'specimenInfoMat.txt')
            mat = np.loadtxt(fileNameM)
            fileNameN = os.path.join(donorPath, 'specimenInfoName.txt')
            f = open(fileNameN, 'r')
            name = f.read()
            f.close()
            specimen = dict.fromkeys(['name', 'alignment3d'])
            specimen['name'] = name
            specimen['alignment3d'] = mat
            self.apidata['specimenInfo'].append(specimen)
            #LOAD SAMPLES
            fileName = os.path.join(donorPath, 'samples.txt')
            f = open(fileName, "r")
            samplesC = json.load(f)
            f.close()
            #LOAD PROBES
            fileName = os.path.join(donorPath, 'probes.txt')
            f = open(fileName, "r")
            probesC = json.load(f)
            f.close()
            #LOAD ZSCORES
            fileName = os.path.join(donorPath, 'zscores.txt')
            zscoresC = np.loadtxt(fileName)
            #LOAD EXPLEVELS
            fileName = os.path.join(donorPath, 'expressionlevels.txt')
            explevelsC = np.loadtxt(fileName)
            apiDataC = dict()
            apiDataC['samples'] = samplesC
            apiDataC['probes'] = probesC
            apiDataC['explevels'] = explevelsC
            apiDataC['zscores'] = zscoresC
            self.apidata['apiinfo'].append(apiDataC)
            print('inside readcachedata ',len(apiDataC['samples']), ' ', apiDataC['zscores'].shape, ' ', apiDataC['explevels'].shape, ' ', len(apiDataC['probes']))
        for s in self.apidata['specimenInfo']:
            print(s['alignment3d'])
            print(s['name'])

    def set_coordinates_region(self, voi, index):
        for i in range(0, len(self.apidata['specimenInfo'])):
            revisedApiDataCombo = dict()
            revisedApiData = self.expressionSpmCorrelation(voi, self.apidata['apiinfo'][i], self.apidata['specimenInfo'][i]) #Maybe an index will work instead of expressionspmcorrelation
            revisedApiDataCombo['explevels'] = revisedApiData['explevels'][:]
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
        revisedApiData = dict()
        revisedApiData['explevels'] = []
        revisedApiData['zscores'] = []
        revisedApiData['coords'] = []
        revisedApiData['samples'] = []
        revisedApiData['probes'] = []
        revisedApiData['specimen'] = []
        Mni = specimen['alignment3d']
        data1 = img.get_data()
        imgMni = img.affine
        invimgMni = inv(imgMni)
        Mni = specimen['alignment3d']
        T = np.dot(invimgMni, Mni)
        coords = transformSamples(apidataind['samples'], T)
        for i in range(0, len(coords)):
            coords[i] = np.rint(coords[i])
        print("After transform samples ")
        for i in range(0, len(coords)):
            coord = coords[i]
            sum = 0
            for j in range(0, len(coord)):
                if(coord[j] > 0):
                    sum += 1
            if(sum != 3 or data1[int(coord[0]),int(coord[1]),int(coord[2])] <= self.mapthreshold/10 or data1[int(coord[0]),int(coord[1]),int(coord[2])] == 0):
                coords[i][0] = -1 #IS IT 1 or -1
                coords[i][1] = -1
                coords[i][2] = -1
        for i in range(0, len(coords)):
            coord = coords[i]
            sum = 0
            for j in range(0, len(coord)) :
                if(coord[j] > 0):
                    sum += 1
            if(sum == 3):
                revisedApiData['explevels'].append(apidataind['explevels'][i])
                revisedApiData['zscores'].append(apidataind['zscores'][i])
                revisedApiData['coords'].append(coord)
        revisedApiData['samples'] = apidataind['samples'][:]
        revisedApiData['probes'] = apidataind['probes'][:]
        revisedApiData['specimen'] = specimen['name']
        for i in range(0, len(revisedApiData['coords'])):
            print(revisedApiData['coords'][i])
            print(revisedApiData['zscores'][i].shape)
        return revisedApiData

    def retrieve_gene_data(self, genelist, cache):
        print('In retrieve_gene_data')
        self.genelist = genelist
        print(len(genelist))

    def set_candidate_genes(self, genelist, cache):
        self.genelist = genelist
        print(len(genelist['entrez_id']))
        donorIds = ['15496','14380','15697','9861','12876','10021'] #HARDCODING DONORIDS
        self.readCachedApiSpecimenData(donorIds, genelist['probe_id'], os.path.dirname(cache))

    def run(self):
        self.performAnova()

    def performAnova(self):
        r = robjects.r
        area1_zscores = []
        area2_zscores = []
        area1_specimen = []
        area2_specimen = []
        area2_zscores = []
        area1_area = []
        area2_area = []
        area1_koords = []
        area2_koords = []
        combined_zscores = []
        print(" ",len(self.main_r)," ",self.main_r[0]['name']," ",self.main_r[1]['name'])
        for i in range(0, len(self.main_r)):
            if(self.main_r[i]['name'] == 'img1'):
                for j in range(0, len(self.main_r[i]['zscores'])):
                    area1_zscores.append(self.main_r[i]['zscores'][j])
                    #combined_zscores.append(main_r[i]['zscores'][j])
                    area1_specimen.append(self.main_r[i]['specimen'])
                    area1_area.append(self.main_r[i]['name'])
                    area1_koords.append(self.main_r[i]['coords'])
            elif(self.main_r[i]['name'] == 'img2'):
                for j in range(0, len(self.main_r[i]['zscores'])):
                    area2_zscores.append(self.main_r[i]['zscores'][j])
                    #combined_zscores.append(main_r[i]['zscores'][j])
                    area2_specimen.append(self.main_r[i]['specimen'])
                    area2_area.append(self.main_r[i]['name'])
                    area2_koords.append(self.main_r[i]['coords'])
        for i in range(0, len(area1_zscores)):
            combined_zscores.append(area1_zscores[i])
        for i in range(0, len(area2_zscores)):
            combined_zscores.append(area2_zscores[i])

        factor_specimen = []
        factor_area = []
        factor_age = []
        factor_race = []
        factor_gender = []
        factor_gender_binary = []
        factor_age_numeric = []
        factor_age_grouped = []
        factor_age_grouped_numeric = []
        for i in range(0, len(area1_specimen)):
            factor_specimen.append(area1_specimen[i])
        for i in range(0, len(area2_specimen)):
            factor_specimen.append(area2_specimen[i])
        for i in range(0, len(area1_area)):
            factor_area.append(area1_area[i])
        for i in range(0, len(area2_area)):
            factor_area.append(area2_area[i])
            #VERIFY THIS
        print(factor_specimen)
        print(factor_area)

        n_samples = len(combined_zscores)
        n_samples_area1 = len(area1_area)
        n_samples_area2 = len(area2_area)

        print("some variables ",n_samples," , ",n_samples_area1," , ",n_samples_area2, " , ", len(factor_specimen))
        specimenFactors = readSpecimenFactors()
        print("number of specimens ", len(specimenFactors), " name: ", len(specimenFactors['name']))
        for counter in range(0, n_samples):
            #info_index=find(strcmp({specimenFactors.name},factor_specimen(counter)));
            info_index = 0 #VERIFY THIS
            for j in range(0, len(specimenFactors['name'])):
                if(specimenFactors['name'][j] == factor_specimen[counter]): #CHECK THIS ENTRY
                    info_index=j
                    break
            factor_age_numeric.append(specimenFactors['age'][info_index])
            factor_race.append(specimenFactors['race'][info_index])
        print('race')
        print(factor_race)
        print('age')
        print(factor_age_numeric)
        print(len(self.genelist))
        allProbeData = switch2gensymbol(self.genelist['entrez_id'],combined_zscores, len(area1_zscores), len(area2_zscores))
        combined_zscores = np.zeros((len(allProbeData['combined_zscores']), len(allProbeData['combined_zscores'])))
        combined_zscores = np.copy(allProbeData['combined_zscores'])
        area1_zscores = np.zeros((len(allProbeData['area1_zscores']), len(allProbeData['area1_zscores'])))
        area1_zscores = np.copy(allProbeData['area1_zscores'])
        area2_zscores = np.zeros((len(allProbeData['area2_zscores']), len(allProbeData['area2_zscores'])))
        area2_zscores = np.copy(allProbeData['area2_zscores'])
        uniqueId = np.copy(allProbeData['uniqueId'])
        geneIds = []
        for i in range(0, len(uniqueId)):
            for j in range(0, len(self.genelist['entrez_id'])):
                if(self.genelist['entrez_id'][j] == uniqueId[i]):
                    geneIds.append(self.genelist['gene_symbol'][j])
                    break
        print(geneIds)
        n_genes = len(combined_zscores[0]) #SHOULD NOT THIS BE 285???
        print(n_genes)
        Reference_Anovan_p = np.zeros(n_genes)
        Reference_Anovan_eta2 = np.zeros(n_genes)
        Reference_Anovan_CI_l = np.zeros(n_genes)
        Reference_Anovan_CI_h = np.zeros(n_genes)
        Reference_Anovan_diff_mean = np.zeros(n_genes)
        F_vec_ref_anovan = np.zeros(n_genes)

        for i in range(0, n_genes):
            print(i+1)
            od = rlc.OrdDict([('Area', robjects.StrVector(factor_area)),
                              ('Specimen', robjects.StrVector(factor_specimen)),
                              ('Age', robjects.IntVector(factor_age_numeric)),
                              ('Race', robjects.StrVector(factor_race)),
                              ('Zscores', robjects.FloatVector(combined_zscores[:,i]))])
            dataf = robjects.DataFrame(od)
            #r['options'](contrasts=r['c']("contr.sum","contr.poly"))
            f = robjects.Formula('Zscores~Area+Specimen+Age+Race')
            a = r['aov'](f, data = dataf)
            #a = r['aov'](f, data = dataf,  type=3)
            summary = r['summary'](a)
            print(summary)
            #print('p is ', )
            F_vec_ref_anovan[i] = summary[0][3][0] #tab{2,6}
            #print(F_vec_ref_anovan[i])
            ss_total = summary[0][1][0]+summary[0][1][1]+summary[0][1][2] #tab{7,2}
            #print(ss_total)
            ss_between_group_area = summary[0][1][0] #tab{2,2}
            Reference_Anovan_eta2[i] = ss_between_group_area/ss_total
            var1 = []
            var2 = []
            row = combined_zscores[:,i]
            for j in range(0, len(row)):
                if factor_area[j] == 'img1':
                    var1.append(combined_zscores[j][i])
                if factor_area[j] == 'img2':
                    var2.append(combined_zscores[j][i])
            mse = (np.var(var1, ddof=1) + np.var(var2, ddof=1))*0.5
            sm1m2 = 2.011*sqrt((2*mse)/n_genes)
            #print('mse = ', mse, ' sm1m2 = ', sm1m2/2.011)
            mean1 = np.mean(var1)
            mean2 = np.mean(var2)
            v = mean1 - mean2
            Reference_Anovan_CI_l[i] = v - sm1m2
            Reference_Anovan_CI_h[i] = v + sm1m2
            Reference_Anovan_diff_mean[i] = v
            Reference_Anovan_p[i] = summary[0][4][0] #p(1)

        n_rep = 1000
        #n_rep = 1
        FWE_corrected_p = np.zeros(n_genes)
        F_mat_perm_anovan = np.zeros((n_rep, n_genes))
        p_mat_perm_anovan = np.zeros((n_rep, n_genes))
        for i in range(0, n_genes):
            F_mat_perm_anovan[0][i] = F_vec_ref_anovan[i]
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
                #print(summary)
                F_vec_perm_anovan[j] = summary[0][3][0]
                p_vec_perm_anovan[j] = summary[0][4][0]
                F_mat_perm_anovan[rep][j] = F_vec_perm_anovan[j]
                p_mat_perm_anovan[rep][j] = p_vec_perm_anovan[j]
        ref = F_mat_perm_anovan.max(1)
        Uncorrected_permuted_p = np.zeros(n_genes)
        for j in range(0, n_genes):
            sum = 0
            for k in range(0, n_rep):
                if(F_mat_perm_anovan[k][j] >= F_vec_ref_anovan[j]):
                    sum = sum+1
            Uncorrected_permuted_p[j] = sum/n_rep
            v = []
            for k in range(0, n_rep):
                if(ref[k] >= F_vec_ref_anovan[j]):
                    v.append(1)
                else:
                    v.append(0)
            FWE_corrected_p[j] = np.mean(v)
        print('FWE_corrected_p is', FWE_corrected_p)

        self.result = {}
        for j in range(0, n_genes):
            print(j,' ',geneIds[j],' ',uniqueId[j],' ',FWE_corrected_p[j])
            self.result[uniqueId[j]] = FWE_corrected_p[j]
        print(self.result.items())

    def pvalues(self):
        print(self.result.keys())
        return self.result;
