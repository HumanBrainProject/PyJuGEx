import numpy as np
import nibabel as nib
import urllib.request
import json
from numpy.linalg import inv
import csv

from numpy import *
import scipy as sp
import scipy.stats.mstats
import rpy2.robjects as robjects
import rpy2.rlike.container as rlc
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from numpy.random import normal


def buildSpecimenFactors():
    url = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"
    specimenFactors = dict()
    response = urllib.request.urlopen(url).read().decode('utf8')
    text = json.loads(response)
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
        '''
        data = dict();
        data['id'] = res[i]['id']
        data['name'] = res[i]['name']
        data['race'] = res[i]['race_only']
        data['gender'] = res[i]['sex']
        data['age'] = res[i]['age']['days']/365
        specimenFactors.append(data)
        '''
    return specimenFactors;

def transformSamples(samples, T):
    nsamples = len(samples)
    coords = []
    #coords = np.zeros((nsamples, 3))
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
#        cx = int(T00*x + T01*y + T02*z + T03)
#        cy = int(T10*x + T11*y + T12*z + T13)
#        cz = int(T20*x + T21*y + T22*z + T23)
        #coord = [round(cx), round(cy), round(cz)]
        coord = [cx, cy, cz]
        coords.append(coord)
    print(len(coords))
    #np.savetxt('Output.txt', coords, fmt='%-7.3f')
    return coords

def readCSVFile(filename):
    rows = dict();
    rows['probe_id'] = []
    rows['gene_symbol'] = []
    rows['entrez_id'] = []
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            rows['probe_id'].append(row['probe_id'])
            rows['gene_symbol'].append(row['gene_symbol'])
            rows['entrez_id'].append(int(row['entrez_id']))
            #print(row['probe_id'], row['gene_symbol'], row['entrez_id'])
    return rows

def extractExpLevel(apiData, vois, specimenInfo, Mni, mapThreshold, searchMode):
    #
    #maps = struct('pmap',pmap,...
    #    'name',name,...
    #   'bar_plot_color',bar_plot_color,...
    #    'validated_zscores','',...
    #    'correlation_coeff','',...
    #    'specimen','');

    #main_r= struct('pmap','',...
    #    'name','',...
     #   'bar_plot_color','',...
     #   'validated_zscores','',...
     #   'correlation_coeff','',...
     #   'specimen','',...
     #   'data2plot','');
    main_r = []
    for i in range(0, len(specimenInfo)):
        revisedApiDataCombo = dict()
        revisedApiData = expressionSpmCorrelation(apiData[i], vois['img1'], specimenInfo[i], Mni, mapThreshold, searchMode)
        revisedApiDataCombo['explevels'] = revisedApiData['explevels'][:]
        revisedApiDataCombo['zscores'] = revisedApiData['zscores'][:]
        revisedApiDataCombo['coords'] = revisedApiData['coords'][:]
        revisedApiDataCombo['samples'] = revisedApiData['samples'][:]
        revisedApiDataCombo['probes'] = revisedApiData['probes'][:]
        revisedApiDataCombo['specimen'] = revisedApiData['specimen']
        revisedApiDataCombo['name'] = 'img1'
        print('extractexplevel img1: ',revisedApiDataCombo['specimen'],' ',len(revisedApiDataCombo['coords']))
        main_r.append(revisedApiDataCombo)
        revisedApiDataCombo = dict()
        revisedApiData = expressionSpmCorrelation(apiData[i], vois['img2'], specimenInfo[i], Mni, mapThreshold, searchMode)
        revisedApiDataCombo['explevels'] = revisedApiData['explevels'][:]
        revisedApiDataCombo['zscores'] = revisedApiData['zscores'][:]
        revisedApiDataCombo['coords'] = revisedApiData['coords'][:]
        revisedApiDataCombo['samples'] = revisedApiData['samples'][:]
        revisedApiDataCombo['probes'] = revisedApiData['probes'][:]
        revisedApiDataCombo['specimen'] = revisedApiData['specimen']
        print('extractexplevel img2: ',revisedApiDataCombo['specimen'],' ',len(revisedApiDataCombo['coords']))
        revisedApiDataCombo['name'] = 'img2'
        main_r.append(revisedApiDataCombo)
    return main_r;


def switch2gensymbol(entrez_id, combined_zscores, area1Len, area2Len):
    #u, ind = np.unique(entrez_id, return_index=True)
    #unique_entrez_id = u[np.argsort(ind)]
    unique_entrez_id = np.unique(entrez_id)
    #print(entrez_id[index] for index in idx)
    #unique_entrez_id = entrez_id[np.sort(idx)]
    print('length of entrez id is ', len(unique_entrez_id))
    print(unique_entrez_id)
    inds =  [[] for _ in range(len(unique_entrez_id))]
    for i in range(0, len(unique_entrez_id)):
        for j in range (0, len(entrez_id)):
            if unique_entrez_id[i] == entrez_id[j]:
                inds[i].append(j);
    #for i in range(0, len(unique_entrez_id)):
    #    print(len(inds[i]))
    winsorzed_mean_zscores = np.zeros((len(combined_zscores), len(unique_entrez_id)))
    for i in range (0, len(unique_entrez_id)):
        tmp = [[] for _ in range(len(combined_zscores))]
        for j in range(0, len(combined_zscores)):
            for k in range (0, len(inds[i])):
                tmp[j].append(combined_zscores[j][inds[i][k]]);
            ncols = len(tmp[j]);
        for j in range(0, len(combined_zscores)):
            #print(scipy.stats.mstats.winsorize(tmp[j], limits=0.1))
            winsorzed_mean_zscores[j][i] = np.mean(scipy.stats.mstats.winsorize(tmp[j], limits=0.1))#maybe add a special case for one element in tmp
    #np.savetxt('switch2gensymbol.txt', winsorzed_mean_zscores, fmt='%-7.3f ')
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


def performAnova(main_r, searchMode):
    #Get data ready for anova
    #factors={factor_area factor_specimen factor_age_numeric factor_race}; varnames = {'Area';'Specimen';'Age';'Race'};
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
    #CHECK WHETHER TO USE APPEND OR COPY LIST
    print(" ",len(main_r)," ",main_r[0]['name']," ",main_r[1]['name'])
    #print("at this point ",len(main_r)," ",len(main_r[0]['coords'])," ",len(main_r[1]['coords'])," ",len(main_r[0]['zscores'][0])," ",len(main_r[1]['zscores'][0])," ",main_r[0]['specimen']," ",main_r[1]['specimen']," ",main_r[0]['name']," ",main_r[1]['name'])
    for i in range(0, len(main_r)):
        if(main_r[i]['name'] == 'img1'):
            for j in range(0, len(main_r[i]['zscores'])):
                area1_zscores.append(main_r[i]['zscores'][j])
                #combined_zscores.append(main_r[i]['zscores'][j])
                area1_specimen.append(main_r[i]['specimen'])
                area1_area.append(main_r[i]['name'])
                area1_koords.append(main_r[i]['coords'])
        elif(main_r[i]['name'] == 'img2'):
            for j in range(0, len(main_r[i]['zscores'])):
                area2_zscores.append(main_r[i]['zscores'][j])
                #combined_zscores.append(main_r[i]['zscores'][j])
                area2_specimen.append(main_r[i]['specimen'])
                area2_area.append(main_r[i]['name'])
                area2_koords.append(main_r[i]['coords'])
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
    #create additional factors (age, gender, race) based on specimen_info
    specimenFactors = buildSpecimenFactors()
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

    if searchMode == 1:
        geneList = readCSVFile('./MDD_Gene_List.csv')
        print(len(geneList))        
        allProbeData = switch2gensymbol(geneList['entrez_id'],combined_zscores, len(area1_zscores), len(area2_zscores))
        combined_zscores = np.zeros((len(allProbeData['combined_zscores']), len(allProbeData['combined_zscores'])))
        combined_zscores = np.copy(allProbeData['combined_zscores'])
        area1_zscores = np.zeros((len(allProbeData['area1_zscores']), len(allProbeData['area1_zscores'])))
        area1_zscores = np.copy(allProbeData['area1_zscores'])
        area2_zscores = np.zeros((len(allProbeData['area2_zscores']), len(allProbeData['area2_zscores'])))
        area2_zscores = np.copy(allProbeData['area2_zscores'])
        uniqueId = np.copy(allProbeData['uniqueId'])
        geneIds = []
        for i in range(0, len(uniqueId)):
            for j in range(0, len(geneList['entrez_id'])):
                if(geneList['entrez_id'][j] == uniqueId[i]):
                    geneIds.append(geneList['gene_symbol'][j])
                    break
        print(geneIds)
        '''
        np.savetxt('switch2gensymbol1.txt', combined_zscores, fmt='%-7.3f ')
        np.savetxt('switch2gensymbol3.txt', area1_zscores, fmt='%-7.3f ')
        np.savetxt('switch2gensymbol4.txt', area2_zscores, fmt='%-7.3f ')
        '''
        n_genes = len(combined_zscores[0]) #SHOULD NOT THIS BE 285???
        print(n_genes)
        Reference_Anovan_p = np.zeros(n_genes)
        Reference_Anovan_eta2 = np.zeros(n_genes)
        Reference_Anovan_CI_l = np.zeros(n_genes)
        Reference_Anovan_CI_h = np.zeros(n_genes)
        Reference_Anovan_diff_mean = np.zeros(n_genes)
        F_vec_ref_anovan = np.zeros(n_genes)
    elif searchMode == 2:
        n_genes = len(combined_zscores[0]) #SHOULD NOT THIS BE 285???
        Reference_Anovan_p = np.zeros(n_genes)
        Reference_Anovan_eta2 = np.zeros(n_genes)
        Reference_Anovan_CI_l = np.zeros(n_genes)
        Reference_Anovan_CI_h = np.zeros(n_genes)
        Reference_Anovan_diff_mean = np.zeros(n_genes)
        F_vec_ref_anovan = np.zeros(n_genes)
        #print(robjects.r['pi'])
        #print(combined_zscores[2]," ",len(combined_zscores))
    #print('some info')
    #np.savetxt('zscores.txt', combined_zscores, fmt='%f ')
    #for i in range(0, n_genes):
        #z = robjects.FloatVector(combined_zscores[:,i])
        #print(z)

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
        #print('var1')
        #print(var1)
        #print('var2')
        #print(var2)
        #print('variances are ',np.var(var1), ' ', np.var(var2))
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
    '''
    np.savetxt('F_vec_ref_anovan.txt', F_vec_ref_anovan, fmt='%f')
    np.savetxt('Reference_Anovan_eta2.txt', Reference_Anovan_eta2, fmt='%f')
    np.savetxt('Reference_Anovan_CI_l.txt', Reference_Anovan_CI_l, fmt='%f')
    np.savetxt('Reference_Anovan_CI_h.txt', Reference_Anovan_CI_h, fmt='%f')
    np.savetxt('Reference_Anovan_diff_mean.txt', Reference_Anovan_diff_mean, fmt='%f')
    np.savetxt('Reference_Anovan_p.txt', Reference_Anovan_p, fmt='%f')
    '''
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
            #f = ['img1', 'img2', 'img2', 'img2', 'img1', 'img2', 'img1', 'img2', 'img1', 'img1', 'img2', 'img2', 'img2']
            #v = robjects.StrVector(f)
            '''
            od = rlc.OrdDict([('Area', robjects.StrVector(factor_area)),
                              ('Specimen', robjects.StrVector(factor_specimen)),
                              ('Age', robjects.IntVector(factor_age_numeric)),
                              ('Race', robjects.StrVector(factor_race)),
                              ('Zscores', robjects.FloatVector(combined_zscores[:,i]))])

            od = rlc.OrdDict([('Area', v),
            '''
            od = rlc.OrdDict([('Area', f),
                              ('Specimen', robjects.StrVector(factor_specimen)),
                              ('Age', robjects.IntVector(factor_age_numeric)),
                              ('Race', robjects.StrVector(factor_race)),
                              ('Zscores', robjects.FloatVector(combined_zscores[:,j]))])
            dataf = robjects.DataFrame(od)
    #        print(dataf)            #f = robjects.Formula('Zscores~Area*Specimen*Age*Race')
            #r['options'](contrasts=r['c']("contr.sum","contr.poly"))
            #f = robjects.Formula('Zscores~Specimen+Area+Age+Race')
            #a = r['aov'](f, data = dataf,  type=3)
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
    res = []
    for j in range(0, n_genes):
        #if(FWE_corrected_p[j] < 0.05):
            #print(j,' ',geneIds[j],' ',FWE_corrected_p[j])
            print(j,' ',geneIds[j],' ',uniqueId[j],' ',FWE_corrected_p[j])
            temp = dict()
            temp['id'] = geneIds[j]
            temp['p'] = FWE_corrected_p[j]
            res.append(temp)
    resJson = {}
    resJson["genes"] = []
    for i in range(0, n_genes):
        if(FWE_corrected_p[i] < 0.05):
            resJson["genes"].append({
                'name':geneIds[i],
                'pval':FWE_corrected_p[i]
            })
    print('resJson',resJson)
    jsonStr = json.dumps(resJson)
    print(jsonStr)
    return jsonStr

def expressionSpmCorrelation(apiData, img1, specimen, Mni, mapThreshold, searchMode):
    revisedApiData = dict()
    revisedApiData['explevels'] = []
    revisedApiData['zscores'] = []
    revisedApiData['coords'] = []
    revisedApiData['samples'] = []
    revisedApiData['probes'] = []
    revisedApiData['specimen'] = []
    Mni = specimen['alignment3d']
    '''
    Mni = np.zeros((4, 4))
    Mni[0][0] = -1
    Mni[0][3] = 91
    Mni[1][2] = -1
    Mni[1][3] = 91
    Mni[2][1] = -1
    Mni[2][3] = 109
    Mni[3][3] = 1
    '''
    #Create a numpy matrix and double check Mni = []
    if searchMode == 2:
        zscores = apiData['zscores']
        mni = img1.affine
        invMni = inv(mni)
        T = np.dot(invMni, specimen['alignment3d'])
        coords = transformSamples(apiData['samples'], T) #I COULD NOT FIND ROUNDING FUNCTION
        revisedApiData['explevels'] = apiData['explevels'][:]
        revisedApiData['zscores'] = apiData['zscores'][:]
        revisedApiData['coords'] = coords[:]
        revisedApiData['samples'] = apiData['samples'][:]
        revisedApiData['probes'] = apiData['probes'][:]
        revisedApiData['specimen'] = specimen['name']
        #coords = int32(round(coords)); #I COULD NOT FIND ROUNDING FUNCTION      
    elif searchMode == 1:
        #print(Mni)
        data1 = img1.get_data()
        imgMni = img1.affine
        #print('spmmni')
        #print(imgMni)
        #print('invmni')
        invimgMni = inv(imgMni)
        #print(invimgMni)
        Mni = specimen['alignment3d']
        #print('MNI')
        #print(Mni)
        T = np.dot(invimgMni, Mni)
        #print('aibstospm')
        #print(T)
        coords = transformSamples(apiData['samples'], T)
        #np.savetxt('oldcoord.txt', coords, fmt='%-7.3f')
        for i in range(0, len(coords)):
            coords[i] = np.rint(coords[i])
        #np.savetxt('roundedcoord.txt', coords, fmt='%-7.3f')
        #coords = int32(round(coords));
        print("After transform samples ")
        #print(len(apiData['explevels']))
        for i in range(0, len(coords)):
            coord = coords[i]
            #print(type(coord[0]))
            sum = 0
            for j in range(0, len(coord)) :
                if(coord[j] > 0):
                    sum += 1
            if(sum != 3 or data1[int(coord[0]),int(coord[1]),int(coord[2])] <= mapThreshold/10 or data1[int(coord[0]),int(coord[1]),int(coord[2])] == 0):
                #print(coord)
                coords[i][0] = -1 #IS IT 1 or -1
                coords[i][1] = -1
                coords[i][2] = -1
        #np.savetxt('newcoord.txt', coords, fmt='%-7.3f')
        for i in range(0, len(coords)):
            coord = coords[i]
            sum = 0
            for j in range(0, len(coord)) :
                if(coord[j] > 0) :
                    sum += 1
            if(sum == 3):
                revisedApiData['explevels'].append(apiData['explevels'][i])
                revisedApiData['zscores'].append(apiData['zscores'][i])
                revisedApiData['coords'].append(coord)
        revisedApiData['samples'] = apiData['samples'][:]
        revisedApiData['probes'] = apiData['probes'][:]
        revisedApiData['specimen'] = specimen['name']
        for i in range(0, len(revisedApiData['coords'])):
            print(revisedApiData['coords'][i])
            print(revisedApiData['zscores'][i].shape)
        #np.savetxt('newzscores.txt', revisedApiData['zscores'], fmt='%-7.3f')

        #np.savetxt('outputcoord.txt', revisedApiData['coords'], fmt='%-7.3f')
        #for i in range(0, len(revisedApiData['coords'])):
        #    text_file.write("%f " %revisedApiData['coords'][i])
        #text_file.close()

    return revisedApiData

def getSpecimenData(info):
    #info = msg[0]
    specimenD = dict()
    '''
    specimenD['cell_prep_sample_id'] = info.cell_prep_sample_id
    specimenD['cell_reporter_id'] = info.cell_reporter_id
    specimenD['data'] = info.data
    specimenD['donor_id'] = info.donor_id
    specimenD['ephys_result_id'] = info.ephys_result_id
    specimenD['external_specimen_name'] = info.external_specimen_name
    specimenD['failed_facet'] = info.failed_facet
    specimenD['hemisphere'] = info.hemisphere
    specimenD['id'] = info.id
    specimenD['is_cell_specimen'] = info.is_cell_specimen
    specimenD['is_ish'] = info.is_ish
    specimenD['name'] = info.name
    specimenD['parent_id'] = info.parent_id
    specimenD['parent_x_coord'] = info.parent_x_coord
    specimenD['parent_y_coord'] = info.parent_y_coord
    specimenD['parent_z_coord'] = info.parent_z_coord
    specimenD['rna_integrity_number'] = info.rna_integrity_number
    specimenD['specimen_id_path'] = info.specimen_id_path
    specimenD['sphinx_id'] = info.sphinx_id
    specimenD['structure_id'] = info.structure_id
    specimenD['tissue_ph'] = info.tissue_ph
    specimenD['treatment_id'] = info.treatment_id
    specimenD['weight'] = info.weight
    '''
    specimenD['name'] = info['name']
    x = info['alignment3d']
    alignment3dMat = np.zeros((4, 4))
    alignment3dMat[0][0] = x['tvr_00']
    alignment3dMat[0][1] = x['tvr_01']
    alignment3dMat[0][2] = x['tvr_02']
    alignment3dMat[0][3] = x['tvr_09']
    alignment3dMat[1][0] = x['tvr_03']
    alignment3dMat[1][1] = x['tvr_04']
    alignment3dMat[1][2] = x['tvr_05']
    alignment3dMat[1][3] = x['tvr_10']
    alignment3dMat[2][0] = x['tvr_06']
    alignment3dMat[2][1] = x['tvr_07']
    alignment3dMat[2][2] = x['tvr_08']
    alignment3dMat[2][3] = x['tvr_11']
    alignment3dMat[3][0] = 0
    alignment3dMat[3][1] = 0
    alignment3dMat[3][2] = 0
    alignment3dMat[3][3] = 1
    specimenD['alignment3d'] =  alignment3dMat
    return specimenD


def downloadSpecimens():
    specimens  = ["H0351.1015", "H0351.1012", "H0351.1016", "H0351.2001", "H0351.1009", "H0351.2002"];
    #specimens  = ["H0351.1015", "H0351.1012"];
    specimenInfo = [];
    for i in range(0, len(specimens)):
        url = "http://api.brain-map.org/api/v2/data/Specimen/query.json?criteria=[name$eq"
        url+="'"
        url += specimens[i]
        url += "']&include=alignment3d"
        print(url)
        response = urllib.request.urlopen(url).read().decode('utf8')
        text = json.loads(response)
        data = text['msg'][0]
        res = getSpecimenData(data);
        specimenInfo.append(res);
    print(len(specimenInfo))
    return specimenInfo
    '''
    for k in range(0, len(specimenInfo)):
         x = specimenInfo[k]['alignment3d']
         for i in range(0, 4):
            for j in range(0, 4):
              print(x[i][j])

    return specimenInfo
    '''
def queryAPI(donorId):
    url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in1059351,1059352,1059353,1058915,1058914,1058916,1029156,1029155,1029150,1029149,1029146,1029142,1029138,1029136,1029135,1029134,1029133,1029126,1029125,1029124,1029123,1029122,1029121,1025133,1023492,1029129,1029128,1029127,1028516,1028515,1028514,1028513,1028512,1028511,1028510,1028509,1028508,1028507,1028506,1028505,1028504,1028503,1028502,1028484,1028483,1028478,1028477,1028474,1028473,1028472,1028471,1028470,1028469,1028466,1028465,1028464,1028462,1028461,1028458,1028457,1028456,1028455,1028454,1028453,1028452,1028451,1028450,1028448,1028447,1028446,1028445,1028444,1028442,1028441,1028440,1028439,1028438,1028437,1028434,1028433,1028432,1028431,1028430,1028429,1028428,1028427,1028426,1028425,1028424,1015329,1028436,1028435,1028501,1028500,1028499,1028498,1028497,1028494,1028493,1028492,1028491,1028490,1028489,1028488,1028487,1028486,1028485,1028548,1028547,1028542,1028536,1028535,1028534,1028527,1028526,1028518,1028517,1057976,1057975,1057974,1059659,1057967,1057966,1057965,1019360,1019310,1013314,1011485,1057962,1057961,1057960,1028346,1028345,1028344,1056548,1056547,1056546,1056545,1056544,1013593,1013489,1013359,1013692,1012504,1012098,1011320,1011137,1010770,1010725,1010443,1010639,1056550,1048801,1048800,1024804,1024803,1024812,1055392,1055381,1055378,1055374,1055373,1055372,1055361,1055357,1055356,1055355,1055354,1055347,1055346,1055345,1055342,1055377,1024617,1024606,1024594,1024593,1024592,1021988,1021983,1021982,1021981,1021979,1021975,1021974,1021973,1021972,1021970,1021969,1021968,1054369,1054368,1054364,1054362,1054360,1054359,1054358,1054357,1054353,1054352,1054350,1054349,1054348,1054347,1054345,1054344,1054343,1054342,1054341,1054340,1054339,1054338,1054337,1054336,1054335,1054334,1054380,1054379,1054375,1054372,1054370,1053289,1053288,1053287,1063261,1023147,1023146,1051465,1051464,1051461,1051460,1051459,1051458,1051457,1051489,1051483,1051472,1051466,1025754,1025749,1025748,1025747,1051004,1029261,1029258,1029257,1029256,1029255,1029254,1029253,1029252,1029251,1029250,1029249,1028377,1050663,1050662,1050661,1035832,1035831,1028612,1028611,1028610,1028609,1028608,1028607,1028606,1028605,1028604,1028603,1028602,1028601,1028634,1028633,1028631,1028627,1028626,1028625,1028620,1028619,1028618,1028617,1028615,1028614,1028613][donors$eq"
    url += donorId;
    url += "]"
    print(url)
    #url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in1059351,1059352,1059353,1058915,1058914,1058916,1029156,1029155,1029150,1029149,1029146,1029142,1029138,1029136,1029135,1029134,1029133,1029126,1029125,1029124,1029123,1029122,1029121,1025133,1023492,1029129,1029128,1029127,1028516,1028515,1028514,1028513,1028512,1028511,1028510,1028509,1028508,1028507,1028506,1028505,1028504,1028503,1028502,1028484,1028483,1028478,1028477,1028474,1028473,1028472,1028471,1028470,1028469,1028466,1028465,1028464,1028462,1028461,1028458,1028457,1028456,1028455,1028454,1028453,1028452,1028451,1028450,1028448,1028447,1028446,1028445,1028444,1028442,1028441,1028440,1028439,1028438,1028437,1028434,1028433,1028432,1028431,1028430,1028429,1028428,1028427,1028426,1028425,1028424,1015329,1028436,1028435,1028501,1028500,1028499,1028498,1028497,1028494,1028493,1028492,1028491,1028490,1028489,1028488,1028487,1028486,1028485,1028548,1028547,1028542,1028536,1028535,1028534,1028527,1028526,1028518,1028517,1057976,1057975,1057974,1059659,1057967,1057966,1057965,1019360,1019310,1013314,1011485,1057962,1057961,1057960,1028346,1028345,1028344,1056548,1056547,1056546,1056545,1056544,1013593,1013489,1013359,1013692,1012504,1012098,1011320,1011137,1010770,1010725,1010443,1010639,1056550,1048801,1048800,1024804,1024803,1024812,1055392,1055381,1055378,1055374,1055373,1055372,1055361,1055357,1055356,1055355,1055354,1055347,1055346,1055345,1055342,1055377,1024617,1024606,1024594,1024593,1024592,1021988,1021983,1021982,1021981,1021979,1021975,1021974,1021973,1021972,1021970,1021969,1021968,1054369,1054368,1054364,1054362,1054360,1054359,1054358,1054357,1054353,1054352,1054350,1054349,1054348,1054347,1054345,1054344,1054343,1054342,1054341,1054340,1054339,1054338,1054337,1054336,1054335,1054334,1054380,1054379,1054375,1054372,1054370,1053289,1053288,1053287,1063261,1023147,1023146,1051465,1051464,1051461,1051460,1051459,1051458,1051457,1051489,1051483,1051472,1051466,1025754,1025749,1025748,1025747,1051004,1029261,1029258,1029257,1029256,1029255,1029254,1029253,1029252,1029251,1029250,1029249,1028377,1050663,1050662,1050661,1035832,1035831,1028612,1028611,1028610,1028609,1028608,1028607,1028606,1028605,1028604,1028603,1028602,1028601,1028634,1028633,1028631,1028627,1028626,1028625,1028620,1028619,1028618,1028617,1028615,1028614,1028613][donors$eq15496,14380,15697,9861,12876,10021][structures$eq4009]"
    response = urllib.request.urlopen(url).read().decode('utf8')
    text = json.loads(response)
    data = text['msg']
    samples = []
    probes = []
    nsamples = len(data['samples'])
    nprobes = len(data['probes'])
    for i in range(0, nsamples):
        samples.append(data['samples'][i])
    for i in range(0, nprobes):
        probes.append(data['probes'][i]);
    #print(len(samples))
    #print(len(probes))
    explevels = np.zeros((nsamples, nprobes))
    zscores = np.zeros((nsamples, nprobes))
    for i in range(0, nprobes):
        for j in range(0, nsamples):
            explevels[j][i] = probes[i]['expression_level'][j]
            zscores[j][i] = probes[i]['z-score'][j]
    apiData = dict()
    apiData['samples'] = samples
    apiData['probes'] = probes
    apiData['explevels'] = explevels
    apiData['zscores'] = zscores
    #print("In queryApi")
    print(len(apiData['samples']), ' ', apiData['zscores'].shape, ' ', apiData['explevels'].shape, ' ', len(apiData['probes']))
    return apiData

def getAPIData(donorIds):
    apiData = []
    for i in range(0, len(donorIds)):
        res = queryAPI(donorIds[i])
        apiData.append(res)
#    return apiData;
    #print(len(apiData['samples']))
    #print(len(apiData['probes']))
    #print(apiData['zscores'][0][2])
    #print(apiData['explevels'][len(apiData['probes'])-1][len(apiData['samples'])-1])
    return apiData;

def readVOIs():
    vois = dict()
    example_filename1 = './ba10m_l_N10_nlin2Stdicbm152casym.nii.gz'
    vois['img1'] = nib.load(example_filename1)
    header1 = vois['img1'].header
    print(header1)
    example_filename2 = './ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'
    vois['img2'] = nib.load(example_filename2)
    header2 = vois['img2'].header
    print(header2)
    return vois

def readVOI(voiName):
    img = nib.load(voiName)
    return img

def performJugex():
    donorIds = ['15496','14380','15697','9861','12876','10021']
    #donorIds = ['15496', '14380']
    vois = readVOIs()
    apiData = getAPIData(donorIds)
    specimenInfo = downloadSpecimens()
    Mni = specimenInfo[0]['alignment3d'] #GET THE CORRECT VALUE HERE
    mapThreshold = 2
    searchMode = 1
    main_r = extractExpLevel(apiData, vois, specimenInfo, Mni, mapThreshold, searchMode)    
    res = performAnova(main_r, searchMode)
    return res
    '''
    res = {}
    res["genes"] = []
    res["genes"].append({
        'name':'abc',
        'pval':'def'
    })
    res["genes"].append({
        'name':'ghi',
        'pval':'jkl'
    })
    res["genes"].append({
        'name':'mno',
        'pval':'pqr'
    })
    jsonStr = json.dumps(res)
    print(jsonStr)
    return jsonStr
    '''
