#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyjugex
import nibabel as nib

refresh_cache = False
dictionaryimages = {'FP1': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz',
'FP2': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'}
parcels = {}
parcels['l0'] = nib.load(pyjugex.ReadParcellations(dictionaryimages, 'FP1'))
parcels['l1'] = nib.load(pyjugex.ReadParcellations(dictionaryimages, 'FP2'))
# specify candidate genes and trigger the workflow
analysis = pyjugex.DifferentialGeneExpression(gene_cache = '/home/hbhattacharya/.pyjugex/cache/')
data = pyjugex.readCSVFile('files/MDD_Gene_List.csv') #CHECK THIS OUT
#gene_symbol = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']
#entrez_id = [150, 553, 1129, 1268, 1385, 1392, 1394, 1395, 2572, 3350, 3351, 3352, 3356, 3359, 3361, 4128, 5136, 6530, 6532, 6750, 6863, 7166, 9248, 23316, 121278]
genelist = {}
genelist['gene_symbols'] = data['gene_symbol']
genelist['entrez_id'] = data['entrez_id']
genelist['probe_ids'] = data['probe_id']
analysis.set_candidate_genes(genelist)

analysis.set_ROI_MNI152(parcels['l0'], 0)
analysis.set_ROI_MNI152(parcels['l1'], 1)

analysis.run()

result = analysis.pvalues()
print([id for id in result if result[id] < .05])
