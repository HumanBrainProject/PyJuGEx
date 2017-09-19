#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyjugex
import nibabel as nib

dictionaryimages = {'FP1': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz',
'FP2': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'}
parcels = {}
parcels['l0'] = nib.load(pyjugex.ReadParcellations(dictionaryimages, 'FP1'))
parcels['l1'] = nib.load(pyjugex.ReadParcellations(dictionaryimages, 'FP2'))
analysis = pyjugex.DifferentialGeneExpression(gene_cache = '/home/hbhattacharya/.pyjugex/cache246/')
#analysis = pyjugex.DifferentialGeneExpression(gene_cache = None)
genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']
analysis.set_candidate_genes(genelist)
analysis.set_ROI_MNI152(parcels['l0'], 0)
analysis.set_ROI_MNI152(parcels['l1'], 1)
analysis.run()
result = analysis.pvalues()
print([id for id in result if result[id] < .05])
