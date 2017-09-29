#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyjugex
import hbp_human_atlas as atlas
import nibabel as nib

'''
with analysis.py
In pyjugex.py
uncomment "import analysis as mymodule" and comment "import analysispyjugex as mymodule"
dictionaryimages = {'FP1': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz',
'FP2': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'}
parcels = {}
parcels['l0'] = nib.load(pyjugex.ReadParcellations(dictionaryimages, 'FP1'))
parcels['l1'] = nib.load(pyjugex.ReadParcellations(dictionaryimages, 'FP2'))
analysis = pyjugex.DifferentialGeneExpression(gene_cache = '~/.pyjugex')
genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']
'''
#WITH ANALYSISPYJUGEX.PY
genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']
roi1 = atlas.jubrain.probability_map('FP1', atlas.MNI152)
roi2 = atlas.jubrain.probability_map('FP2', atlas.MNI152)
jugex = pyjugex.PyJugex(cache="~/.pyjugex")
result = jugex.DifferentialAnalysis(genelist, roi1, roi2)
print([id for id in result if result[id] < .05])


'''
ISSUES:
with analysis.py
Run with genelist = ['GAD2', 'MAOA']
Then run with genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']
This gives different result than running without any cache and
genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']
'''
