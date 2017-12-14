#!/usr/bin/env python
# -*- coding: utf-8 -*-
import hbp_human_atlas as atlas
import analysispyjugex
#import hbp_human_atlas_from_metadata

genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']

roi1 = {}
roi2 = {}
roi1['name'] = 'FP1'
roi1['data'] = atlas.jubrain.probability_map('FP1', atlas.MNI152)
roi2['name'] = 'FP2'
roi2['data'] = atlas.jubrain.probability_map('FP2', atlas.MNI152)

jugex = analysispyjugex.Analysis(gene_cache_dir='.pyjugex', verbose=True)
result = jugex.DifferentialAnalysis(genelist, roi1, roi2)
print(result)
