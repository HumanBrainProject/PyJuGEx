#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyjugex
import nibabel as nib

refresh_cache = False
gene_cache = '/home/hbhattacharya/.pyjugex/cache/'
#CHECK THIS PART WITH TIMO
parcels = dict()
parcels['l0'] = nib.load('files/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz')
parcels['l1'] = nib.load('files/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz')
# specify candidate genes and trigger the workflow
analysis = pyjugex.DifferentialGeneExpression()
genelist = pyjugex.readCSVFile('files/MDD_Gene_List.csv') #CHECK THIS OUT
if refresh_cache is True:
    analysis.retrieve_gene_data(genelist, gene_cache)
analysis.set_candidate_genes(genelist, gene_cache)
analysis.set_coordinates_region(parcels['l0'], 0)
analysis.set_coordinates_region(parcels['l1'], 1)
analysis.run()
result = analysis.pvalues()
print([id for id in genelist['entrez_id'] if result[id] < .05])
