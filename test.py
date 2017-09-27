#!/usr/bin/env python
# -*- coding: utf-8 -*-
import requests
import nibabel as nib
import os
dictionary_images = {'FP1': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz',
'FP2': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'}
url = dictionary_images['FP2']
r = requests.get(url, verify=False)
with open('output.nii.gz', 'wb') as f:
    f.write(r.content)
file1 = nib.load('output.nii.gz')
print(file1.header)

url = dictionary_images['FP2']
r = requests.get(url, verify=False)
with open('output.nii.gz', 'wb') as f:
    f.write(r.content)
file1 = nib.load('output.nii.gz')
print(file1.header)
