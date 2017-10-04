#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nibabel as nib
import requests
from future.standard_library import install_aliases
install_aliases()

dictionaryimages = {'FP1': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz',
'FP2': 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'}
MNI152 = True

class jubrain:
    def probability_map(regionname, coordspace):
        url = dictionaryimages[regionname]
        r = requests.get(url, verify=False)
        filename = 'output'+regionname+'.nii.gz'
        with open(filename, 'wb') as f:
            f.write(r.content)
        return nib.load(filename)
