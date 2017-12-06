#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Inpython2
import nibabel as nib
import requests
import tempfile
import os
dictionaryimages = {}
'''
with open('roinames.txt') as f:
    content = f.readlines()
content = [x.strip() for x in content]
url = 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/'
for c in content:
    names = c.split('_')
    roiname = names[-1].split('.')[0]
    if len(names) > 2:
        roiname = names[-2] + '_' + names[-1].split('.')[0]
    else:
        roiname = names[-1].split('.')[0]
    dictionaryimages[roiname] = url + c
'''

dictionaryimages['FP1'] = 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz'
dictionaryimages['FP2'] = 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'

MNI152 = True

class jubrain:        
    @classmethod
    def probability_map(cls, regionname, coordspace):
        if coordspace is False:
            print('Only MNI152 template space is supported')
            exit()
        if regionname not in dictionaryimages:
            print('Regionname',regionname,'is not present in the database')
            exit()
        url = dictionaryimages[regionname]
        last = url.split('.')[-1]
        try:
            r = requests.get(url, verify=False)
        except requests.HTTPError as e:
            print(e)
            raise
        try:
            if last not in ('nii', 'gz'):
                raise OSError
        except OSError:
            print('Not an acceptable file format for rois')
            raise
        fp, fp_name = tempfile.mkstemp(suffix='.'+last if last == 'nii' else '.nii.'+last)
        os.write(fp, r.content)
        img_array = nib.load(fp_name)
        os.close(fp)
        return img_array

