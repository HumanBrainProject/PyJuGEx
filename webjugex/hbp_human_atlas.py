#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Inpython2
import nibabel as nib
import requests
import tempfile
import os
import logging

dictionaryimages = {}

dictionaryimages['FP1'] = 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10m_l_N10_nlin2Stdicbm152casym.nii.gz'
dictionaryimages['FP2'] = 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/ba10p_l_N10_nlin2Stdicbm152casym.nii.gz'

MNI152 = True

class jubrain:        
    @classmethod
    def probability_map(cls, url, regionname, coordspace):
        if coordspace is False:
            raise ValueError('Only MNI152 template space is supported')
        last = url.split('.')[-1]
        try:
            r = requests.get(url, verify=False)
        except requests.HTTPError as e:
            logging.basicConfig(level=logging.INFO)
            logging.getLogger(__name__).error(e)
            raise
        if last not in ('nii', 'gz'):
            raise OSError('Not an acceptable file format for rois')
        fp, fp_name = tempfile.mkstemp(suffix='.'+last if last == 'nii' else '.nii.'+last)
        os.write(fp, r.content)
        img_array = nib.load(fp_name)
        os.close(fp)
        return img_array

    @classmethod
    def probability_map_v2(cls, url, regionname, body, coordspace):

        # v2 uses post instead of get, and expects optionally a body to be passed as body 
        if coordspace is False:
            raise ValueError('Only MNI152 template space is supported')
        try:
            r = requests.post(url, json=body)
        except requests.HTTPError as e:
            logging.basicConfig(level=logging.INFO)
            logging.getLogger(__name__).error(e)
            raise
        if r.status_code >= 400:
            raise ValueError(r.status_code)
        # no need to specify correct extension
        fp, fp_name = tempfile.mkstemp(suffix=".nii")
        os.write(fp, r.content)
        img_array = nib.load(fp_name)
        os.close(fp)
        return img_array

