# -*- coding: utf-8 -*-
import querydb
import nibabel as nib
import requests

class jubrain:
    def __init__(self):
        self.q = querydb.querydb()

    def probability_map(self, regionname, coordspace):
        if coordspace != 'MNI152':
            print('Only MNI152 template space is supported')
            exit()
        id = self.q.getidfromname(regionname)
        if not id or id[0] is None:
            print('Regionname is not present in the database')
            exit()
        pmapurl = self.q.getpmapurlfromid(id[0])
        last = pmapurl.split('.')[-1]
        try:
            r = requests.get(pmapurl, verify=False)
        except requests.exceptions.RequestException as e:
            print(e)
            exit()
        if last == 'nii':
            filename = 'output'+regionname+'.'+last
        elif last == 'gz':
            filename = 'output'+regionname+'.nii.gz'
        else:
            print('Not an acceptable file format for rois')
            exit()
        with open(filename, 'wb') as f:
            f.write(r.content)
        return nib.load(filename)
