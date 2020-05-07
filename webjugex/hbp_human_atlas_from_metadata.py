# -*- coding: utf-8 -*-

# Copyright 2020 Forschungszentrum Jülich
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
