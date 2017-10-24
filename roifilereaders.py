#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
with open('roinames.txt') as f:
    content = f.readlines()
content = [x.strip() for x in content]
roiurls = {}
url = 'https://hbp-unic.fz-juelich.de:7112/UFTP/rest/access/JUDAC/2f054eee-7fa5-4ed3-b046-2ddc1315fe9c/'
for c in content:
    names = c.split('_')
    roiname = names[-1].split('.')[0]
    if len(names) > 2:
        roiname = names[-2] + '_' + names[-1].split('.')[0]
    else:
        roiname = names[-1].split('.')[0]
    roiurls[roiname] = url + c
with open('roidict.txt', 'w') as f:
    f.write(json.dumps(roiurls))
