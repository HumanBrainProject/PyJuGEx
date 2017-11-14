#!/usr/bin/env python
# -*- coding: utf-8 -*-
import querydb

q = querydb.querydb()
#q.printrois()
id = q.getidfromname('Area 3a (PostCG)')
print(id)
pmapurl = q.getpmapurlfromid(id)
print(pmapurl)
