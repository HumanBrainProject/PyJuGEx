#!/usr/bin/env python
# -*- coding: utf-8 -*-
import querydb

q = querydb.querydb()
#q.printrois()
id = q.getidfromname('Area 5M (SPL)')
print(id)
pmapurl = q.getpmapurlfromid(id)
print(pmapurl)
