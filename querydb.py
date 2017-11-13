# -*- coding: utf-8 -*-
#import pyArango
from pyArango.connection import *
import json
import sys

class querydb:
    def __init__(self):
        self.conn = Connection(arangoURL='http://cudaws02.ime.kfa-juelich.de:8529', username="haimasree", password="haimasree123")
        self.db = self.conn["Metadata"]
        self.rois = self.db["rois"]
        with open('files/filteredJuBrainJson.json') as file:
            self.metadata = json.load(file)
        self.idpmapdict = {}



    def fun(self, d):
        if isinstance(d, dict):
            if 'PMapURL' in d:
                self.idpmapdict[d['name']] = d['PMapURL']
                self.fun(d['children'])
            if 'PMapURL' not in d:
                self.fun(d['children'])
        elif isinstance(d, list):
            for k in d:
                self.fun(k)
        else:
            print(type(d))
            exit()

    def printrois(self):
        if sys.version_info[0] < 3:
            for key, value in self.metadata[0].iteritems():
                print(key,' ',value)
        else:
            self.fun(self.metadata[0])
            for key, values in self.idpmapdict.items():
                print(key,' values = ',values)
