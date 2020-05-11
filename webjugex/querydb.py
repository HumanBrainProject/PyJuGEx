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

#import pyArango
from pyArango.connection import *
import json
import sys

class querydb:
    def __init__(self):
        with open('config.json', 'r') as fp:
            configvals = json.load(fp)
        self.conn = Connection(arangoURL=configvals['arangoURL'], username=configvals['username'], password=configvals['password'])
        self.db = self.conn["Metadata"]
        self.rois = self.db["rois"]
        with open('files/filteredJuBrainJsonmod.json') as file:
            self.metadata = json.load(file)
        self.idpmapdict = {}
        self.createidpmap(self.metadata[0])

    def createidpmap(self, d):
        if isinstance(d, dict):
            if 'PMapURL' in d:
                if 'ontologyMetadata' in d and '_Id' in d['ontologyMetadata']:
                    self.idpmapdict[d['ontologyMetadata']['_Id']] = d['PMapURL']
                else:
                    self.idpmapdict[d['name']] = d['PMapURL']
                self.createidpmap(d['children'])
            if 'PMapURL' not in d:
                self.createidpmap(d['children'])
        elif isinstance(d, list):
            for k in d:
                self.createidpmap(k)
        else:
            print(type(d),' is not supported ')
            exit()

    def getidfromname(self, name):
        return [item['_Id'] for item in self.rois.fetchAll() if item['Display name'] == name]
        print('Roi name not found in the database')
        return

    '''
    def getpmap(self, id):
        for item in self.rois.fetchAll():
            if item['_Id'] is id:
                print(item['PmapURL']) #Ideally this should return the pmaps once they are part of the database
    '''

    def getpmapurlfromid(self, id):
        return self.idpmapdict[id]
        print('Not a valid id')
        return

    def printrois(self):
        if sys.version_info[0] < 3:
            for key, values in self.idpmapdict.iteritems():
                print(key,' values = ',values)
        else:
            for key, values in self.idpmapdict.items():
                print(key,' values = ',values)
