#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2020 Human Brain Project/EBRAINS
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

import csv
import sys
"""
Read a list of gene symbols from a .txt file  and saves the list into a .csv file, to be used by the import functionality of webjugex ui
Args:
      sys.argv[1] (str): Name of .txt file which contains gene symbols, each gene symbol in one line.
      sys.argv[2] (str): name of .csv file which will contain the gene symbols in a csv format
Usage:
    python makegenesreadyforimport.py genelistwebjugex.txt genelistwebjugex.csv
"""
f = open(sys.argv[2], 'w')
w = csv.writer(f, delimiter = ',')
genelist = open(sys.argv[1], 'r').read()
genelist = genelist.split('\n')
w.writerows([x.split(',') for x in genelist])
f.close()
