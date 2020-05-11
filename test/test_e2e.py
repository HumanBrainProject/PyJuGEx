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

import sys
sys.path.append('..')

import os
nii_dir = os.path.dirname(os.path.abspath(__file__))

import pyjugex
import nibabel as nib

def test_e2e():

  gene_list=['MAOA','TAC1']

  nii1 = nib.load(os.path.join(nii_dir, 'data/hoc1_th_l.nii.gz'))
  nii2 = nib.load(os.path.join(nii_dir, 'data/hoc2_th_l.nii.gz'))
  analysis = pyjugex.analysis(
    gene_list=gene_list,
    roi1 = nii1,
    roi2 = nii2
  )

  filtered_coord = analysis.get_filtered_coord()
  assert(len(filtered_coord['roi1']) == 12)
  assert(len(filtered_coord['roi2']) == 11)

  analysis.run()

  maoa = analysis.anova.result.get('MAOA')
  tac1 = analysis.anova.result.get('TAC1')

  assert(0.95 <= maoa <= 1.0)
  assert(0.35 <= tac1 <= 0.52)
  
