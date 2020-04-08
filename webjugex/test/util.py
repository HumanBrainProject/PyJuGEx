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

from webjugex import util
import os
import json

test_nii_url = 'https://neuroglancer.humanbrainproject.eu/precomputed/JuBrain/17/icbm152casym/pmaps/OFC_Fo1_l_N10_nlin2MNI152ASYM2009C_3.4_publicP_b76752e4ec43a64644f4a66658fed730.nii.gz'

def test_get_pmap():
  resp = util.get_pmap(test_nii_url)
  assert resp.ok

def test_read_byte_via_nib():
  resp = util.get_pmap(test_nii_url)
  img_array = util.read_byte_via_nib(resp.content)
  assert img_array.get_shape() == (193, 229, 193)

def test_is_gzipped():
  not_gzipped = 'test.nii'
  not_gzipped_1 = 'https://blabla.com/post?test.nii'
  is_gzipped = 'test.nii.gz'
  is_gzipped_1 = 'https://blabla.com/post?test.nii.gz'

  assert util.is_gzipped(not_gzipped) == False
  assert util.is_gzipped(not_gzipped_1) == False
  assert util.is_gzipped(is_gzipped)
  assert util.is_gzipped(is_gzipped_1)

def test_from_brainmap_retrieve_gene():
  resp_dist = util.from_brainmap_retrieve_gene('MAOA')
  total_rows = int(resp_dist['@total_rows'])
  assert total_rows > 0

def test_from_brainmap_retrieve_specimen():
  resp_dict = util.from_brainmap_retrieve_specimen('H0351.1015')
  assert resp_dict['success']

def test_from_brainmap_retrieve_microarray_filterby_donorids_probeids():
  gene = util.from_brainmap_retrieve_gene('MAOA')
  probes = gene['Response']['probes']['probe']
  probe_ids = list(map(lambda item: item['id'], probes))
  donor_id = '15496'
  resp = util.from_brainmap_retrieve_microarray_filterby_donorids_probeids(donor_id=donor_id, probe_ids=probe_ids)
  assert resp['success']
