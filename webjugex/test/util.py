from webjugex import util

test_nii_url = 'https://neuroglancer.humanbrainproject.eu/precomputed/JuBrain/17/icbm152casym/pmaps/OFC_Fo1_l_N10_nlin2MNI152ASYM2009C_3.4_publicP_b76752e4ec43a64644f4a66658fed730.nii.gz'

def test_get_pmap():
  resp = util.get_pmap(test_nii_url)
  assert resp.ok

def test_read_byte_via_nib():
  resp = util.get_pmap(test_nii_url)
  img_array = util.read_byte_via_nib(resp.content)
  assert img_array.get_shape() = (193, 229, 193)

def test_is_gzipped():
  not_gzipped = 'test.nii'
  not_gzipped_1 = 'https://blabla.com/post?test.nii'
  is_gzipped = 'test.nii.gz'
  is_gzipped_1 = 'https://blabla.com/post?test.nii.gz'

  assert util.is_gzipped(not_gzipped) == False
  assert util.is_gzipped(not_gzipped_1) == False
  assert util.is_gzipped(is_gzipped)
  assert util.is_gzipped(is_gzipped_1)
