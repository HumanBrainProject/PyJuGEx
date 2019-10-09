import requests
import tempfile
import nibabel as nib
import os
import re

# given url as either a string or obj, interpretes, and performs get/post request
# returns resp
# may raise HTTP exception

def get_pmap(url, body=None):
  resp = requests.get(url) if body is None else requests.post(url, body=body)
  if resp.ok:
    return resp
  else:
    resp.raise_for_status()

# input is byte
# save to cache, then generate a random hashed file name
# specify gzip to append .gz

def read_byte_via_nib(content, gzip=False):
  fp, fp_name = tempfile.mkstemp(suffix='.nii.gz' if gzip else '.nii')
  os.write(fp, content)
  img_array = nib.load(fp_name)
  os.close(fp)
  return img_array

def is_gzipped(filename):
  re.search("\.gz$", filename)
