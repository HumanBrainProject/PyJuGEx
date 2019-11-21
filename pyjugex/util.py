import requests
import tempfile
import nibabel as nib
import os
import re
import logging
import xmltodict
import numpy as np
import scipy as sp

from .cache import MemoryCache

cache = MemoryCache()

def get_filename_from_resp(resp):
  # determine the type of the file. look at the disposition header, use PMapURL as a fallback
  content_disposition_header = resp.headers.get('content-disposition')
  filename = re.search(r'filename=(.*?)$', content_disposition_header).group(1) if content_disposition_header is not None and re.search(r'filename=(.*?)$', content_disposition_header) is not None else resp.url
  return filename

# given url as either a string or obj, interpretes, and performs get/post request
# returns resp
# may raise HTTP exception

def get_pmap(url, json=None):
  if json is None:
    resp = requests.get(url)
  else:
    resp = requests.post(url, json=json)
  
  resp.raise_for_status()
  return resp

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
  return re.search(r"\.gz$", filename) is not None

def from_brainmap_retrieve_gene(gene, verbose=False):

  """
  Retrieve probe ids for the given gene lists, update self.probe_ids which will be used by download_and_save_zscores_samples() or download_and_save_zscores_samples_partial() to
  form the url and update self.gene_symbols to be used by get_mean_zscores()
  """

  _key = f'from_brainmap_retrieve_gene__{gene}'
  if cache.get_from_key(_key) is not None:
    return cache.get_from_key(_key)

  base_retrieve_probe_ids = "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Probe,rma::criteria,[probe_type$eq'DNA'],products[abbreviation$eq'HumanMA'],gene[acronym$eq"
  end_retrieve_probe_ids = "],rma::options[only$eq'probes.id']"

  url = '{}{}{}'.format(base_retrieve_probe_ids, gene, end_retrieve_probe_ids)

  if verbose:
    logging.getLogger(__name__).info('url: {}'.format(url))

  resp = requests.get(url)
  resp.raise_for_status()

  resp_dict = xmltodict.parse(resp.text)
  cache.store_key_value(_key, resp_dict)
  return resp_dict

def from_brainmap_retrieve_specimen(specimen_id, verbose=False):

  """
  Download names and transformation matrix for each specimen/donor from Allen Brain Api and save them on disk as specimenName.txt
  and specimenMat.txt respectively, load.
  """

  _key=f'from_brainmap_retrieve_specimen__{specimen_id}'
  if cache.get_from_key(_key) is not None:
    return cache.get_from_key(_key)

  base_url_download_specimens = "http://api.brain-map.org/api/v2/data/Specimen/query.json?criteria=[name$eq"+"'"
  end_url_download_specimens = "']&include=alignment3d"

  url = '{}{}{}'.format(base_url_download_specimens, specimen_id, end_url_download_specimens)
  resp = requests.get(url)
  resp.raise_for_status()

  resp_json = resp.json()
  cache.store_key_value(_key, resp_json)
  return resp_json

def from_brainmap_retrieve_microarray_filterby_donorids_probeids(donor_id, probe_ids, verbose=False):

  """
  Query Allen Brain Api for given set of genes for the donor given by donor_id
  Args:
  donor_id (int): Id of a donor which is used to query Allen Brain API.
  Returns:
  dict: A dictionary representing the just downloaded samples, probes and zscores for the given donor_id and the probes given by self.probe_ids.
  """
  _key = f'from_brainmap_retrieve_microarray_filterby_donorids_probeids__{donor_id}__{"_".join(probe_ids)}'
  if cache.get_from_key(_key) is not None:
    return cache.get_from_key(_key)

  base_query_api = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
  probes = ','.join(probe_ids)
  end_query_api = "][donors$eq{}]".format(donor_id)
  url = '{}{}{}'.format(base_query_api, probes, end_query_api)

  resp = requests.get(url)

  resp.raise_for_status()
  return resp.json()

# TODO cleanup
# TODO need tests
def transform_samples_MRI_to_MNI152(samples, transformation_mat):
    """
    Convert the MRI coordinates of samples to MNI152 space
    Args:
          samples (dict): Contains mri coordinates, well and polygon id for each sample used in Allen Brain.
          transformation_mat (numpy.ndarray): A 4x4 numpy array to convert the above mentioned MRI coordinates to MNI152 space.
    Returns:
          dict: A dictionary containing three keys -
            mnicoords - two dimensional numpy array where each row represents a three dimensional coordinate in MNI152 space, for all the samples.
            well - list of well id for the sample the respective coordinate belongs to
            polygon -  list of polygon id for the sample the respective coordinate belongs to
    """
    np_T = np.array(transformation_mat[0:3, 0:4])
    mri = np.vstack(s['sample']['mri'] for s in samples)
    add = np.ones((len(mri), 1), dtype=np.int)
    mri = np.append(mri, add, axis=1)
    mri = np.transpose(mri)
    coords = np.matmul(np_T, mri)
    coords = coords.transpose()
    well = [s['sample']['well'] for s in samples]
    polygon = [s['sample']['polygon'] for s in samples]
    return {'mnicoords' : coords, 'well' : well, 'polygon' : polygon}
    #return coords

# TODO write test
def filter_coordinates_and_zscores(roi_nii, index_to_samples_zscores_and_specimen_dict, specimen, index, roi_name='Unnamed ROI', filter_threshold=0.2):
  """
  Populate self.filtered_coords_and_zscores with zscores and coords for samples which belong to a particular specimen and spatially represented in the given roi.
  Args:
    roi (nib.nifti1.Nifti1Image): probability map of a region of interest.
    index_to_samples_zscores_and_specimen_dict (dict): Index into samples_zscores_and_specimen_dict
    specimen (dict): dictionary representing a specimen with its name and transformation matrix as keys
    index (int): 0 or 1, representing which region of interest it is.
  Returns:
    dict : Contains the following keys -
      a) zscores - Lists of zscore corresponding to the Allen Brain coordinates (in MNI152 space) which are spatially represented in region of interest given by roi parameter.
      b) coords - Lists of Allen Brain coordinates (in MNI152 space) which are spatially represented in region of interest given by roi parameter.
      c) coord_well - Lists of well id for all the samples which are spatially represented in region of interest given by roi parameter.
      d) coord_polygon - Lists of polygon id for all the samples which are spatially represented in region of interest given by roi parameter.
      e) specimen - same as specimen['name'].
      f) name - 'img1' representing first region of interest, 'img2' representing second region of interest.
  """
  revised_samples_zscores_and_specimen_dict = dict.fromkeys(['zscores', 'coords', 'coord_well', 'coord_polygon', 'specimen', 'name'])
  revised_samples_zscores_and_specimen_dict['realname'] = roi_name
  revised_samples_zscores_and_specimen_dict['name'] = 'img{}'.format(str(index+1))
  img_arr = roi_nii.get_data()
  invroiMni = np.linalg.inv(roi_nii.affine)
  T = np.dot(invroiMni, specimen['alignment3d'])
  '''
  coords = transform_samples_MRI_to_MNI152(index_to_samples_zscores_and_specimen_dict['samples'], T)
  coords = (np.rint(coords)).astype(int)
  coords = [np.array([-1, -1, -1]) if (coord > 0).sum() != 3 or img_arr[coord[0],coord[1],coord[2]] <= self.filter_threshold or img_arr[coord[0],coord[1],coord[2]] == 0 else coord for coord in coords]
  revised_samples_zscores_and_specimen_dict['coords'] = [coord for coord in coords if (coord > 0).sum() == 3]
  '''
  coords_dict = transform_samples_MRI_to_MNI152(index_to_samples_zscores_and_specimen_dict['samples'], T)
  coords = (np.rint(coords_dict['mnicoords'])).astype(int)
  coords = [np.array([-1, -1, -1]) if (coord > 0).sum() != 3 or img_arr[coord[0],coord[1],coord[2]] <= filter_threshold or img_arr[coord[0],coord[1],coord[2]] == 0 else coord for coord in coords]
  revised_samples_zscores_and_specimen_dict['coords'] = [coord for coord in coords if (coord > 0).sum() == 3]
  revised_samples_zscores_and_specimen_dict['coord_well'] = [well for (coord, well) in zip(coords, coords_dict['well']) if (coord > 0).sum() == 3]
  revised_samples_zscores_and_specimen_dict['coord_polygon'] = [polygon for (coord, polygon) in zip(coords, coords_dict['polygon']) if (coord > 0).sum() == 3]
  revised_samples_zscores_and_specimen_dict['zscores'] = [zscore for (coord, zscore) in zip(coords, index_to_samples_zscores_and_specimen_dict['zscores']) if (coord > 0).sum() == 3]
  revised_samples_zscores_and_specimen_dict['specimen'] = specimen['name']
  return revised_samples_zscores_and_specimen_dict

# TODO write test
def from_brainmap_retrieve_specimen_factors():
  """
  Download various factors such as age, name, race, gender of the six specimens from Allen Brain Api, save them at cache/specimenFactors.txt and create a dict.

  """
  _key='from_brainmap_retrieve_specimen_factors'
  if cache.get_from_key(_key) is not None:
    return cache.get_from_key(_key)

  url_build_specimen_factors = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]"

  resp = requests.get(url_build_specimen_factors)

  resp.raise_for_status()
  resp_json = resp.json()

  res = resp_json['msg']
  specimen_factors={}
  specimen_factors['id'] = [r['id'] for r in res]
  specimen_factors['name'] = [r['name'] for r in res]
  specimen_factors['race'] = [r['race_only'] for r in res]
  specimen_factors['gender'] = [r['sex'] for r in res]
  specimen_factors['age'] = [r['age']['days']/365 for r in res]

  cache.store_key_value(_key, specimen_factors)
  return specimen_factors

def from_brainmap_on_genes_retrieve_data(genes=[]):
  """
  Based on selected genes, return data from Allen Institute
  """

  _key=f'from_brainmap_on_genes_retrieve_data__{"_".join(genes)}'
  if cache.get_from_key(_key) is not None:
    return cache.get_from_key(_key)

  donor_ids = ['15496', '14380', '15697', '9861', '12876', '10021'] #HARDCODING donor_ids
  specimens  = ['H0351.1015', 'H0351.1012', 'H0351.1016', 'H0351.2001', 'H0351.1009', 'H0351.2002']
  
  samples_zscores_and_specimen_dict = {
    "specimen_info":[],
    "samples_and_zscores":[]
  }

  if not len(genes) > 0:
    raise ValueMissingError('genes are required')

  probe_ids = []
  for gene in genes:
    data=from_brainmap_retrieve_gene(gene=gene)

    if int(data['Response']['@num_rows']) <= 0:
      raise ValueError('Please check the spelling of {}. No such gene exists in Allen Brain API.'.format(gene))
    probe_ids = probe_ids + [ donor['id'] for donor in data['Response']['probes']['probe'] ]

  for specimen in specimens:
    text = from_brainmap_retrieve_specimen(specimen)
    samples_zscores_and_specimen_dict['specimen_info'] = samples_zscores_and_specimen_dict['specimen_info'] + [get_specimen_data(text['msg'][0])]

  for donor_id in donor_ids:
    text = from_brainmap_retrieve_microarray_filterby_donorids_probeids(probe_ids=probe_ids, donor_id=donor_id)
    data = text['msg']
    
    # @TODO clean this
    zscores = np.array([[float(data['probes'][i]['z-score'][j]) for i in range(len(data['probes']))] for j in range(len(data['samples']))])

    samples_zscores_and_specimen_dict['samples_and_zscores'].append({
      'samples' : data['samples'],
      'probes' : data['probes'],
      'zscores' : zscores
      })

  cache.store_key_value(_key, samples_zscores_and_specimen_dict)
  return samples_zscores_and_specimen_dict

# TODO write test
def get_mean_zscores(gene_symbols=None, combined_zscores=None):
  """
  Compute Winsorzed mean of zscores over all probes associated with a given gene. combined_zscores have zscores for all the probes and all the valid coordinates.
  As a gene_id_and_pvalues you get a numpy array of size len(self.filtered_coords_and_zscores)xlen(self.gene_list). self.genesymbol_and_mean_zscores['combined_zscores'][i][j] returns the     winsorzed mean of jth gene taken over all the probes corresponding to the ith valid sample.
  Args:
  combined_zscores (list): lists of zscores corresponding to each region of interest, populated from filtered_coords_and_zscores
  """
  
  if gene_symbols is None:
    raise ValueMissingError('gene_symbols needs to be defined')
  if combined_zscores is None:
    raise ValueMissingError('combined_zscores is required')
  
  unique_gene_symbols = np.unique(gene_symbols)
  
  '''
  A = [a,a,a,b,b,b,c,c]
  B = [a,b,c]
  Following line of code will give  indices = [[0,1,2],[3,4,5],[6,7]]
  '''
  indices = [np.where(np.in1d(gene_symbols, genesymbol))[0] for genesymbol in unique_gene_symbols]

  # this is more legible
  '''
  for i in range (len(unique_gene_symbols)):
    for j in range(len(combined_zscores)):
      for k in range(len(indices[i])):
        tmp[j] = combined_zscores[j][indices[i][k]][:]
    winsorzed_mean_zscores[j][i] = np.mean(sp.stats.mstats.winsorize(tmp[j], limits=0.1))
  '''
  winsorzed_mean_zscores =  np.array([[np.mean(sp.stats.mstats.winsorize([combined_zscores[j][indices[i][k]] for k in range(0, len(indices[i]))], limits=0.1)) for i in range (len(unique_gene_symbols))] for j in range(len(combined_zscores))])
  return {
    "uniqueId": unique_gene_symbols,
    "combined_zscores": winsorzed_mean_zscores
  }

def get_specimen_data(specimen_metadata):
  """
  For each specimen, extract the name and alignment matrix and stores a dict
  Args:
  specimen_metadata (dict): Contains metadata for specimens used in Allen Brain. They can be downloaded through Allen Brain API.
  Returns:
  dict: specimen dict with two keys -
  a) name: name of the specimen
  b) alignment3d: transformation matrix to convert from MRI to MNI152 space.
  """
  specimen = dict()
  specimen['name'] = specimen_metadata['name']
  x = specimen_metadata['alignment3d']
  specimen['alignment3d'] = np.array([
  [x['tvr_00'], x['tvr_01'], x['tvr_02'], x['tvr_09']],
  [x['tvr_03'], x['tvr_04'], x['tvr_05'], x['tvr_10']],
  [x['tvr_06'], x['tvr_07'], x['tvr_08'], x['tvr_11']],
  [0, 0, 0, 1]])
  return specimen