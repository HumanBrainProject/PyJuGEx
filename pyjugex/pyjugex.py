from .util import (
  get_mean_zscores,
  filter_coordinates_and_zscores,
  from_brainmap_retrieve_gene,
  from_brainmap_retrieve_specimen,
  from_brainmap_on_genes_retrieve_data,
  from_brainmap_retrieve_specimen_factors,
  from_brainmap_retrieve_microarray_filterby_donorids_probeids
)

from .error import (
  NotYetImplementedError,
  ValueMissingError
)
from .anova import anova
import json
import numpy as np

class analysis:
  """
  Usage:

  ```python
  import nibabel as nib
  nii1 = nib.load('nii1.nii.gz')
  nii2 = nib.load('nii2.nii.gz')

  # either set variable during init
  new_analysis = analysis(n_rep=2000, gene_list=['MAOA', 'TAC1'], roi1=nii1)

  # or set them after init
  new_analysis.roi2 = nii2
  new_analysis.single_proble.mode = True

  # carry out analysis, this may take awhile
  result = new_analysis.run()

  # the result is also stored in the analysis instance, but will be overwritten if a new run is called
  print(new_analysis.result)
  ```

  """

  def __init__(self, n_rep=1000, filter_threshold=0.2, single_probe_mode=False, roi1=None, roi2=None, gene_list=[], verbose=False):
    """
    TODO doc
    """
    self.n_rep=n_rep
    self.filter_threshold=filter_threshold
    self.single_probe_mode=single_probe_mode
    self.verbose=verbose
    self.gene_list=gene_list
    self.roi1=roi1
    self.roi2=roi2

    self.anova = None

  def _check_prereq(self):
    """
    Checks that the members of all the prereq has been met to carry out run
    """
    error_message = []
    if self.roi1 is None:
      error_message.append(' roi1 is missing ')
    if self.roi2 is None:
      error_message.append(' roi2 is missing ')
    if len(self.gene_list) == 0:
      error_message.append(' gene_list is empty ')
    if len(error_message) > 0:
      raise ValueMissingError(','.join(error_message))
    return True
  
  def get_filtered_coord(self):
    """
    TODO write doc
    """

    gene_symbols = get_gene_symbols(self.gene_list)
    samples_zscores_and_specimen_dict = from_brainmap_on_genes_retrieve_data(genes=self.gene_list)
    
    filtered_coords_and_zscores = self.get_filtered_coords_and_zscores(samples_zscores_and_specimen_dict)
    combined_zscores = [roi_coord_and_zscore['zscores'][i] for roi_coord_and_zscore in filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
    genesymbol_and_mean_zscores = get_mean_zscores(gene_symbols, combined_zscores)

    areainfo = {}
    
    for roi_coord_zscore in filtered_coords_and_zscores:
      key = roi_coord_zscore['realname']
      if key not in areainfo:
        areainfo[key] = []
      i = 0
      for c in roi_coord_zscore['coords']:
        areainfo[key].append({'xyz' : np.transpose(np.matmul(self.roi1.affine,np.transpose(np.append(c,1))))[0:3].tolist(), 'winsorzed_mean' : genesymbol_and_mean_zscores['combined_zscores'][i].tolist()})
        #areainfo[key].append({'xyz' : c.tolist(), 'winsorzed_mean' : self.genesymbol_and_mean_zscores['combined_zscores'][i].tolist()})
        i = i+1

    return areainfo

  def get_filtered_coords_and_zscores(self, samples_zscores_and_specimen_dict):

    filtered_coords_and_zscores = []

    for o_index, specimen_info in enumerate(samples_zscores_and_specimen_dict['specimen_info']):
      for index, roi in enumerate([self.roi1, self.roi2]):
        filtered_coords_and_zscores.append(
          filter_coordinates_and_zscores(
            roi_nii=roi,
            roi_name=f'roi{index + 1}',
            index=index,
            index_to_samples_zscores_and_specimen_dict=samples_zscores_and_specimen_dict['samples_and_zscores'][o_index],
            specimen=specimen_info,
            filter_threshold=self.filter_threshold
          )
        )
    return filtered_coords_and_zscores

  def run(self):
    """
    Start differential analysis 
    """
    self._check_prereq()

    samples_zscores_and_specimen_dict = from_brainmap_on_genes_retrieve_data(genes=self.gene_list)

    filtered_coords_and_zscores = self.get_filtered_coords_and_zscores(samples_zscores_and_specimen_dict)

    specimen_factors = from_brainmap_retrieve_specimen_factors()
    specimen=[roi_coord_and_zscore['specimen'] for roi_coord_and_zscore in filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]
    
    # Both Age and Race should have len(self.filtered_coords_and_zscores) entries. The following three lines are used to get the correct values from specimenFactors['Age'] and specimenFactors['Race'] using specimenFactors['name'] and repeat them the required number of times as given by self.anova_factors['Specimen']
    combined_zscores = [roi_coord_and_zscore['zscores'][i] for roi_coord_and_zscore in filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))]

    gene_symbols = get_gene_symbols(self.gene_list)

    self.anova = anova(
      area=[roi_coord_and_zscore['name'] for roi_coord_and_zscore in filtered_coords_and_zscores for i in range(len(roi_coord_and_zscore['zscores']))],
      specimen=specimen,
      age=[specimen_factors['age'][specimen_factors['name'].index(specimen_name)] for ind, specimen_name in enumerate(specimen)],
      race=[specimen_factors['race'][specimen_factors['name'].index(specimen_name)] for ind, specimen_name in enumerate(specimen)],
      gene_symbols=gene_symbols,
      combined_zscores=combined_zscores,
      n_rep=self.n_rep
      )

    # in webjugex, accumulate_roicoords_and_name gets called to return the relevant coord
    self.anova.run()
  

def get_gene_symbols(genes=[]):
  gene_symbols = []
  for gene in genes:
    data=from_brainmap_retrieve_gene(gene=gene)
    gene_symbols = gene_symbols + [gene for donor in data['Response']['probes']['probe']]
  return gene_symbols
