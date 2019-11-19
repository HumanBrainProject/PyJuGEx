from util import ValueMissingError
import json
with open('data/genesymbols.json', 'r') as f:
  gene_symbol_arr = f.read()

class PyjugexAnalysis:
  """
  Rewritten PyjugexAnalysis class

  Usage:

  import nibabel as nib
  nii1 = nib.load('nii1.nii.gz')
  nii2 = nib.load('nii2.nii.gz')

  # either set variable during init
  new_analysis = PyjugexAnalysis(n_rep=2000, gene_list=['MAOA', 'TAC1'], roi1=nii1)

  # or set them after init
  new_analysis.roi2 = nii2
  new_analysis.single_proble.mode = True

  # carry out analysis, this may take awhile
  result = new_analysis.differential_analysis()

  # the result is also stored in the analysis instance, but will be overwritten if a new differential_analysis is called
  print(new_analysis.result)
  """

  def __init__(self, n_rep=1000, filter_threshold=0.2, single_probe_mode=False, roi1=None, roi2=None, gene_list=[], verbose=False):
    self.n_rep=n_rep
    self.filter_threshold=filter_threshold
    self.single_probe_mode=single_probe_mode
    self.verbose=verbose
    self.gene_list=gene_list
    self.roi1=roi1
    self.roi2=roi2

    self.result = None

  @staticmethod
  def get_gene_list():
    """
    Returns list of gene symbols.
    """
    return json.loads(gene_symbol_arr)

  def _check_prereq(self):
    """
    Checks that the members of all the prereq has been met to carry out differential_analysis
    """
    error_message = []
    if self.roi1 is None:
      error_message.append(' roi1 is missing ')
    if self.roi2 is None:
      error_message.append(' roi2 is missing ')
    if len(self.gene_list) == 0:
      error_message.append(' gene_list is empty ')
    if error_message != '':
      raise ValueMissingError(','.join(error_message))
    return true
    
  
  def differential_analysis(self):
    """
    Start differential analysis 
    """
    self._check_prereq()
    return 0