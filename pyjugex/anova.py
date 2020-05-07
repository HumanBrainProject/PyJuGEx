from statsmodels.formula.api import ols
import statsmodels.api as sm
import numpy as np
import multiprocessing
from .util import get_mean_zscores
from .error import ValueMissingError

def _unwrap_self_do_anova_with_permutation_rep(*args, **kwargs):
  """
  Helper function to enable usage of the multiprocessing module inside a class.
  Args:
  *args: Variable length argument list.
  **kwargs: Arbitrary keyword arguments.
  Returns:
  do_anova_with_permutation_rep()
  """
  return anova.do_anova_with_permutation_rep(*args, **kwargs)


class anova:
  """
  Class for carrying out anova analysis with initialised parameters.

  Required parameters include:
  - factors
  - n_genes
  - n_rep
  """
  def __init__(self, z_scores=None, area=None, specimen=None, age=None, race=None, combined_zscores=None, gene_symbols=None, n_rep=1000, verbose=True):
    """
    init function for anova
    """
    if combined_zscores is None:
      raise ValueMissingError('combined_zscores is required')
  
    if gene_symbols is None:
      raise ValueMissingError('gene_symbols is required')

    self.factors = {
      "Area": area,
      "Specimen": specimen,
      "Age": age,
      "Race": race,
      "Zscores": z_scores
    }

    self.genesymbol_and_mean_zscores = get_mean_zscores(gene_symbols, combined_zscores)
    self.n_genes = len(self.genesymbol_and_mean_zscores['combined_zscores'][0])

    self.F_vec_ref_anovan = None
    self.n_rep = n_rep
    self.F_mat_perm_anovan = None

    self.verbose = verbose
    self.result = None

  def run(self):
    self._first_iteration()
    self._fwe_correction()
    self._collate_result()

  def _collate_result(self):
    def div_func(arr):
      return np.count_nonzero(arr)/self.n_rep

    FWE_corrected_p = np.apply_along_axis(div_func, 0, self.F_mat_perm_anovan.max(1)[:, np.newaxis] >= np.array(self.F_vec_ref_anovan))
    self.result = dict(zip(self.genesymbol_and_mean_zscores['uniqueId'], FWE_corrected_p))

  def _first_iteration(self):
    """
    Perform one iteration of ANOVA. Use output of this to populate F_vec_ref_anovan which becomes initial estimate of n_rep passes of FWE.
    """
    self.F_vec_ref_anovan = np.zeros(self.n_genes)
    for i in range(self.n_genes):
      self.factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,i]
          
      mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.factors).fit()
      aov_table = sm.stats.anova_lm(mod, typ=1)

      #F_vec_ref_anovan is used as an initial condition to F_mat_perm_anovan in _fwe_correction
      self.F_vec_ref_anovan[i] = aov_table['F'][0]

  def _fwe_correction(self):
    """
    Perform n_rep passes of FWE using gene_id_and_pvalues of _first_iteration() as an initial guess
    """
    initial_guess_F_vec = self.F_vec_ref_anovan
    pool = multiprocessing.Pool()

    # double check math here
    self.F_mat_perm_anovan = np.array(pool.map(_unwrap_self_do_anova_with_permutation_rep, [self]*(self.n_rep-1)))
    self.F_mat_perm_anovan = np.insert(self.F_mat_perm_anovan, 0, initial_guess_F_vec, axis=0)

  def do_anova_with_permutation_gene(self, index_to_gene_list):
    """
    Perform one repetition of anova for each gene
    Args:
    index_to_gene_list (int) : Index into the genesymbol_and_mean_zscores['combined_zscores'] array, representing mean zscore of a gene.
    Returns:
    float: F value extracted from the anova table
    """
    self.factors['Area'] = np.random.permutation(self.factors['Area'])
    self.factors['Zscores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,index_to_gene_list]
    
    mod = ols('Zscores ~ Area + Specimen + Age + Race', data=self.factors).fit()
    aov_table = sm.stats.anova_lm(mod, typ=1)
    return aov_table['F'][0]
  
  def do_anova_with_permutation_rep(self):
    """
    Perform one repetition of anova for all genes
    Returns:
    list: a list of F_values, one for each gene.
    """
    return list(map(self.do_anova_with_permutation_gene, range(0,self.n_genes)))