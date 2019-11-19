from statsmodels.formula.api import ols
import numpy as np
import multiprocessing

class Anova:
  """
  Class for carrying out anova analysis with initialised parameters.

  Required parameters include:
  - factors
  - n_genes
  - n_rep
  """
  def __init__(self):
    self.factors = {
      area: None,
      specimen: None,
      age: None,
      race: None,
      z_scores: None
    }

    self.n_genes = []
    self.result = []
    self.F_vec_ref_anovan = None
    self.n_rep = 1000

    # combined_zscores for single_probe_mode
    # genesymbol_and_mean_zscores['combined_zscores'] for multi probe mode

  def first_iteration(self):
    """
    Perform one iteration of ANOVA. Use output of this to populate F_vec_ref_anovan which becomes initial estimate of n_rep passes of FWE.
    """
    self.F_vec_ref_anovan = np.zeros(self.n_genes)
    for i in range(self.n_genes):
      if self.single_probe_mode:
        self.factors['z_scores'] = self.combined_zscores[:,i]
      else:
        self.factors['z_scores'] = self.genesymbol_and_mean_zscores['combined_zscores'][:,i]
          
      mod = ols('z_scores ~ Area + Specimen + Age + Race', data=self.factors).fit()
      aov_table = sm.stats.anova_lm(mod, typ=1)
      if self.verbose:
        logging.getLogger(__name__).info('aov table: {}'.format(aov_table))

      #F_vec_ref_anovan is used as an initial condition to F_mat_perm_anovan in fwe_correction
      self.F_vec_ref_anovan[i] = aov_table['F'][0]

  def fwe_correction(self):
    """
    Perform n_rep passes of FWE using gene_id_and_pvalues of first_iteration() as an initial guess
    """
    invn_rep = 1/self.n_rep
    initial_guess_F_vec = self.F_vec_ref_anovan
    pool = multiprocessing.Pool()
    pass