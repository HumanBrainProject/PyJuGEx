import sys
sys.path.append('..')

import os
nii_dir = os.path.dirname(os.path.abspath(__file__))

from pyjugex import PyjugexAnalysis
import nibabel as nib

def test_e2e():

  gene_list=['MAOA','TAC1']

  nii1 = nib.load(os.path.join(nii_dir, 'data/hoc1_th_l.nii.gz'))
  nii2 = nib.load(os.path.join(nii_dir, 'data/hoc2_th_l.nii.gz'))
  analysis = PyjugexAnalysis(
    gene_list=gene_list,
    roi1 = nii1,
    roi2 = nii2
  )

  filtered_coord = analysis.get_filtered_coord()
  assert(len(filtered_coord['roi1']) == 12)
  assert(len(filtered_coord['roi2']) == 11)

  analysis.differential_analysis()

  maoa = analysis.anova.result.get('MAOA')
  tac1 = analysis.anova.result.get('TAC1')

  assert(0.95 <= maoa <= 1.0)
  assert(0.42 <= tac1 <= 0.52)
  