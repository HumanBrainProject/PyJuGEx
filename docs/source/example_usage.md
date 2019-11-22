# Example Usage

```python

from pyjugex import PyjugexAnalysis
import nibabel as nib

gene_list=['MAOA','TAC1']

# loading masked nifti volumes

nii1 = nib.load('data/hoc1_th_l.nii.gz')
nii2 = nib.load('data/hoc2_th_l.nii.gz')
analysis = PyjugexAnalysis(
  n_rep=1000,
  gene_list=gene_list,
  roi1 = nii1,
  roi2 = nii2
)

filtered_coord = analysis.get_filtered_coord()
assert(len(filtered_coord['roi1']) == 12)
assert(len(filtered_coord['roi2']) == 11)

analysis.differential_analysis() # Go grab a coffee

maoa = analysis.anova.result.get('MAOA')
tac1 = analysis.anova.result.get('TAC1')

assert(0.95 <= maoa <= 1.0)
assert(0.40 <= tac1 <= 0.52)
  
# same parameter, except n_rep

analysis.n_rep = 10000
analysis.differential_analysis() # Really go grab a coffee

maoa = analysis.anova.result.get('MAOA')
tac1 = analysis.anova.result.get('TAC1')

assert(0.95 <= maoa <= 1.0)
assert(0.40 <= tac1 <= 0.52)
```