from pyjugex import PyjugexAnalysis
import nibabel as nib

gene_list=['MAOA','TAC1']

nii1 = nib.load('hoc1_th_l.nii.gz')
nii2 = nib.load('hoc2_th_l.nii.gz')
analysis = PyjugexAnalysis(
  gene_list=gene_list,
  roi1 = nii1,
  roi2 = nii2
)
analysis.differential_analysis()

print(analysis.get_filtered_coord())
print(analysis.anova.result)