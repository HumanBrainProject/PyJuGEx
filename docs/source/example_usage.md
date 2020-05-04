# Example Usage

```python

# (optional) download probabilistic maps of interest using ebrains_atlascore

from ebrains_atlascore import regions
from ebrains_atlascore.util.hemisphere import Hemisphere
from ebrains_atlascore.region import Region
pmap_hoc1 = regions.get_probability_map_for_region(Region('Area-hOc1', parcellation='JuBrain Cytoarchitectonic Atlas', referencespace='MNI152'), Hemisphere.LEFT.value, 0.2)
pmap_hoc2 = regions.get_probability_map_for_region(Region('Area-hOc2', parcellation='JuBrain Cytoarchitectonic Atlas', referencespace='MNI152'), Hemisphere.LEFT.value, 0.2)

with open('hoc1.nii', 'wb') as fp:
  fp.write(pmap_hoc1.data)
with open('hoc2.nii', 'wb') as fp:
  fp.write(pmap_hoc2.data)


# import pyjugex and nibabel
import pyjugex
import nibabel as nib

# setup parameters
gene_list=['MAOA','TAC1']
nii1 = nib.load('hoc1.nii')
nii2 = nib.load('hoc2.nii')

# load parameters and setup analysis
analysis = pyjugex.analysis(
  n_rep=1000,
  gene_list=gene_list,
  roi1 = nii1,
  roi2 = nii2
)

# prior to analysis, one can retrieve the coordinates of the probes in MNI152 space
filtered_coord = analysis.get_filtered_coord()
assert(len(filtered_coord['roi1']) == 12)
assert(len(filtered_coord['roi2']) == 11)

analysis.run() # Go grab a coffee

maoa = analysis.anova.result.get('MAOA')
tac1 = analysis.anova.result.get('TAC1')

assert(0.95 <= maoa <= 1.0)
assert(0.40 <= tac1 <= 0.52)

# results of the differential analysis is saved in the result object
analysis.n_rep = 10000
analysis.run() # Really go grab a coffee

maoa = analysis.anova.result.get('MAOA')
tac1 = analysis.anova.result.get('TAC1')

assert(0.95 <= maoa <= 1.0)
assert(0.40 <= tac1 <= 0.52)
```