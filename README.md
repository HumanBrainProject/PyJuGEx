# WebJugex
Find a set of differentially expressed genes between two user defined volumes of interest based on JuBrain maps. The tool downloads expression values of user specified sets of genes from Allen Brain API. Then, it uses zscores to find which genes are expressed differentially between the user specified regions of interests. This tool has been integrated into the 3D atlas viewer, where the user can select two JuBrain regions, enlist desired set of genes either by using an autocomplete feature or uploading a csv file. After the analysis is finished, the tool displays the genes and their calculated p values. The user also has the option of downloading the gene names and their p values and the ROI coordinates used in the analysis. For the python package, please take a look at the master branch.

## Dependencies
* numpy
* scipy
* statsmodels
* requests
* nibabel
* xmltodict

### Installing
```
git clone https://github.com/haimasree/Jugex.git
cd Jugex
python setup.py install
```
### Usage
A typical usage is as follows -
```
import hbp_human_atlas as atlas
import analysispyjugex
#import hbp_human_atlas_from_metadata

genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH', 'CRHR1', 'CRHR2', 'GAD2', 'HTR1A', 'HTR1B', 'HTR1D', 'HTR2A', 'HTR3A', 'HTR5A', 'MAOA', 'PDE1A', 'SLC6A2', 'SLC6A4', 'SST', 'TAC1', 'TPH1', 'GPR50', 'CUX2', 'TPH2']

roi1 = {}
roi2 = {}
roi1['name'] = 'FP1'
roi1['data'] = atlas.jubrain.probability_map('FP1', atlas.MNI152)
roi2['name'] = 'FP2'
roi2['data'] = atlas.jubrain.probability_map('FP2', atlas.MNI152)

jugex = analysispyjugex.Analysis(gene_cache_dir='.pyjugex', verbose=True)
result = jugex.DifferentialAnalysis(genelist, roi1, roi2)
print(result)
```

## Versioning
0.1

## Authors

* Big Data Analytics Group, INM-1, Research Center Juelich
## Acknowledgments

* Dr. Sebastian Bludau
* Dr. Timo Dickscheid and other members of BDA-INM1 

