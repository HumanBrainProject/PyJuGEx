# PyJugex
Find a set of differentially expressed genes between two user defined volumes of interest based on JuBrain maps. The tool downloads expression values of user specified sets of genes from Allen Brain API. Then, it uses zscores to find which genes are expressed differentially between the user specified regions of interests. This tool is available as a Python package. 

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

genelist = ['ADRA2A', 'AVPR1B', 'CHRM2', 'CNR1', 'CREB1', 'CRH']

roi1 = atlas.jubrain.probability_map('FP1', atlas.MNI152)
roi2 = atlas.jubrain.probability_map('FP2', atlas.MNI152)

jugex = analysispyjugex.Analysis(gene_cache='.pyjugex', verbose=True)
result = jugex.DifferentialAnalysis(genelist, roi1, roi2)
if len([id for id in result if result[id] < .05]) > 0:
    print('Differentially expressed genes are : ')
    print([id for id in result if result[id] < .05])
else:
    print('There are no differentially expressed genes in the given regions')
```

## Versioning
0.1

## Authors

* **Haimasree Bhattacharya** - *h.bhattacharya@fz-juelich.de* 

## Acknowledgments

* Dr. Sebastian Bludau
* Dr. Timo Dickscheid and other members of BDA-INM1 

