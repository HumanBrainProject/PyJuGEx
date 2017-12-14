# WebJugex
Find a set of differentially expressed genes between two user defined volumes of interest based on JuBrain maps. The tool downloads expression values of user specified sets of genes from Allen Brain API. Then, it uses zscores to find which genes are expressed differentially between the user specified regions of interests. This tool has been integrated into the 3D atlas viewer, where the user can select two JuBrain regions, enlist desired set of genes either by using an autocomplete feature or uploading a csv file. After the analysis is finished, the tool displays the genes and their calculated p values. The user also has the option of downloading the gene names and their p values and the ROI coordinates used in the analysis. For the python package, please take a look at the master branch.

## Dependencies
* node 6+
* numpy
* scipy
* statsmodels
* requests
* nibabel
* xmltodict

### Installing
```
git clone -b webjugex https://github.com/haimasree/Jugex.git
git clone https://github.com/xgui3783/ws_filewatcher.git
cd ws_filewatcher
npm i
```

### Instructions
Following are the steps to get WebJugex up and running on 3D Atlas Viewer -
```
cd Jugex
python aioserver.py
cd ws_filewatcher
npm start
Go to the url and run your plugin
```
## Versioning
0.1

## Authors

* Big Data Analytics Group, INM-1, Research Center Juelich
## Acknowledgments

* Dr. Sebastian Bludau
* Dr. Timo Dickscheid and other members of BDA-INM1 

