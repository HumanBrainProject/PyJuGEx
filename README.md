# pyjugex

[![Travis (.com)](https://img.shields.io/travis/com/humanbrainproject/pyjugex)](https://travis-ci.com/github/HumanBrainProject/PyJuGEx)
[![PyPI](https://img.shields.io/pypi/v/pyjugex)](https://pypi.org/project/pyjugex/)
[![Read the Docs (version)](https://img.shields.io/readthedocs/pyjugex/latest)](https://pyjugex.readthedocs.io/en/latest/index.html)


Find a set of differentially expressed genes between two user defined volumes of interest based on JuBrain maps. The tool downloads expression values of user specified sets of genes from Allen Brain API. Then, it uses zscores to find which genes are expressed differentially between the user specified regions of interests. This tool is available as a Python package.

JuGEx homepage: <http://www.fz-juelich.de/inm/inm-1/DE/Forschung/_docs/JuGex/JuGex_node.html>

repository: <https://github.com/HumanBrainProject/PyJuGEx>


## Installation

Via pip:

```bash
pip install pyjugex
```

or from source

```bash
git clone https://github.com/HumanBrainProject/PyJuGEx
cd pyjugex
pip install -r requirements.txt
pip install .
```

## Usage

[see example usage](example_usage.md)

## Dependencies
* numpy
* scipy
* statsmodels
* requests
* nibabel
* xmltodict

## Versioning
1.0.1alpha

## Authors
* Big Data Analytics Group, INM-1, Research Center Juelich

## Acknowledgments
* Haimasree Bhattacharya
* Dr. Sebastian Bludau, Dr. Thomas Mühleisen
* Dr. Timo Dickscheid and other members of BDA-INM1

## Reference
Sebastian Bludau, Thomas W. Mühleisen, Simon B. Eickhoff, Michael J. Hawrylycz, Sven Cichon, Katrin Amunts. Integration of transcriptomic and cytoarchitectonic data implicates a role for MAOA and TAC1 in the limbic-cortical network. 2018, Brain Structure and Function. <https://doi.org/10.1007/s00429-018-1620-6>

<!-- ## License

Apache 2.0 -->