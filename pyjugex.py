import analysis as mymodule
import nibabel as nib
import urllib.request
import json
import csv
import requests


def DifferentialGeneExpression(gene_cache):
    analysisjugex = mymodule.Analysis(gene_cache)
    return analysisjugex

def ReadParcellations(dictionaryimages, regionname):
    """
    Download and save a local copy of the masks
    for the respective regions
    """
    url = dictionaryimages[regionname]
    r = requests.get(url, verify=False)
    filename = 'output'+regionname+'.nii.gz'
    with open(filename, 'wb') as f:
        f.write(r.content)
    return filename

