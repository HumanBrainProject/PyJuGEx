import analysis as mymodule
import nibabel as nib
import urllib.request
import json
import csv

def DifferentialGeneExpression(gene_cache):
    analysisjugex = mymodule.Analysis(gene_cache)
    return analysisjugex

def readCSVFile(filename):
    rows = dict();
    rows['probe_id'] = []
    rows['gene_symbol'] = []
    rows['entrez_id'] = []
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            rows['probe_id'].append(row['probe_id'])
            rows['gene_symbol'].append(row['gene_symbol'])
            rows['entrez_id'].append(int(row['entrez_id']))
    return rows

