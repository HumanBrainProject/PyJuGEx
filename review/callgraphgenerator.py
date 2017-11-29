#!/usr/bin/env python
#sudo pip2 install pycallgraph
#USED PYTHON2
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
import pyjugex
import hbp_human_atlas as atlas

def main():
    graphviz = GraphvizOutput()
    graphviz.output_file = 'callgraphwithoutcache.png'
    genelist = ['ADRA2A', 'AVPR1B', 'CHRM2']
    roi1 = atlas.jubrain.probability_map('FP1', atlas.MNI152)
    roi2 = atlas.jubrain.probability_map('FP2', atlas.MNI152)	
    jugex = pyjugex.PyJugex(cache='.pyjugex', verbose=True)
    with PyCallGraph(output=graphviz):
        result = jugex.DifferentialAnalysis(genelist, roi1, roi2)
	if len([id for id in result if result[id] < .05]) > 0:
            print('Differentially expressed genes are : ')
	    print([id for id in result if result[id] < .05])
	else:
	    print('There are no differentially expressed genes in the given regions')

if __name__ == '__main__':
    main()
