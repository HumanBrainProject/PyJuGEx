import os
import time
import json
import logging
import brainscapes
import jwt_handler

from default import default_param

print("start")
token_handler = jwt_handler.jwt_handler()
time.sleep(5)
print("token service started")

def run_pyjugex_analysis(jsonobj):
    print(jsonobj)
    filter_threshold = jsonobj.get('threshold', default_param['threshold'])
    n_rep = jsonobj.get('nPermutations', default_param['nPermutations'])

    os.environ['HBP_AUTH_TOKEN'] = token_handler.token["access_token"]

    brainscapes.logger.setLevel("INFO") # we want to see some messages!

    atlas = brainscapes.atlases.MULTILEVEL_HUMAN_ATLAS

    area1_julich_brain_version = jsonobj["area1"]["areas"][0]["atlas"]["version"]
    area2_julich_brain_version = jsonobj["area2"]["areas"][0]["atlas"]["version"]

    print(area1_julich_brain_version)

    if not area1_julich_brain_version == area2_julich_brain_version:
        print("version mismatch")

    area1_julich_brain_version.replace(".", "_")
    print(area1_julich_brain_version)

    atlas.select_parcellation(brainscapes.parcellations.JULICH_BRAIN_PROBABILISTIC_CYTOARCHITECTONIC_MAPS_V2_5_)
    # as in the original JuGEx, we prefer thresholded probability maps # over the labelled region in the maximum probability map
    #atlas.enable_continuous_map_thresholding(filter_threshold)

    jugex = brainscapes.analysis.DifferentialGeneExpression(atlas)

    for gene in jsonobj['selectedGenes']:
        jugex.add_candidate_genes(gene)

    jugex.define_roi1(jsonobj["area1"]["areas"][0]["name"] + " " + jsonobj["area1"]["areas"][0]["hemisphere"])
    jugex.define_roi2(jsonobj['area2']["areas"][0]["name"] + " " + jsonobj["area2"]["areas"][0]["hemisphere"])

    #from nilearn import plotting import numpy as np
    #for region,samples in zip(['v1 right','v2 right'],[jugex.samples1,jugex.samples2]):
    #atlas.select_region(region)
    #mask = atlas.get_mask(bs.spaces.MNI_152_ICBM_2009C_NONLINEAR_ASYMMETRIC) display = plotting.plot_roi(mask)
    #display.add_markers([k for k,v in samples.items()])

    #for item in jugex._samples1.items():
    #    print(item[1])
    #    for item_1 in item[1]:
    #        print(item_1)
    #        print(item[1][item_1])
    #        print(type(item[1][item_1]))

    jugex.run(permutations=n_rep)
    jugex_result = jugex.result()
    print(jugex_result)

    #logger.log("info", {"jugex_result": str(jugex_result)})

    zscores = jugex_result["zscores"][jsonobj['selectedGenes'][0]]

    for gene in jsonobj['selectedGenes'][1:]:
        zscores = zip(zscores, jugex_result["zscores"][gene])

    #"probes":[{"probe_properties": {"MAOA": 0.0101010, "TAC1": -1232434}, "position": [x,y,z]

    probes_area1 = []
    probes_area2 = []

    brainscapes_areas = list(set(jugex_result["area"]))

    for i in range(len(list(jugex_result["zscores"].values())[0])):
        tmp_probe = {"probe_properties": {}}
        for gene in jsonobj['selectedGenes']:
            tmp_probe["probe_properties"][gene] = jugex_result["zscores"][gene][i]

        tmp_probe["position"] = jugex_result["mnicoord"][i]

        if jugex_result["area"][i] == brainscapes_areas[0]:
            print("Area 1")
            probes_area1.append(tmp_probe)
        else:
            print("Area 2")
            probes_area2.append(tmp_probe)

    print(probes_area1)
    print(probes_area2)

    result = {"result": jugex_result["p-values"]}
    #result["Version"] = os.environ["OPENSHIFT_BUILD_COMMIT"]
    result["Areas"] = []
    result["Areas"].append(
                        {
                            "name": jsonobj["area1"]["areas"][0]["name"],
                            "hemisphere": jsonobj["area1"]["areas"][0]["hemisphere"],
                            "probes": probes_area1
                        })
    result["Areas"].append(
                    {
                        "name": jsonobj["area2"]["areas"][0]["name"],
                        "hemisphere": jsonobj["area2"]["areas"][0]["hemisphere"],
                        "probes": probes_area2
                    })

    print(result)

    #logger.log("info", {"returned_result": str(result)})

    return json.dumps(result)

jsonobj = {"area1":{"threshold":0.2,"areas":[{"name":"Area 3b (PostCG)","hemisphere":"right","atlas":{"name":"Julich Brain","version":"v2.4","id":"minds/core/parcellationatlas/v1.0.0/94c1125b-b87e-45e4-901c-00daee7f2579-25"}}]},"area2":{"threshold":0.2,"areas":[{"name":"Area PFt (IPL)","hemisphere":"left","atlas":{"name":"Julich Brain","version":"v2.4","id":"minds/core/parcellationatlas/v1.0.0/94c1125b-b87e-45e4-901c-00daee7f2579- 25"}}]},"selectedGenes":["MAOA","TAC1"],"nPermutations":2,"threshold":"0.2","singleProbeMode":False,"ignoreCustomProbe":False}

run_pyjugex_analysis(jsonobj)
print("test finished")
