# Copyright 2020 Forschungszentrum Jülich
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import aiohttp_cors
import json
import os
import requests
import siibra as brainscapes
import siibra_jugex

import HBPLogger

from aiohttp import web
from default import default_param

_fluent_host = os.getenv('FLUENT_HOST', None)
_fluent_protocol = os.getenv('FLUENT_PROTOCOL', None)
_logger_url = '{protocol}://{hostname}/'.format(protocol=_fluent_protocol, hostname=_fluent_host) if _fluent_host is not None and _fluent_protocol is not None else None
_application_name = os.getenv('APPLICATION_NAME', 'webjugex-backend')
_deployment = os.getenv('DEPLOYMENT', 'local')

logger = HBPLogger.HBPLogger(_logger_url,_application_name,_deployment)

with open("files/genesymbols.txt", "r") as f:
    dictAutocompleteString = f.read()


def _transform_brainscapes_response(brainscapes_resp, jsonobj):
    zscores = brainscapes_resp["zscores"][jsonobj['selectedGenes'][0]]

    for gene in jsonobj['selectedGenes'][1:]:
        zscores = zip(zscores, brainscapes_resp["zscores"][gene])

    probes_area1 = []
    probes_area2 = []

    brainscapes_areas = list(set(brainscapes_resp["area"]))

    for i in range(len(list(brainscapes_resp["zscores"].values())[0])):
        tmp_probe = {"probe_properties": {}}
        for gene in jsonobj['selectedGenes']:
            tmp_probe["probe_properties"][gene] = brainscapes_resp["zscores"][gene][i]

        tmp_probe["position"] = brainscapes_resp["mnicoord"][i]

        if brainscapes_resp["area"][i] == brainscapes_areas[0]:
            print("Area 1")
            probes_area1.append(tmp_probe)
        else:
            print("Area 2")
            probes_area2.append(tmp_probe)

    print(probes_area1)
    print(probes_area2)

    result = {"result": brainscapes_resp["p-values"]}
    #result["Version"] = os.environ["OPENSHIFT_BUILD_COMMIT"]
    result["Areas"] = []
    result["Areas"].append(
                        {
                            "name": jsonobj["area1"]["areas"][0]["name"],
                            "hemisphere": jsonobj["area1"]["areas"][0].get('hemisphere'),
                            "probes": probes_area1
                        })
    result["Areas"].append(
                    {
                        "name": jsonobj["area2"]["areas"][0]["name"],
                        "hemisphere": jsonobj["area2"]["areas"][0].get('hemisphere'),
                        "probes": probes_area2
                    })

    print(result)

    return result


def run_pyjugex_analysis(jsonobj):
    print(jsonobj)
    filter_threshold = jsonobj.get('threshold', default_param['threshold'])
    n_rep = jsonobj.get('nPermutations', default_param['nPermutations'])


    brainscapes.logger.setLevel("INFO") # we want to see some messages!

    atlas = brainscapes.atlases.MULTILEVEL_HUMAN_ATLAS

    if len(jsonobj["area1"]["areas"]) == 1 and len(jsonobj["area2"]["areas"]) == 1:
        area1_julich_brain_version = jsonobj["area1"]["areas"][0]["atlas"]["version"]
        area2_julich_brain_version = jsonobj["area2"]["areas"][0]["atlas"]["version"]

        print(area1_julich_brain_version)

        if not area1_julich_brain_version == area2_julich_brain_version:
            print("version mismatch")

        area1_julich_brain_version.replace(".", "_")
        print(area1_julich_brain_version)

        atlas.select_parcellation(brainscapes.parcellations['julich 2.9'])
        atlas.threshold_continuous_maps(float(filter_threshold))
        jugex = siibra_jugex.DifferentialGeneExpression(atlas)

        for gene in jsonobj['selectedGenes']:
            jugex.add_candidate_genes(gene)

        roi1=jsonobj["area1"]["areas"][0]["name"]
        roi2=jsonobj['area2']["areas"][0]["name"]
        jugex.define_roi1(roi1)
        jugex.define_roi2(roi2)

        jugex.run(permutations=n_rep)
        jugex_result = jugex.result()
        print(jugex_result)

        logger.log("info", {"jugex_result": str(jugex_result)})

        result = _transform_brainscapes_response(jugex_result, jsonobj)

        logger.log("info", {"returned_result": str(result)})

        return json.dumps(result)
    else:
        raise ValueError('Multiple regions are not yet supported')


async def handle_post(request):
    if request.can_read_body:
        jsonobj = await request.json()
    else:
        return web.Response(status=400)
    try:
        data = run_pyjugex_analysis(jsonobj)
        return web.Response(status=200,content_type="application/json",body=data)
    except Exception as e:
        logger.log('error', {"error":str(e)})
        return web.Response(status=400,body=str(e))

async def return_auto_complete(request):
     return web.Response(status=200,content_type="application/json",body=dictAutocompleteString)

async def handle_post2(request):
    if request.can_read_body:
        jsonobj = await request.json()
    else:
        return web.Response(status=400)
    if "cbUrl" in jsonobj:
        web.Response(status=200)
        try:
            data = run_pyjugex_analysis(jsonobj)
            requests.post(jsonobj["cbUrl"], json=json.loads(data))
        except Exception as e:
            error = {
                "error": str(e),
                "detail": {
                    "jsonobj": jsonobj
                }
            }
            logger.log('error', error)
            requests.post(jsonobj["cbUrl"], json=error)
    else:
        try:
            data = run_pyjugex_analysis(jsonobj)
            return web.Response(status=200, content_type="application/json", body=data)
        except Exception as e:
            logger.log('error', {"error":str(e)})
            return web.Response(status=400,body=str(e))

def main():
    app = web.Application()
    cors = aiohttp_cors.setup(app)
    cors.add(app.router.add_post("/jugex",handle_post), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})

    cors.add(app.router.add_post("/jugex_v2", handle_post2), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
    cors.add(app.router.add_get("/",return_auto_complete), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
    logger.log('info', {"message": "webjugex backend started"})
    web.run_app(app,host="0.0.0.0",port=8003)

if __name__=='__main__':
    main()
