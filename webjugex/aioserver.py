from aiohttp import web
import aiohttp_cors
import json
import nibabel as nib
import hbp_human_atlas as atlas
import webjugex
import os
import requests

from default import default_param


with open("files/genesymbols.txt", "r") as f:
    dictAutocompleteString = f.read()
#dictAutocomplete = ["ADRA2A", "AVPR1B", "CHRM2", "CNR1", "CREB1", "CRH", "CRHR1", "CRHR2", "GAD2", "HTR1A", "HTR1B", "HTR1D", "HTR2A", "HTR3A", "HTR5A", "MAOA", "PDE1A", "SLC6A2", "SLC6A4", "SST", "TAC1", "TPH1", "GPR50", "CUX2", "TPH2"]

# get cache dir from environment variable
if os.getenv('GENE_CACHE_DIR') is not None:
    gene_cache_dir = os.getenv('GENE_CACHE_DIR')
else:
    gene_cache_dir = '.pyjugex'

def get_roi_img_array(obj):
    pmap_resp = webjugex.util.get_pmap(obj['PMapURL'], body=obj['body'] if 'body' in obj else None)
    return webjugex.util.read_byte_via_nib(pmap_resp.content, gzip=webjugex.util.is_gzipped(obj['PMapURL']))

def run_pyjugex_analysis(jsonobj):
    roi1 = {}
    roi2 = {}

    roi1_obj = jsonobj['area1']
    roi1['data'] = get_roi_img_array(roi1_obj)
    roi1['name'] = jsonobj['area1']['name']

    roi2_obj = jsonobj['area2']
    roi2['data'] = get_roi_img_array(roi2_obj)
    roi2['name'] = jsonobj['area2']['name']

    single_probe_mode = jsonobj.get('mode', default_param['mode'])
    filter_threshold = jsonobj.get('threshold', default_param['threshold'])
    n_rep = jsonobj.get('nPermutations', default_param['nPermutations'])

    jugex = webjugex.Analysis(gene_cache_dir=gene_cache_dir, filter_threshold=filter_threshold, single_probe_mode = single_probe_mode, verbose=True, n_rep=n_rep)

    result = jugex.DifferentialAnalysis(jsonobj['selectedGenes'], roi1, roi2)
    return result

async def handle_post(request):
    if request.can_read_body:
        jsonobj = await request.json()
    else:
        return web.Response(status=400)
    try:
        data = run_pyjugex_analysis(jsonobj)
        return web.Response(status=200,content_type="application/json",body=data)
    except Exception as e:
        print(e)
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
            error = {}
            error['error'] = e
            requests.post(jsonobj["cbUrl"], json=json.loads(error))
            print("result callback error", e)
    else:
        try:
            data = run_pyjugex_analysis(jsonobj)
            return web.Response(status=200,content_type="application/json",body=data)
        except Exception as e:
            print(e)
            return web.Response(status=400,body=str(e))

def main():
    app = web.Application()
    cors = aiohttp_cors.setup(app)
    cors.add(app.router.add_post("/jugex",handle_post), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})

    cors.add(app.router.add_post("/jugex_v2", handle_post2), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
    cors.add(app.router.add_get("/",return_auto_complete), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
    cors.add(app.router.add_static('/public/',path=str('./public/')))
    web.run_app(app,host="0.0.0.0",port=8003)

if __name__=='__main__':
    main()
