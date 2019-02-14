from aiohttp import web
import aiohttp_cors
import json
import nibabel as nib
import hbp_human_atlas as atlas
import webjugex
import os

with open("files/genesymbols.txt", "r") as f:
    dictAutocompleteString = f.read()
#dictAutocomplete = ["ADRA2A", "AVPR1B", "CHRM2", "CNR1", "CREB1", "CRH", "CRHR1", "CRHR2", "GAD2", "HTR1A", "HTR1B", "HTR1D", "HTR2A", "HTR3A", "HTR5A", "MAOA", "PDE1A", "SLC6A2", "SLC6A4", "SST", "TAC1", "TPH1", "GPR50", "CUX2", "TPH2"]

# get cache dir from environment variable
if os.getenv(['GENE_CACHE_DIR']) is not None:
    gene_cache_dir = os.getenv['GENE_CACHE_DIR']
else:
    gene_cache_dir = '.pyjugex'

def pyjugex_analysis(jsonobj):
    roi1 = {}
    roi2 = {}
    roi1['data'] = atlas.jubrain.probability_map(jsonobj['area1']['PMapURL'], jsonobj['area1']['name'], atlas.MNI152)
    roi1['name'] = jsonobj['area1']['name']
    roi2['data'] = atlas.jubrain.probability_map(jsonobj['area2']['PMapURL'], jsonobj['area2']['name'], atlas.MNI152)
    roi2['name'] = jsonobj['area2']['name']
    jugex = webjugex.Analysis(gene_cache_dir=gene_cache_dir, filter_threshold=jsonobj['threshold'], single_probe_mode = jsonobj['mode'], verbose=True)
    result = jugex.DifferentialAnalysis(jsonobj['selectedGenes'], roi1, roi2)
    return result

async def handle_post(request):
    if request.can_read_body:
        jsonobj = await request.json()
    else:
        return web.Response(status=400)
    try:
        data = pyjugex_analysis(jsonobj)
        return web.Response(status=200,content_type="application/json",body=data)
    except Exception as e:
        print(e)
        return web.Response(status=400,body=str(e))

async def return_auto_complete(request):
     return web.Response(status=200,content_type="application/json",body=dictAutocompleteString)


def main():
    app = web.Application()
    cors = aiohttp_cors.setup(app)
    cors.add(app.router.add_post("/jugex",handle_post), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
    cors.add(app.router.add_get("/",return_auto_complete), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
    cors.add(app.router.add_static('/public/',path=str('./public/')))
    web.run_app(app,host="0.0.0.0",port=8003)

if __name__=='__main__':
    main()
