from aiohttp import web
import aiohttp_cors
import json
# import jugex

with open("files/genesymbols.txt", "r") as f:
  dictAutocompleteString = f.read()

async def handle_post(request):
  '''
  required params:

  roi1 : url of masked nifti
  roi2 : url of masked nifti
  selectedGenes : list of genes symbols
  singleProbeMode : boolean
  ingoreCustomProbe : boolean
  nPermutations : positive integer
  '''
  if request.has_body:
    body = await request.json()
  else:
    return web.Response(status=400,text="cannot parse body as json")
  try:
    analysis_result = pyjugex.analysis(
      roi1=body['roi1'],
      roi2=body['roi2'],
      selected_genes=body['selectedGenes'],
      single_probe_mode=body['singleProbeMode'],
      ignore_custom_probe=body['ignoreCustomProbe'],
      n_permutations=body['nPermutations'])
    return web.Response(status=200,content_type="application/json",body=analysis_result)
  except Exception as e:
    return web.Response(status=400,body=str(e))


async def return_auto_complete(request):
  return web.Response(status=200,content_type="application/json",body=dictAutocompleteString)

def main():
  app = web.Application()
  cors = aiohttp_cors.setup(app)
  cors.add(app.router.add_post("/jugex", handle_post), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
  cors.add(app.router.add_get("/", return_auto_complete), {"*": aiohttp_cors.ResourceOptions(expose_headers="*", allow_headers="*")})
  cors.add(app.router.add_static('/public/',path=str('./public/')))
  web.run_app(app,host="0.0.0.0",port=8003)
  
if __name__=='__main__':
    main()
