from collections import namedtuple
import os
import json

_script_dir = os.path.dirname(__file__)

with open(os.path.join(_script_dir, 'data/genesymbols.json')) as handle:
  _file_content = handle.read()
  _json_array = json.loads(_file_content)

_Allen = namedtuple('Allen', _json_array)

allen = _Allen()