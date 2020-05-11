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

import logging
import requests

class HBPLogger:
  def __init__(self, url=None, application_name=None, deployment='local'):
    self.url = url
    self.application_name = application_name
    self.deployment = deployment

  def log(self, tag, payload):
    if self.url is None:
      print(payload)
    else:
      r = requests.post(url='{url}{application_name}.{deployment}.{tag}'.format(url=self.url, deployment=self.deployment, application_name=self.application_name, tag=tag), json=payload)
