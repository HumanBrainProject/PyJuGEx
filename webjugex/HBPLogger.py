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