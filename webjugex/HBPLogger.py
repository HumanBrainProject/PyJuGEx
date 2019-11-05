import logging
import requests

class HBPLogger:
  def __init__(self, url=None, application_name=None, deployment='local'):
    self.url = url
    self.application_name = application_name

  def log(self, tag, payload):
    if self.url is None:
      print(payload)
    else:
      r = requests.post(url='{url}{application_name}.{deployment}.{tag}'.format(url=self.url, application_name=self.application_name, tag=tag), json=payload)