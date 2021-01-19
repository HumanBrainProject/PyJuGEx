import time
import json
import requests
import threading

from os import environ
from string import Template

REFRESH_TOKEN_env = 'JUGEX_REFRESH_TOKEN'
CLIENT_ID_env = 'JUGEX_CLIENT_ID'
CLIENT_SECRET_env = 'JUGEX_CLIENT_SECRET'
HBP_OIDC_ENDPOINT_env = 'HBP_OIDC_ENDPOINT'

def build_request_object():
    request_template = Template('grant_type=refresh_token&refresh_token=${REFRESH_TOKEN}&client_id=${CLIENT_ID}&client_secret=${CLIENT_SECRET}')
    result = request_template.substitute(
                REFRESH_TOKEN = environ[REFRESH_TOKEN_env],
                CLIENT_ID = environ[CLIENT_ID_env],
                CLIENT_SECRET = environ[CLIENT_SECRET_env]
            )
    return result

def request_token():
    result = requests.post(
                environ[HBP_OIDC_ENDPOINT_env],
                data = build_request_object(),
                headers = {'content-type': 'application/x-www-form-urlencoded'}
            )
    return result

class jwt_handler:
    def _update_token(self):
        while True:
            self.token = json.loads(request_token().content.decode("utf-8"))
            time_until_next_refresh = int(self.token["expires_in"]) - 300
            time.sleep(time_until_next_refresh)

    def __init__(self):
        self.token = None
        update_thread = threading.Thread(target=self._update_token, name="update_token")
        update_thread.start()

def main():
    token_handler = jwt_handler()
    time.sleep(5)
    print(token_handler.token)

if __name__ == "__main__":
    main()
