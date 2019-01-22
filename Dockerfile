FROM python:3

COPY . /webjugex
WORKDIR /webjugex

RUN pip install -r requirements.txt
RUN pip install .

WORKDIR /webjugex/webjugex

ENTRYPOINT [ "python", "aioserver.py" ]