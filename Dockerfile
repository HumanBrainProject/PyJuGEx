FROM python:3.7

COPY . /webjugex
WORKDIR /webjugex

RUN pip install -r requirements.txt
RUN pip install .

WORKDIR /webjugex/webjugex

EXPOSE 8003

ENTRYPOINT [ "python", "aioserver.py" ]
