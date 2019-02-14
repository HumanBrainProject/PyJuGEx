FROM python:3

ARG GENE_CACHE_DIR
ENV GENE_CACHE_DIR=$GENE_CACHE_DIR

COPY . /webjugex
WORKDIR /webjugex

RUN pip install -r requirements.txt
RUN pip install .

WORKDIR /webjugex/webjugex

EXPOSE 8003

ENTRYPOINT [ "python", "aioserver.py" ]
