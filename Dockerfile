FROM python:3

COPY . /webjugex
WORKDIR /webjugex

RUN pip install -r requirements.txt
RUN pip install .

RUN groupadd -r appuser && useradd --no-log-init -r -g appuser appuser
USER appuser

WORKDIR /webjugex/webjugex

EXPOSE 8003

ENTRYPOINT [ "python", "aioserver.py" ]
