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

FROM python:3.7

COPY . /webjugex
WORKDIR /webjugex

RUN pip install -r requirements.txt
RUN pip install .

RUN git clone https://jugit.fz-juelich.de/v.marcenko/brainscapes.git
RUN cd brainscapes && git checkout koehnen_Dictionary_merge_compatability\
  && pip install -e .


WORKDIR /webjugex/webjugex

EXPOSE 8003

ENTRYPOINT [ "python", "aioserver.py" ]
