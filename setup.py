# Copyright 2020 Human Brain Project/EBRAINS
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

from setuptools import setup
setup(name='pyjugex',
      version='1.0.1alpha',
      packages=['pyjugex'],
      license='apache-2.0',
      description='Perform web based differential gene expression on two chosen brain regions',
      url='https://github.com/HumanBrainProject/PyJuGEx',
      author='Big Data Analytics Group, INM-1, Research Center Juelich',
      author_email='inm1-bda@fz-juelich.de',
      install_requires=[
            'numpy',
            'scipy==1.2',
            'statsmodels',
            'requests',
            'nibabel',
            'xmltodict'
      ])
