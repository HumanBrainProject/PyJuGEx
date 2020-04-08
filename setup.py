from setuptools import setup
setup(name='pyjugex',
      version='1.0.1alpha',
      packages=['pyjugex'],
      license='BSD',
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