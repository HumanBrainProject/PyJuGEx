from distutils.core import setup
setup(name='pyjugex',
      version='0.1',
      description='Perform differential gene expression on two chosen brain regions' ,
      author='Haimasree Bhattacharya',
      author_email='h.bhattacharya@fz-juelich.de',
      py_modules=['pyjugex', 'analysispyjugex', 'hbp_human_atlas'],
      install_requires=['numpy', 'scipy', 'statsmodels', 'requests', 'nibabel', 'xmltodict',], 
      )
