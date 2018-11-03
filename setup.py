
import sys
import shutil
from subprocess import call
from setuptools import setup

if sys.version_info.major != 3:
    raise RuntimeError('SCTree test requires Python 3')

setup(name='SCTree',
      description='SCTree test can statistical test the continuous tree sturecture for single-cell data',
      author='Xiangqi Bai',
      author_email='xqbai@amss.ac.cn',
      package_dir={'': 'src'},
      packages=['SCTree'],
      install_requires=[
          'numpy>=1.12.0',
          'pandas>=0.19.2',
          'scipy>=0.18.1',
          'Cython',
          'bhtsne',
          'matplotlib>=2.0.0',
          'seaborn>=0.7.1',
          'sklearn',
          'networkx>=1.11',
          'fcsparser>=0.1.2',
          'statsmodels>=0.8.0'],
      )

# install wishbone
if shutil.which('pip3'):
    call(['pip3', 'install', 'git+https://github.com/manusetty/wishbone.git'])

