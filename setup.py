#!/usr/bin/env python3
from glob import glob
from distutils.core import setup

setup(name='dct-redux',
      version='0.2.0',
      description='DCT/LMI reduction and analysis scripts.',
      author='Michael S. P. Kelley',
      author_email='msk@astro.umd.edu',
      url='https://github.com/mkelley/dct-redux',
      scripts=glob('bin/*.py') + glob('bin/*.sh')
      )
