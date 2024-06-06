# Persistence-capable sheaf manipulation library setup script
#
# Copyright (c) 2017, Michael Robinson
# This software comes with no warrantees express or implied

from setuptools import setup
#from Cython.Build import cythonize

def readme():
    with open('README') as f:
        return f.read()

setup(name='pysheaf',
      version='0.3.1',
      description='Python applied sheaf computation library',
      url='http://github.com/kb1dds/pysheaf',
      author='Michael Robinson and Brenda Praggastis and Janelle Henrich',
      author_email='michaelr@american.edu',
      license='MIT',
      packages=['pysheaf'],
      install_requires=['scipy >= 1.0','numpy','networkx','deap','pdoc'],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Mathematics',
      ],
#      ext_modules=cythonize('pysheaf/pysheaf.py'),
      zip_safe=False)
