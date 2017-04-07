# Persistence-capable sheaf manipulation library setup script
#
# Copyright (c) 2017, Michael Robinson
# This software comes with no warrantees express or implied 

from setuptools import setup

def readme():
    with open('README') as f:
        return f.read()
    
setup(name='pysheaf',
      version='0.1',
      description='Python applied sheaf computation library',
      url='http://github.com/kb1dds/pysheaf',
      author='Michael Robinson',
      author_email='michaelr@american.edu',
      license='MIT',
      packages=['pysheaf'],
      install_requires=['numpy','matplotlib','networkx','scipy'],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Mathematics',
      ],
      zip_safe=False)
