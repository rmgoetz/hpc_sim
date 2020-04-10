# -*- coding: utf-8 -*-
from setuptools import setup

setup(name='hpc',
      version='0.1',
      description='Simulation package for heterodyne phase camera',
      long_description=open('README.txt').read(),
      install_requires=['numpy','pykat','matplotlib'],
      classifiers=[
            'Development Status :: 1 - Planning',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering :: Physics'    
      ],
      url='https://github.com/Mo-Rice/hpc_sim',
      author='Maurico Diaz-Ortiz',
      author_email='mdiazort@ufl.edu',
      license='MIT',
      packages=['hpc'],
      zip_safe=False)

