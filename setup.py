# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 16:00:07 2023

@author: Martin.Mihaylov
"""

from setuptools import setup, find_packages



with open("requirements.txt") as requirement_file:
    requirements = requirement_file.read().split()



setup(
      name='LNOI-Constructor',                 # This is the name of your PyPI-package.
      description="A SCT Group Lumerical Constructor package.",
      author="Martin.Mihaylov",
      version='1.0',                          # Update the version number for new releases
      install_requires=requirements,
      packages=find_packages() 
      
 )