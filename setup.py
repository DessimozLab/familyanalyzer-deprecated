from setuptools import setup, find_packages
import os

name = 'familyanalyzer'

with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as fh:
    long_description = fh.read()


setup(
    name=name,
    version='0.7.2',
    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='A tool to analyse gene family evolution from orthoxml',
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license='MIT',
    classifiers = [
         'Development Status :: 3 - Alpha',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'License :: OSI Approved :: MIT License',
         'Programming Language :: Python :: 2',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3',
         'Programming Language :: Python :: 3.5',
         'Programming Language :: Python :: 3.6',
         'Programming Language :: Python :: 3.7',
         'Programming Language :: Python :: 3.8',
         ],
    packages=find_packages(),
    install_requires=['lxml', 'progressbar-latest', 'future'],
    scripts=['bin/familyanalyzer']
)
