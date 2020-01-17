from setuptools import setup, find_packages
import os

name = 'familyanalyzer'

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name=name,
    version='0.6.0',
    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='A tool to analyse gene family evolution from orthoxml',
    long_description=read('README.rst'),
    license='MIT',
    classifiers = [
         'Development Status :: 3 - Alpha',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'License :: OSI Approved :: MIT licence',
         'Programming Language :: Python :: 2',
         'Programming Language :: Python :: 2.7',
         'Programming Language :: Python :: 3',
         'Programming Language :: Python :: 3.3',
         'Programming Language :: Python :: 3.4',
         'Programming Language :: Python :: 3.5',
         ],
    packages=find_packages(),
    install_requires=['lxml', 'progressbar-latest', 'future'],
    scripts=['bin/familyanalyzer']
)
