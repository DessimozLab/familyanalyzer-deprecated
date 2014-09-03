from setuptools import setup, find_packages

name = 'familyanalyzer'

setup(
    name=name,
    version='0.5.0',
    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='todoc',
    packages=find_packages(),
    install_requires=['lxml', 'progressbar-latest', 'future'],
    scripts=['bin/familyanalyzer']
)
