"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='arc',
    version='0.1',
    description='ARC (Alkali Rydberg Calculator) ',
    long_description=long_description,
    url='https://github.com/nikolasibalic/ARC-Alkali-Rydberg-Calculator',
    author='',
    author_email='',
    license='',
    classifiers=[],
    keywords=['atom', 'rydberg', 'alkali'],
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
)
