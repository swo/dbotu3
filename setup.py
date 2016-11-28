from setuptools import find_packages, setup
import pypandoc

with open('LICENSE') as f:
    license=f.read()

setup(name='dbotu',
    version='1.1',
    description='Distribution-based OTU calling',
    long_description=pypandoc.convert('README.md', 'rst'),
    author='Scott Olesen',
    author_email='swo@alum.mit.edu',
    platforms=['any'],
    license=license,
    url='http://github.com/swo/dbotu3',
    data_files=[('input', ['data/input/counts.txt', 'data/input/seq.fa']),
        ('output', ['data/output/log.txt', 'data/output/membership.txt', 'data/output/otu.txt'])],
    classifiers=['Development Status :: 3', 'Environment :: Console',
        'Programming Language :: Python :: 3', 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages=find_packages()
    )
