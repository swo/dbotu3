from setuptools import find_packages, setup
import pypandoc

with open('LICENSE') as f:
    license=f.read()

setup(name='dbotu',
        version='1.0',
        description='Distribution-based OTU calling',
        long_description=pypandoc.convert('README.md', 'rst'),
        author='Scott Olesen',
        author_email='swo@alum.mit.edu',
        platforms=['any'],
        license=license,
        url='http://github.com/swo/dbotu3',
        data_files=[('input', ['data/input/counts.txt', 'data/input/seq.fa']),
            ('output', ['data/output/log.txt', 'data/output/membership.txt', 'data/output/otu.txt'])],
        packages=find_packages(),
        )
