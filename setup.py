from setuptools import find_packages, setup

with open('LICENSE') as f:
    license=f.read()

setup(name='dbotu',
        version='0.1',
        description='Distribution-based OTU calling',
        author='Scott Olesen',
        author_email='swo@mit.edu',
        platforms=['any'],
        license=license,
        url='http://github.com/swo/dbotu2',
        packages=find_packages(),
        )
