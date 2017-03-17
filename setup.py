from setuptools import find_packages, setup

with open('README.rst') as f:
    readme=f.read()

setup(name='dbotu',
    version='1.1',
    description='Distribution-based OTU calling',
    long_description=readme,
    author='Scott Olesen',
    author_email='swo@alum.mit.edu',
    platforms=['any'],
    license='License :: OSI Approved :: MIT License',
    url='http://github.com/swo/dbotu3',
    scripts=['dbotu.py'],
    package_date={'dbotu': ['data/input/counts.txt', 'data/input/seq.fa',
        'data/output/log.txt', 'data/output/membership.txt', 'data/output/otu.txt']},
    classifiers=['Development Status :: 3', 'Environment :: Console',
        'Programming Language :: Python :: 3', 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages=find_packages()
    )
