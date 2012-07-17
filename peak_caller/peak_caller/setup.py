from distutils.core import Extension
from setuptools import setup, find_packages

peaks = Extension("src/peaks", sources = ['src/peaksmodule.cc'])
setup(
    name = "FindPeaks",
    version = "0.1",
    packages = find_packages(),
    ext_modules = [peaks],
    package_data = {
        '' : ['*.lengths', '*.gz', '*.bam']
        },
    #Project Dependences go here

    #data and scripts go here
    install_requires = ['setuptools >= 0.6', 
                        'pysam >= 0.6',
                        'numpy >= 1.5.1 ',
                        'scipy >= 0.11.0b1',
                        'matplotlib >= 0.99.3',
                        'deap >= 0.8',
                        'pybedtools >= 0.6',
                        ],
    #metadata for upload to PyPI
    author = "Michael Lovci and Gabriel Pratt",
    author_email = "lovci@ucsd.edu",
    description = "A set of scripts for calling peaks on CLIP-seq data",
    license = "TBD",
    keywords = "CLIP-seq, peaks, bioinformatics",
    url = "http://yeolab.ucsd.edu/",
    
    #Other stuff I feel like including here
    include_package_data = True,
    #zip_safe = True #True I think
)
