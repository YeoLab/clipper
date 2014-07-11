from setuptools import setup
from setuptools import find_packages

from distutils.extension import Extension
from Cython.Distutils import build_ext
peaks = Extension("clipper.src.peaks", sources = ['clipper/src/peaksmodule.cc'],
#for debugging
#                  extra_compile_args = ['-O0'] 
)                 

readsToWiggle = Extension("clipper.src.readsToWiggle", ['clipper/src/readsToWiggle.pyx'])

long_description = "CLIPPER - clip peak enrichment"
setup(
    name = "clipper",
    long_description = long_description,
    version = "0.2.0",
    packages = find_packages(),
    cmdclass = {'build_ext' : build_ext},
    ext_modules = [readsToWiggle, peaks],

    package_data = {
        'clipper' : ['data/*.gff', 'data/regions/*.bed']
        },

    install_requires = ['setuptools', 
                        'pysam >= 0.6',
                        'numpy >= 1.5.1 ',
                        'scipy >= 0.11.0',
                        'matplotlib >= 1.1.0',
                        'pybedtools >= 0.5',
                        'scikit-learn >= 0.13.0',
                        ],
      
    setup_requires = ["setuptools_git >= 0.3",],
    
    entry_points = {
                    'console_scripts': [
                                        'clipper = clipper.src.peakfinder:call_main',
                                        'clip_analysis = clipper.src.CLIP_analysis:call_main',],
                    },

    #metadata for upload to PyPI
    author = "Michael Lovci and Gabriel Pratt",
    author_email = "mlovci@ucsd.edu",
    description = "A set of scripts for calling peaks on CLIP-seq data",
    license = "GPL2",
    keywords = "CLIP-seq, peaks, bioinformatics",
    url = "https://github.com/YeoLab/clipper",
    
    #Other stuff I feel like including here
    #include_package_data = True
    #zip_safe = True #True I think
)
