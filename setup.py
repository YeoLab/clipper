from setuptools import setup, find_packages, Extension

peaks = Extension("src/peaks", sources = ['src/peaksmodule.cc'],
#for debugging
#                  extra_compile_args = ['-O0'] 
)                 

with open("README") as file:
    long_description = file.read()

setup(
    name = "clipper",
    long_description = long_description,
    version = "0.1",
    packages = find_packages(),
    ext_modules = [peaks],
    package_data = {
        '' : ['*.lengths', '*.gz', '*.bam', '*.bai']
        },
    
    install_requires = ['setuptools >= 0.6', 
                        'pysam >= 0.6',
                        'numpy >= 1.5.1 ',
                        'scipy >= 0.8.0',
                        'matplotlib >= 1.1.0',
                        'pp >= 1.6.2',
                        'pybedtools >= 0.6',
                        ],
      
    setup_requires = ["setuptools_git >= 0.3",],
    
    entry_points = {
                    'console_scripts': [
                                        'clipper = src.peakfinder:call_main',],
                    },
    #metadata for upload to PyPI
    author = "Michael Lovci and Gabriel Pratt",
    author_email = "mlovci@ucsd.edu",
    description = "A set of scripts for calling peaks on CLIP-seq data",
    license = "GPL2",
    keywords = "CLIP-seq, peaks, bioinformatics",
    url = "https://github.com/YeoLab/clipper",
    
    #Other stuff I feel like including here
    include_package_data = True,
    #zip_safe = True #True I think
)
