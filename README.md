# CLIPper - CLIP peak enrichment recognition

A tool to detect CLIP-seq peaks.

Please visit our wiki page to learn more about usage of clipper: https://github.com/YeoLab/clipper/wiki/CLIPper-Home

## Installation

```shell script
# recreate PYTHON3 conda environment
cd clipper
conda env create -f environment3.yml
conda activate clipper3
pip install .
```
### Alternative installation
1. Thanks @rekado for making clipper available at GNU Guix `guix install clipper`
2. We notice installation might be failing in some platform. Dockerized clipper is in the `eclip` repository [here](https://github.com/YeoLab/eclip).

## Command Line Usage

```shell script
# shows all the options
clipper -h 

# minimal command
clipper -b YOUR_BAM_FILE.bam -o YOUR_OUT_FILE.bed -s hg19
````

## Run test
```shell script
cd clipper/clipper/test
python -m unittest discover
```
Right now the test coverage is still not 100%.
And some subprocess warnings are not handled.


## Frequently Asked Questions:
1. How do use additional reference genome?
[See here for instructions: Supporting additional species](https://github.com/YeoLab/clipper/wiki/Supporting-additional-species)

2. Where can I use specify the Input bam file?
Currently CLIPper does include input normalization. The input normalization pipeline is in another repository: [Merge Peaks](https://github.com/YeoLab/merge_peaks)

## Questions and suggestions
please open an issue in the repo. or email Charlene `hsher@ucsd.edu`

## Reference
Yeo GW, Coufal NG, Liang TY, Peng GE, Fu XD, Gage FH. An RNA code for the FOX2 splicing regulator revealed by mapping RNA-protein interactions in stem cells. Nat Struct Mol Biol. 2009;16(2):130-137. doi:10.1038/nsmb.1545

Lovci MT, Ghanem D, Marr H, et al. Rbfox proteins regulate alternative mRNA splicing through evolutionarily conserved RNA bridges. Nat Struct Mol Biol. 2013;20(12):1434-1442. doi:10.1038/nsmb.2699
