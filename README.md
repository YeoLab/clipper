# CLIPper - CLIP peak enrichment recognition

A tool to detect CLIP-seq peaks.

Please visit our wiki page to learn more about usage of clipper: https://github.com/YeoLab/clipper/wiki/CLIPper-Home

## Installation

```shell script
# recreate conda environment
cd clipper
conda env create -f environment.yml
python setup.py install
```
## Command Line Usage

```shell script
clipper -h
````

## Run test
```shell script
cd clipper/clipper/test
python -m unittest discover
```
Right now the test coverage is still not 100%.
And `test_spline.py` don't pass. 

## Reference
Yeo GW, Coufal NG, Liang TY, Peng GE, Fu XD, Gage FH. An RNA code for the FOX2 splicing regulator revealed by mapping RNA-protein interactions in stem cells. Nat Struct Mol Biol. 2009;16(2):130-137. doi:10.1038/nsmb.1545

Lovci MT, Ghanem D, Marr H, et al. Rbfox proteins regulate alternative mRNA splicing through evolutionarily conserved RNA bridges. Nat Struct Mol Biol. 2013;20(12):1434-1442. doi:10.1038/nsmb.2699