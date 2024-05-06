<div align="center">
	<img src="docs/source/_static/logo.png">
</div>

----

[![Code of Conduct](https://img.shields.io/badge/code%20of-conduct-ff69b4.svg?style=flat)](CODE_OF_CONDUCT.md)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![Read the Docs](https://readthedocs.org/projects/yt2mp3/badge/?version=latest)](https://lbf.readthedocs.io/)
![python versions](https://img.shields.io/badge/python->=3.10-blue.svg)

## New: singularity image installation 

have a look at the documentation for the singularity installation [here](https://lbf.readthedocs.io/)

## Introduction.


LBFextract is a Python package for extracting features for all genomic intervals described in a Browser Extensible Data (BED) file or multiple BED files, from a Binary Alignment Map (BAM) file and identifying condition-specific or cluster-specific differentially active Transcription Factors (TF).
It focuses on liquid biopsy related features, transcription factor binding sites (TFBSs) and Transcription Start Sites (TSSs), but can be generalized to any kind of genomic intervals with similar properties. 
The package is built as a plugin interface, in which each plugin is a feature. It is composed by a core package, which contains the main logic, and a set of
plugins, which represent the features extraction methods. The core package (lbfextract) describes the workflow and how different hooks will be executed to extract the features. 
The plugins implement the hooks. Default coverage-based and fragmentoimics-based feature extraction methods are provided as lbfextract subpackages. 

The following feature extraction methods are available:

- coverage
- coverage-in-batch
- central-60b (Peter Ulz coverage)
- sliding-window-coverage
- sliding-window-coverage-in-batch
- wps-coverage
- coverage-around-dyads
- coverage-around-dyads-in-batch
- middle-point-coverage
- middle-point-coverage-in-batch
- middle-n-points-coverage
- middle-n-points-coverage-in-batch
- entropy
- entropy-in-batch 
- fragment-length-distribution ( per position )
- fragment-length-distribution-in-batch ( per position )
- fragment-length-ratios ( per position )
- relative-entropy-to-flanking
- relative-entropy-to-flanking-in-batch
- extract-signal

## Installation

For the installation of LBFextract, the following is required:

    - python>=3.10
    - build (python -m pip install build)
    - conda 
    - setuptools~=62.0.0

conda is used to create a separate environment used by LBFextract, such that the user conda environment is 
not influenced. If samtools is already available, conda is not necessary. When LBFextract does not find its specific 
conda environment, will look for samtools in the current env. Be aware that samtools should be 1.14.

To be able to run the tests also the following python package is required:

    - pytest~=8.1.1



LBFextract can be installed as follows:

```bash
git clone https://github.com/Isy89/LBF.git && cd LBF
python -m build --wheel . && python -m pip install --force-reinstall dist/*.whl # "python -m pip install ." should also work
```

After the installation, the command line interface lbfextract should be available.
Using it, a conda environment isolated from the current one containing samtools can be created.
This step can be omitted if samtools~=1.14 is already present.
This installation of this conda env is done as follows:

```bash
lbfextract setup create-conda-envs # creates a separate conda env used for filtering the bam files and other steps
```

if using bash, autocompletion can be enabled (at the moment it only support autocompletion for bash, other shell will be added in the feature...):

```bash
lbfextract setup enable-autocompletion
```

and disabled:

```bash
lbfextract setup disable-autocompletion
```

## Documentation

[https://lbf.readthedocs.io/en/latest/](https://lbf.readthedocs.io/en/latest/)


## Copyright

Original work on lbfextract package and subpackages accessory code Copyright (c) 2023 Isaac Lazzeri

## Licence

GNU General Public License v3.0

## Contact

For any questions please contact:

* <LBFextract@gmail.com>

If you find any bugs please report them here:

* <https://github.com/Isy89/LBF/issues> 
