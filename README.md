```
         __       _______  ________                     __                                  __     
        |  \     |       \|        \                   |  \                                |  \    
        | ▓▓     | ▓▓▓▓▓▓▓\ ▓▓▓▓▓▓▓▓ ______  __    __ _| ▓▓_    ______   ______   _______ _| ▓▓_   
        | ▓▓     | ▓▓__/ ▓▓ ▓▓__    /      \|  \  /  \   ▓▓ \  /      \ |      \ /       \   ▓▓ \  
        | ▓▓     | ▓▓    ▓▓ ▓▓  \  |  ▓▓▓▓▓▓\\▓▓\/  ▓▓\▓▓▓▓▓▓ |  ▓▓▓▓▓▓\ \▓▓▓▓▓▓\  ▓▓▓▓▓▓▓\▓▓▓▓▓▓  
        | ▓▓     | ▓▓▓▓▓▓▓\ ▓▓▓▓▓  | ▓▓    ▓▓ >▓▓  ▓▓  | ▓▓ __| ▓▓   \▓▓/      ▓▓ ▓▓       | ▓▓ __ 
        | ▓▓_____| ▓▓__/ ▓▓ ▓▓     | ▓▓▓▓▓▓▓▓/  ▓▓▓▓\  | ▓▓|  \ ▓▓     |  ▓▓▓▓▓▓▓ ▓▓_____  | ▓▓|  \
        | ▓▓     \ ▓▓    ▓▓ ▓▓      \▓▓     \  ▓▓ \▓▓\  \▓▓  ▓▓ ▓▓      \▓▓    ▓▓\▓▓     \  \▓▓  ▓▓
         \▓▓▓▓▓▓▓▓\▓▓▓▓▓▓▓ \▓▓       \▓▓▓▓▓▓▓\▓▓   \▓▓   \▓▓▓▓ \▓▓       \▓▓▓▓▓▓▓ \▓▓▓▓▓▓▓   \▓▓▓▓ 
    
        O       o O         o O       o
        | O   o | | O     o | | O   o |
        | | O | | | | LBF | | | | O | |
        | o   O | | o     O | | o   O |
        o       O o         O o       O
        
        For any questions please contact: 
        - <LBFextract@gmail.com>

        If you find any bugs please report them here:
        - <https://github.com/Isy89/LBF/issues>

        Copyright: Original work on LBFextract and accessory code Copyright (c) 2023 Isaac Lazzeri
        Licence: GNU General Public License v3.0
```

## A plugin implementation of feature extraction from fastq files and bed files.

Fextract defines a series of hooks to carry out the feature extraction process from bam files. It extracts signals from
the intervals defined in the bed files. The fextract default package calculates coverage across intervals defined in the
bed file. Plugins for fextract can be implemented defining the hooks that fextract exposes as entry points during the
program runtime

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
