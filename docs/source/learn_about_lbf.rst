Learn about LBFextract
======================
introduction to lbf
----------------
LBFextract is a Python package for extracting features from a a bam file,
with a special focus on liquid biopsy related features and transcription factors. 
The package is built as a plugin interface, where each plugin is a feature.
It is composed by a core package, which contains the main logic, and a set of
plugins which are the features. The core package is responsible for describing the
workflow and how different hooks will be 
executed to extract the features. The plugins implement the hooks.

.. image:: _static/LBF_structure.png
    :alt: LBF hook system and plugins architecture

The current hooks that can be implemented there are:

* ***fetch_reads***: extract the feature from a bam file
* ***load_reads***: load reads in case they were already extracted
* ***save_fetched_reads***: save the fetched reads specific to the regions of interest
* ***transform_reads***: apply a transformation to each read extracted
* ***transform_single_intervals***: extract the signal of one region
* ***transform_all_intervals***: apply a transformation which requires all the regions
* ***plot_signal***: plot the final signal
* ***save_signal***: save the final signal

LBFextracts provides also CLIhooks which allow the automatic integration of all 
the plugins with LBFextract cli and tui interfaces.

installation
------------
For the installation of LBFextract, the following is required:

- python>=3.10
- build (`python -m pip install build`)
- conda 
- setuptools~=62.0.0

conda is used to create a separate environment used by LBFextract, such that the user conda environment is not influenced. If samtools is already available, conda is not necessary. When LBFextract does not find its specific conda environment, it will look for samtools in the current environment. Be aware that samtools should be version 1.14.

To be able to run the tests, the following Python package is also required:

- pytest~=8.1.1

LBFextract can be installed as follows:

.. code-block:: bash

    git clone https://github.com/Isy89/LBF.git && cd LBF
    python -m build --wheel . && python -m pip install --force-reinstall dist/*.whl # "python -m pip install ." should also work

After the installation, the command line interface `lbfextract` should be available. Using it, a conda environment isolated from the current one containing samtools can be created. This step can be omitted if samtools~=1.14 is already present. This installation of this conda env is done as follows:

.. code-block:: bash

    lbfextract setup create-conda-envs # creates a separate conda env used for filtering the bam files and other steps


Singularity Image isntallation
-------------------------------

To install LBFextract using the Singularity image, the following steps are required:
.. code-block:: bash

    singularity pull lbfextract_v0.1.0a1.sif library://lbfextract/lbfextract/lbfextract_v0.1.0a1.sif:0.1.0a1
    singularity run lbfextract_v0.1.0a1.sif --help

Using the run command you will have access to the lbfextract command line interface.
When using the singularity image it may be necessary to bind the directory containing the BAM files and BED files and
the output directory to the singularity container. This can be done using the following command:

.. code-block:: bash

    singularity run --bind /path/to/data_bam:/data_bam --bind /path/to/data_bed:/data_bed --bind /path/to/output_dir:/output_dir lbfextract_v0.1.0a1.sif --help

example:

.. code-block:: bash

    singularity run --bind /path/to/data_bam:/data_bam --bind /path/to/data_bed:/data_bed --bind /path/to/output_dir:/output_dir lbfextract_v0.1.0a1.sif feature_extraction_commands extract-coverage --path_to_bam /data_bam/example.bam --path_to_bed /data_bed/example.bed --output_path /output_dir



Coming Soon: Installation via pip (PyPI)
-----------------------------------------

We are currently working on making LBFextract installable directly from the Python Package Index (PyPI) using pip. This feature will allow for easier installation and distribution across different platforms.

Stay tuned for updates on when this feature will be available. In the meantime, please refer to the installation instructions provided above.



usage
-----

LBFextract can be used through the command line interface (CLI), through the
terminal user interface (TUI) or through the python API.

The CLI offers four major set of commands:

1. feature_extraction_commands
2. post_extraction_analysis_commands
3. setup
4. start-tui

The first set of commands are used to extract the features from the bam file.
The second set of commands are used to analyze the extracted features.
The third set of commands are used to setup the conda environments required
for the features present in LBFextract to work.
The fourth command is used to start the TUI interface.