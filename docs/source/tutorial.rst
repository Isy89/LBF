LBFextract tutorial
-------------------

In this tutorial, we show how to use the CLI and python interfaces provided by LBFextract. In LBFextract, we provide 
fragmentomics-based and coverage-based feature extraction methods. In this tutorial we will focus on transcription 
factors (TFs).
This will be done for CTCF and AR, but can be done for any TF for which a BED file is provided.  A TF is 
represented by a BED file containing regions corresponding to its TFBSs as provided by the GTRD database. 
When a TF binds to the chromatin it may cause nucleosome displacement and a chromatin conformation change.
This results in a TF footprint on the chromatin, which can be observed from liquid biopsy data. In liquid biopsy, we have cfDNA,
which is the DNA shed into body fluids by cells undergoing apoptosis, necrosis, necroptosis or NETosis. DNA molecules bound to
nucleosomes have higher protection against degradation. In regulatory regions, TF binding  
displaces nucleosomes in a TFBS  and causes a phasing effect of nucleosomes within regions bound by the same TF within 
a cell and between cells.
This causes a signal to emerge when these regions are analyzed from coverage prospective or fragmentation perspective.
Here we will see how these information can be extracted from liquid biopsy data using LBFextract.

1. Installation
---------------
To run this tutorial, it is suggested to create a folder in which the code and the results will be stored and set this as
the current working directory. 

Open a terminal and run the following command:

.. code-block:: bash
    
    mkdir -p tutorial_lbf && cd tutorial_lbf

To run this tutorial, conda is required. If you do not have conda installed, the following can be used to install it:

.. code-block:: bash
    
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
    rm Miniconda3-latest-Linux-x86_64.sh

After the download is complete and conda is installed, the prompt should show a `(base)` prefix. Running `conda init`
may be necessary to initialize conda.
A new conda environment for this tutorial can be created and activated as follow:

.. code-block:: bash

    conda create -n LBF_tutorial python=3.10 -y
    conda activate LBF_tutorial

Inside this environment we will install the necessary packages required
for the tutorial (jupyter, ipykernel):

.. code-block:: bash

    python -m pip install jupyter
    python -m pip install ipykernel
    python -m ipykernel install --name LBF_tutorial --user

We will now clone the LBF repository and enter the folder containing the code (LBF folder):

.. code-block:: bash

    git clone https://github.com/Isy89/LBF.git && cd LBF

To install the LBFextract package the following command can be run:

.. code-block:: bash

    python -m pip install .

To be sure that the installation succeded, we can run the following command:

.. code-block:: bash

    lbfextract --help

This should show the help message of the CLI.
To install lbfextract specific conda environment to have a separate running samtools with the correct version, run:

.. code-block:: bash

    lbfextract setup create-conda-envs

this step can be skipped if samtools is already installed in the system and the correct version is available.

LBF extract can be used in 3 ways:

1. Through the command line interface (CLI) 
2. Through the terminal user interface (TUI) 
3. Through the python interface

We should be now inside the LBF folder and the conda environment LBF_tutorial should be active. This can be seen by the 
prompt being preceded by (LBF_tutorial).

2. CLI inerface 
---------------

2.1. Help

After running the command `lbfextract --help`, in the help message we can see the available subcommands:

1. feature_extraction_commands
2. post_extraction_analysis_commands
3. setup
4. start-tui

to have a look at the available feature extraction methods one can run:

.. code-block:: bash

    lbfextract feature_extraction_commands --help

This will show the help message of the feature_extraction_commands subcommand. The feature_extraction_commands subcommand 
are reported with a truncated help message. To see for example the full explanation of the extract-coverage subcommand
one can run the following command:

.. code-block:: bash

    lbfextract feature_extraction_commands extract-coverage --help

2.2. Extract coverage

We will now extract coverage from the sample provided in the test dataset folder.
To do this, we need to provide the following arguments to the `lbfextract feature_extraction_commands extract-coverage` 
command:

1. --path_to_bam: the path to the BAM file containing the reads
2. --path_to_bed: the path to the BED file containing the regions of interest
3. --output_path: the path to the output file

We will use the datasets provided in the `tests/test_datasets` folder. We will use the CTCF and AR BED files and the 
fextract_anonymized_test.bam file. 

To extract coverage for CTCF, the following command can be run:

.. code-block:: bash

    lbfextract feature_extraction_commands extract-coverage \
    --path_to_bam "tests/test_dataset/bam/fextract_anonymized_test.bam" \
    --path_to_bed "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed"\
    --output_path "tests/test_out/test_cli" \
    --summarization_method mean \
    --cores 4

The result will be a pickle (.pkl ) file containing the Signal object, which contains:

- the coverage signal in the array attribute, 
- the metadata in the metadata attribute
- the tags in the tags attribute and a pdf containing the plot

The coverage signal is a 1D numpy vector containing the coverage signal for the regions of interest. In case the value of 
the summarization method was set to skipp the array will be a 2D matrix.
We can extract also the same signal from more than one BED file using the feature extraction methods containing an 
``in-batch` suffix.

To extract coverage for CTCF and AR, run the following command:

.. code-block:: bash

    lbfextract feature_extraction_commands extract-coverage-in-batch \
    --path_to_bam "tests/test_dataset/bam/fextract_anonymized_test.bam" \
    --path_to_bed "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/" \
    --output_path "tests/test_out/test_cli" \
    --summarization_method mean \
    --cores 4

The output will contain:

1. a plot for each BED file
2. a batch_signals plot containing the top and bottom 5 signals one next to the other 
3. a csv file containing one signal per BED file in each row 
4. a pkl file containing the Signal object
5. a correlation matrix 
6. the summary containing:
    - the amplitude for the top and bottom 5 signals
    - the correlation for the top and bottom 5 signals
    - the heatmap for the top and bottom 5 signals
    - the PCA plot for all BED files


2.3. Extract fragmentomics-based features

The fragment length distribution per position can be extracted in the following way:

.. code-block:: bash

    lbfextract feature_extraction_commands extract-fragment-length-distribution \
    --path_to_bam "tests/test_dataset/bam/fextract_anonymized_test.bam" \
    --path_to_bed "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed" \
    --output_path "tests/test_out/test_cli" \
    --cores 4

The results is a plot showing a heatmap of the distribution per position. Each column corresponds to a position relative
to the center of the binding site. The rows correspond to the fragment length.
At the bottom the averaged fragment length distribution is shown for chunks of position over the +- n positions considered 
(+-2000 by default). The same can be done using the in-batch version of this feature extraction method

.. code-block:: bash

    lbfextract feature_extraction_commands extract-fragment-length-distribution-in-batch \
    --path_to_bam "tests/test_dataset/bam/fextract_anonymized_test.bam" \
    --path_to_bed "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/" \
    --output_path "tests/test_out/test_cli" \
    --cores 4

The output in this case will be a plot of the fragment length distribution per position for each BED file and a .npz file 
containing each array in a different key.

3. Python interface
-------------------
To better enjoy the python interface, we suggest to use a jupyter notebook. To start a jupyter notebook or lab, run the following:

.. code-block:: bash

    jupyter notebook

or:

.. code-block:: bash

    jupyter lab

This will open a new tab in the browser. In the new tab, navigate to the folder where the tutorial is located and create a
new notebook. In the notebook, one can copy paste the following commands.
    
The python interface can be used to extract the same features as the CLI.

To use the python interface the `FeatureExtractor` object needs to be imported from the `lbfextract.feature_extractor` 
module.

.. code-block:: python

    import lbfextract
    from lbfextract.feature_extractor import FeatureExtractor
    import pathlib
    
    fe = FeatureExtractor()

we can now extract the coverage signal for CTCF using the following code:

.. code-block:: python

    path_to_bam = "tests/test_dataset/bam/fextract_anonymized_test.bam"
    path_to_bed = "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed"
    output_path = "tests/test_out/test_python"
    summarization_method = "mean"
    
    fe.extract(
        "extract_coverage",
         path_to_bam = pathlib.Path(path_to_bam),
         path_to_bed = pathlib.Path(path_to_bed),
         output_path = pathlib.Path(output_path),
         summarization_method = summarization_method,
         cores = 4
    )

The same can be done for the in-batch version of the feature extraction method:

.. code-block:: python

    path_to_bam = "tests/test_dataset/bam/fextract_anonymized_test.bam"
    path_to_bed = "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/"
    output_path = "tests/test_out/test_python"
    summarization_method = "mean"
    
    fe.extract(
        "extract_coverage_in_batch",
         path_to_bam = pathlib.Path(path_to_bam),
         path_to_bed = pathlib.Path(path_to_bed),
         output_path = pathlib.Path(output_path),
         summarization_method = summarization_method,
         cores = 4
    )

The fragment length distribution can be extracted using the following code:

.. code-block:: python

    path_to_bam = "tests/test_dataset/bam/fextract_anonymized_test.bam"
    path_to_bed = "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed"
    output_path = "tests/test_out/test_python"
    
    fe.extract(
        "extract_fragment_length_distribution",
         path_to_bam = pathlib.Path(path_to_bam),
         path_to_bed = pathlib.Path(path_to_bed),
         output_path = pathlib.Path(output_path),
         cores = 4
    )


The same can be done for the in-batch version of the feature extraction method:

.. code-block:: python

    path_to_bam = "tests/test_dataset/bam/fextract_anonymized_test.bam"
    path_to_bed = "tests/test_dataset/multi_bed_ar_ctcf_for_dyads/"
    output_path = "tests/test_out/test_python"
    
    fe.extract(
        "extract_fragment_length_distribution_in_batch",
         path_to_bam = pathlib.Path(path_to_bam),
         path_to_bed = pathlib.Path(path_to_bed),
         output_path = pathlib.Path(output_path),
         cores = 4
    )


4. Differentially active Transcription Factors (TFs)
----------------------------------------------------

In this section we will show how to use the CLI and python interface to extract differentially active TFs. We will use
a generated dataset containing 2 groups of samples having 20% of differentially active TFs.

4.1. Generate the dataset

We can generate the dataset using the `DatasetGenerator` object. This can be done in the following way:

.. code-block:: python

    from lbfextract.data.dummy_dataset_generator import DummyDatasetGenerator
    
    dg = DummyDatasetGenerator(n_tf=100, groups=2, dimensions=2000, n_samples=20, n_diff_active=0.2, noise_level=1,sign=-1)
    dataset = dg.create_dataset()
    dg.save_dataset_by_sample(dataset, output_path="tests/test_dataset/generated_dataset")
    
    
This creates a datasets having two groups having 20 samples each having 100 TF with signals having 4000 dimensions 
with 20% of differentially active TFs. Each sample is saved as a csv file in the generated_dataset folder.

4.2. Generation of the sample sheet

To run the rest of the analysis we can go back to the previous terminal. Be sure to be in the `LBF` folder and the
`LBF_tutorial` conda environment is active.
LBFextract offers a way to automatically generate a sample sheet containing a column with the paths to each generated
.csv file and a empty column (group), which will contain the group membership of each sample.
This is done in the following way: 

.. code-block:: bash

    lbfextract post_extraction_analysis_commands \
    --path_to_res_summary "tests/test_dataset/generated_dataset" \
    --output_path "tests/test_out/" \
    generate-sample-sheet

Once the sample sheet is generated the information about the group membership can be filled in and the information about the 
samples name can be corrected. When sample specific datasets are generated through a normal lbfextract feature extraction method
they will be imported with the correct name. This may need to be changed when the structure of the results of the feature 
extraction method dose not follow the standard structure.

We can attach the group information to the sample sheet using the following command:

.. code-block:: python

    import pandas as pd
    import pathlib
    
    df = pd.read_csv("tests/test_out/fextract_diff_signal_results/sample_sheet.csv", index_col=0)
    df["sample_name"] = df["paths_to_sample_result"].apply(lambda x: pathlib.Path(x).stem)
    df["group"] = df["sample_name"].apply(lambda x: x.split("_")[0])
    df.to_csv("tests/test_out/fextract_diff_signal_results/sample_sheet_filled.csv")

4.3. Extract differentially active TFs

Once the groups are given and the sample name corrected we can start the differentially active analysis in the following way:

.. code-block:: bash

    lbfextract post_extraction_analysis_commands \
    --save_indivitual_plots \
    --commit_hash "TEST" \
    --center_signal_indices "1985,2015" \
    --remove_outliars \
    --flanking_signal_indices "1000,3000" \
    --correction_method "fdr_bh" \
    --max_iter 1 \
    --path_to_res_summary "tests/test_dataset/generated_dataset" \
    --alpha 0.05 \
    --output_path "tests/diff_active_tfs" \
    --path_to_sample_sheet "tests/test_out/fextract_diff_signal_results/sample_sheet_filled.csv" \
    --outer_group_column "group" \
    get-differentially-active-genomic-intervals

The output will be located under `tests/diff_active_tfs/`. The output will contain:

- a folder containing the heatmap of the differentially active transcription factors per group
- a barplot of the number of differentially active TFs per group
- a subfolder for each pair of groups containing the plots for each differentially active TFs, the table of the enrichment
  analysis
- the table of the differentially active transcription factors.
- a .pkl file containing the results as python object
- a metadata.json file containing the metadata of the analysis
