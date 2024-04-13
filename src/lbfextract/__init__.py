"""

LBFextract
==========

LBFextract is a Python package designed for extracting features from a BAM file, with a particular emphasis on liquid 
biopsy-related features and transcription factors. The package employs a plugin interface where each plugin represents a
feature that can be extracted from a BAM file. It consists of a core package housing the main logic, supplemented by a 
collection of plugins encompassing various features. Within the ``lbfextract`` package, the core components include:

- **Core Module:** Describes the workflow and underlying logic.
- **CLI Module:** Constructs the Command-Line Interface (CLI).
- **Feature Extractor:** Builds the Python interface.
- **Hookspecs:** Defines the interface for hooks.
- **Plugin Manager:** Manages the loading of plugins.

Third-party plugins can implement the following hooks:

- **fetch_reads:** Extracts features from a BAM file.
- **load_reads:** Loads reads if previously extracted.
- **save_fetched_reads:** Saves the fetched reads specific to regions of interest.
- **transform_reads:** Applies transformations to each extracted read.
- **transform_single_intervals:** Extracts the signal of a single region.
- **transform_all_intervals:** Applies transformations requiring data from all regions.
- **plot_signal:** Plots the final signal.
- **save_signal:** Saves the final signal.

LBFextract utilizes Pluggy to provide an interface for third-party software. The general interface for hooks, located 
in the ``hookspecs.py`` file, outlines the expected signatures for plugins. Default hook implementations are found in 
``lib.py`` within the ``fextract`` subpackage.

Additionally, CLI-specific hooks are specified in ``hookspecs.py``, with default implementations found in ``cli_lib.py``
and other subpackage plugin modules. These hooks supply CLI commands to execute workflows associated with each plugin.

For seamless integration with Python, the ``FeatureExtractor`` object allows users to invoke different feature 
extraction methods programmatically. Upon installing new plugins containing ``cli_hooks``, the associated commands are
automatically incorporated into the ``FeatureExtractor`` object as accessible methods.

To maintain consistency, each plugin should include a ``schema.py`` file defining the configuration schema for its 
feature extraction method. Refer to existing plugins, particularly those generating default coverage and 
fragmentomics-based signals, for guidance. Schemas should adhere to the Voluptuous schema syntax. 
(`voluptuous documentation <https://pypi.org/project/voluptuous/>`_).


"""

import logging.config
import os
import pathlib
import subprocess
import tempfile
import click
import pluggy
import yaml
from rich.traceback import install

from lbfextract.setup_conda_env import check_in_conda_env, check_if_conda_env_is_base
from lbfextract.utils_classes import ExternalModuleFilter

install(show_locals=True, suppress=[click, pluggy])

__version__ = "0.1.0a1"

os.environ['MPLCONFIGDIR'] = os.path.join(tempfile.gettempdir(), "matplotlib")

with open(pathlib.Path(__file__).parent / "config" / 'logger_config.yml', 'rt') as f:
    logger_config = yaml.safe_load(f.read())

with open(pathlib.Path(__file__).parent / "config" / 'profiler_config.yml', 'rt') as f:
    PROFILER_DEBUG = yaml.safe_load(f.read())["debug"]

# Configure the logging system
logging.config.dictConfig(logger_config)
root_logger = logging.getLogger()

if PROFILER_DEBUG:
    root_logger.setLevel(logging.DEBUG)

# attaching the filter to all handlers to filter out messages not coming from lbfextract
external_module_filter = ExternalModuleFilter()
for handler in root_logger.handlers:
    handler.addFilter(external_module_filter)
logger = logging.getLogger(__name__)

hookimpl = pluggy.HookimplMarker("lbfextract")
hookimpl_cli = pluggy.HookimplMarker("lbfextract_cli")

if check_in_conda_env():
    if check_if_conda_env_is_base():

        env_samtools = (
                pathlib.Path(os.environ.get("CONDA_EXE")).parent.parent /
                "envs" /
                "@lbfextract" /
                "bin" /
                "samtools"
        )
    else:
        if (pathlib.Path(os.environ.get("CONDA_PREFIX")) / "bin" / "samtools").exists():
            logging.warning(f"found samtools in current env, but samtools from fextract "
                            f"conda env will be used for functions of the fextract package ")
        env_samtools = (
                pathlib.Path(os.environ.get("CONDA_EXE")).parent.parent / "envs" /
                "@lbfextract" /
                "bin" /
                "samtools"
        )
    if not env_samtools.exists() and not (
            pathlib.Path(os.environ.get("CONDA_PREFIX")) / "bin" / "samtools").exists():
        logging.warning("@lbfextract conda env was not generated. run 'setup_conda_env' to create it")
        logging.warning(
            "neither was samtools found in the @lbfextract conda env nor a samtools is available"
            "in the current conda env")
        logging.warning("No samtools available! lbfextract functions using samtools will fail!"
                        "run setup_conda_env to fix this!")
        env_samtools = None
else:
    try:
        path_to_samtools = subprocess.run(["which", "samtools"], check=True, capture_output=True).stdout.decode()
        logging.warning(f"using samtools at: {path_to_samtools}")
        env_samtools = pathlib.Path(path_to_samtools)

    except subprocess.CalledProcessError:
        logging.warning("No samtools available! lbfextract functions using samtools will fail!"
                        "run setup_conda_env to fix this!")
        env_samtools = None

if env_samtools:
    PATH_TO_SAMTOOLS = str(env_samtools)
else:
    logger.warning("samtools file was not found. Be sure to have run lbfextract setup create-conda-envs and that the env"
                   "fextract is present. Or install samtools version 1.14. ")
