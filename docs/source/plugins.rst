Plugins 
=======

.. warning::

    **work in progress**

Plugins offer a way to exchange the behaviour of some part of LBFextract's workflow. 
This is achieved through a hook mechanism. LBFextract workflow for a feature extraction method is described in the core 
module and is composed by a sequence of steps. Each step is a hook. Hooks are implemented by plugins. The default 
behaviour of LBFextract is itself implemented in the fextract plugin, which provides the implementation of the 
coverage-based methods. 
Each hook is essentially a decorated python function. In LBFextract there are two types of hooks: CliHook and FextractHooks.
The former is responsible for the implementation of the different steps of the workflow. The latter is responsible for the 
CLI commands needed to run the feature extraction method.
The implementation in LBFextract is simplified through the "setup new-plugin" command. i.e. to create
a new plugin extracting a  "cool_signal_name" (the complete implementation of this example can be found 
`here <https://github.com/Isy89/fextract_cool_signal_name/tree/main>`_), one can run the following:

.. code-block:: bash
    
    lbfextract setup new-plugin --name_of_the_signal cool_signal_name --name_of_the_cli_command extract_cool_signal_name --out_dir OUTPUTDIR

This will generate most of the boilerplate code needed to generate a LBFextract's pluing. An LBFextract's plugin includes 
the following files:

.. code-block:: bash

    fextract-cool_signal_name
    ├── setup.py
    ├── src
    │   ├── cool_signal_name
    │   │   ├── __init__.py
    │   └── └── plugin.py
    └── tests

Once it has been generated, the plugin can be installed as follow:

.. code-block:: bash
    
    cd OUTPUTDIR &&  python -m pip install .

after a successful installation, the `extract_cool_signal_name` command will be visible in the
lbfextract feature_extraction_commands --help. 
At this moment the package won't do much since the hooks are not implemented. To implement them
one can open the generated plugin.py file.
Inside, there is the definition of all hooks that can be implemented and the default CLI command. 
By default the CLI command is equal to the extract_coverage command, because this is the original definition of the hooks. 
Depending on the hooks that are implemented by the plugin, the CLI command has to be updated. 

After generating the plugin the plugin.py file will look like this:

.. code-block:: python
    :linenos:
    
    import logging
    import pathlib
    from typing import Any
    from typing import List
    from typing import Optional
    
    import click
    import matplotlib
    import pandas as pd
    from lbfextract.utils_classes import Signal
    
    import lbfextract.fextract
    from lbfextract.core import App
    from lbfextract.fextract.schemas import Config, AppExtraConfig, ReadFetcherConfig, SingleSignalTransformerConfig, \
        SignalSummarizer
    
    logger = logging.getLogger(__name__)
    
    class FextractHooks:
    
        @lbfextract.hookimpl
        def fetch_reads(self,
                        path_to_bam: pathlib.Path,
                        path_to_bed: pathlib.Path,
                        config: Any,
                        extra_config: Any) -> pd.DataFrame:
            """
            :param path_to_bam: path to the bam file
            :param path_to_bed: path to the bed file with the regions to be filtered
            :param config: configuration file containing the configuration object required by the fetch_reads function
            :param extra_config: extra configuration that may be used in the hook implementation
            :return: ReadsPerIntervalContainer object containing all the ReadsPerInterval objects in all the intervals
                     contained in the bed file
            """
            return None
    
        @lbfextract.hookimpl
        def save_fatched_reads(self, reads_per_interval_container: pd.DataFrame,
                               config: Any,
                               extra_config: Any
                               ) -> None:
            """
            Hook implementing the strategy to save the reads fetched from the intervals
            :param reads_per_interval_container: ReadsPerIntervalContainer containing information about the genomic region
                and the reads mapping to it
            :param output_path: path to the location where the data should be stored
            :param id: run id
            :param time_stamp: time stamp
            :param extra_config: extra configuration that may be used in the hook implementation
    
            :return: None
            """
            return None
    
        @lbfextract.hookimpl
        def load_fetched_reads(self, config: Any, extra_config: AppExtraConfig) -> pd.DataFrame:
            """
            :param config: config specific to the function
            :param extra_config: extra configuration that may be used in the hook implementation
            """
            return None
    
        @lbfextract.hookimpl
        def transform_reads(self, reads_per_interval_container: pd.DataFrame, config: Any,
                            extra_config: Any) -> pd.DataFrame:
            """
            :param reads_per_interval_container: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
                basically lists with information about start and end of the interval
            :param config: config specific to the function
            :param extra_config: extra configuration that may be used in the hook implementation
            """
            return None
    
        @lbfextract.hookimpl
        def transform_single_intervals(self, transformed_reads: pd.DataFrame, config: Any,
                                       extra_config: Any) -> Signal:
            """
            :param transformed_reads: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
                basically lists with information about start and end of the interval
            :param config: config specific to the function
            :param extra_config: config containing context information plus extra parameters
            """
            return None
    
        @lbfextract.hookimpl
        def transform_all_intervals(self, single_intervals_transformed_reads: Signal, config: Any,
                                    extra_config: Any) -> Signal:
            """
            :param single_intervals_transformed_reads: Signal object containing the signals per interval
            :param config: config specific to the function
            :param extra_config: extra configuration that may be used in the hook implementation
            """
            return None
    
        @lbfextract.hookimpl
        def plot_signal(self, signal: Signal, config: Any, extra_config: Any) -> matplotlib.figure.Figure:
            """
            :param signal: Signal object containing the signals per interval
            :param extra_config: extra configuration that may be used in the hook implementation
            """
            return None
    
        @lbfextract.hookimpl
        def save_signal(self,
                        signal: Signal,
                        config: Any,
                        extra_config: Any) -> None:
            """
            :param signal: Signal object containing the signals per interval
            :param extra_config: extra configuration that may be used in the hook implementation
            """
    
            return None
    
    
    class CliHook:
        @lbfextract.hookimpl_cli
        def get_command(self) -> click.Command | List[click.Command]:
    
            @click.command()
            @click.option('--path_to_bam', type=click.Path(exists=False,
                                                           file_okay=True,
                                                           dir_okay=True,
                                                           writable=False,
                                                           readable=True,
                                                           resolve_path=False,
                                                           allow_dash=True,
                                                           path_type=pathlib.Path,
                                                           executable=False),
                          help='path to the bam file to be used')
            @click.option('--path_to_bed', type=click.Path(exists=False,
                                                           file_okay=True,
                                                           dir_okay=True,
                                                           writable=False,
                                                           readable=True,
                                                           resolve_path=False,
                                                           allow_dash=True,
                                                           path_type=pathlib.Path,
                                                           executable=False),
                          help='path to the bed file to be used')
            @click.option('--output_path', type=click.Path(exists=False,
                                                           file_okay=False,
                                                           dir_okay=True,
                                                           writable=True,
                                                           readable=True,
                                                           resolve_path=False,
                                                           allow_dash=True,
                                                           path_type=pathlib.Path,
                                                           executable=False),
                          help='path to the output directory')
            @click.option("--skip_read_fetching", is_flag=True, show_default=True,
                          help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest'
                               'timestamp of this run (identified by the id) is retrieved')
            @click.option("--exp_id", default=None, type=str, show_default=True,
                          help="run id")
            @click.option("--window", default=1000, type=int, show_default=True,
                          help="Integer describing the number of bases to be extracted around the middle point of an"
                               " interval present in the bedfile")
            @click.option("--flanking_window", default=1000, type=int, show_default=True,
                          help="Integer describing the number of bases to be extracted after the window")
            @click.option("--extra_bases", default=2000, type=int, show_default=True,
                          help="Integer describing the number of bases to be extracted from the bamfile when removing the "
                               "unused bases to be sure to get all the proper paires, which may be mapping up to 2000 bs")
            @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                          help="number of intervals to be used to extract the signal, if it is higher then the provided"
                               "intervals, all the intervals will be used")
            @click.option("--summarization_method", default="mean",
                          type=click.Choice(["mean", "median", "max", "min", "skip"]),
                          show_default=True,
                          help=f"method to be used to summarize the signal: (Undefined, Undefined, Undefined, Undefined)")
            @click.option("--cores", default=1, type=int, show_default=True,
                          help="number of cores to be used")
            @click.option("--flip_based_on_strand", is_flag=True,
                          show_default=True,
                          default=False,
                          help="flip the signal based on the strand")
            @click.option('--gc_correction_tag', type=str,
                          default=None, help='tag to be used to extract gc coefficient per read from a bam file')
    
            # Here you can add the options you want to be available in the command line (they should match the ones in the
            # function)
            def extract_cool_signal_name(
                 path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                 output_path: pathlib.Path,
                 skip_read_fetching,
                 window: int,
                 flanking_window: int,
                 extra_bases: int,
                 n_binding_sites: int,
                 summarization_method: str,
                 cores: int,
                 exp_id: Optional[str],
                 flip_based_on_strand: bool,
                 gc_correction_tag: Optional[str]
                ):
                """
                here yuo can add the code that will be executed when the command is called.
                You can use the following code to customize the hooks
                each hook receives a dictionary of configuration parameters
                if the hook is implemented in the plugin, the parameters should be added in the correct hook.
                Here you can find the parameters of the default hook (extract_coverage). These can be modified
                depending on which hook you are using.
                Hook specified above will be used in place of the default one. The first hook with a non None return value
                will be used. Therfore the order of the plugins specified in the plugins_name list is important.
                If you are inheriting from multiple signals, be sure to pass the correct parameters to the correct hooks.
                Generally, parameters to be pass to hooks are shown in the schemas.py file of each plugin.
                If you define a new hook for best practice you could inherit from lbfextract.fextract.schemas and specify
                the voluptuous schema for the parameters of the hook.
                """
                read_fetcher_config = {
                    "window": window,
                    "flanking_region_window": flanking_window,
                    "extra_bases": extra_bases,
                    "n_binding_sites": n_binding_sites
                }
                reads_transformer_config = {}
                single_signal_transformer_config = {
                    "signal_transformer": "coverage",
                    "flip_based_on_strand": flip_based_on_strand,
                    "gc_correction": True if gc_correction_tag else False,
                    "tag": gc_correction_tag
    
                }
                transform_all_intervals_config = {
                    "summarization_method": summarization_method
                }
                plot_signal_config = {}
                save_signal_config = {}
                extra_config = {
                    "cores": cores
                }
                output_path.mkdir(parents=True, exist_ok=True)
                output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
                output_path_interval_spec.mkdir(parents=True, exist_ok=True)
                res = App(plugins_name=["cool_signal_name", ],
                          path_to_bam=path_to_bam,
                          path_to_bed=path_to_bed,
                          output_path=output_path_interval_spec or pathlib.Path.cwd(),
                          skip_read_fetching=skip_read_fetching,
                          read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                          reads_transformer_config=Config(reads_transformer_config),
                          single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                          transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                          plot_signal_config=Config(plot_signal_config),
                          save_signal_config=Config(save_signal_config),
                          extra_config=AppExtraConfig(extra_config),
                          id=exp_id).run()
                return res
    
            return extract_cool_signal_name
    
    hook = FextractHooks()
    hook_cli = CliHook()


Assuming the cool signal is just the number of reads in each interval, one could implement this by 
replacing the transform_single_intervals and transform_all_intervals hooks as follow:

.. code-block:: python
    :linenos:
    
    @lbfextract.hookimpl
    def transform_single_intervals(self, transformed_reads: pd.DataFrame, config: Any,
                                   extra_config: Any) -> Signal:
        """
        :param transformed_reads: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
            basically lists with information about start and end of the interval
        :param config: config specific to the function
        :param extra_config: config containing context information plus extra parameters
        """
        # transformed_reads is a pandas DataFrame with the following columns: 
        # 'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'reads_per_interval'
        # and reads_per_interval contains pysam.libcalignmentfile.IteratorRowRegion for each interval

        def count_rads_in_interval(x: pysam.libcalignmentfile.IteratorRowRegion) -> int:
            return len(list(x))

        transformed_reads.index = transformed_reads.Chromosome.astype("str") + "-" + transformed_reads.Start.astype(
        "str") + "-" + transformed_reads.End.astype("str")

        reads_per_interval_pd_series = transformed_reads.reads_per_interval.apply(lambda x: count_rads_in_interval(x))

        return Signal(reads_per_interval_pd_series.values,
                      metadata=reads_per_interval_pd_series.index.copy(),
                      tags=tuple(["cool_signal", ]))

    @lbfextract.hookimpl
    def transform_all_intervals(self, single_intervals_transformed_reads: Signal, config: Any,
                                extra_config: Any) -> Signal:
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        return single_intervals_transformed_reads

In this way when executing:

.. code-block:: bash

    lbfextract feature_extraction_commands extract-cool-signal-name --path_to_bam <path_to_bam> --path_to_bed <path_to_bed> --output_path <output_path>

after fetching the reads, lbfextract will execute the newly implemented hooks in place of the default ones. 
One can use already implemented hooks in other plugins by registering them in the App.plugins_name attribute. 
This accepts a list of plugins. Hooks are parsed from left to right, and only the first not None hook will be run.
i.e. if there are 2 plugins implementing the plot_signal hook, which are registered in the App object as follow: 
App.plugins_name(["plug2", "plug1"]), only plug2 plot_signal hook will be run. New parameters can be provided in form of 
config. Config can be `dict` object but LBFextract provides also a Config object, which offers a way of validating the provided 
parameters with a voluptuous schema. One can subclass the Config class providing a voluptuous schema. One can inspect the 
lbfextract.fextract.schema module for an example. 