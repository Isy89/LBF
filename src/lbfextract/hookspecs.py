"""
In this module we introduced the specification for both the CLI hooks as well as the LBFextract hooks, which are the 
ones responsible for the feature extraction workflow.

"""

import logging
import pathlib
from typing import Any, Union, List

import click
import matplotlib.figure
import pandas as pd
import pluggy

from lbfextract.utils_classes import Signal

logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

hookspec = pluggy.HookspecMarker("lbfextract")
hookspec_cli = pluggy.HookspecMarker("lbfextract_cli")


class FextractHooksSpecs:
    @hookspec(firstresult=True)
    def save_fetched_reads(self, reads_per_interval_container: pd.DataFrame,
                           config: Any,
                           extra_config: Any
                           ) -> None:
        """
        Hook implementing the strategy to save the reads fetched for the intervals
        :param reads_per_interval_container: pd.DataFrame containing the following columns:
                'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'reads_per_interval'. 
            The reads_per_interval columns contains for each interval the pysam.libcalignment.IterRowRegion object 
            from which the reads are fetched or a least of reads
        :param config: configuration file containing the configuration object required by the save_fetched_reads
                       function
        :param extra_config: extra configuration that may be used in the hook implementation

        :return: None
        """
        pass

    @hookspec(firstresult=True)
    def load_fetched_reads(self, config: Any, extra_config: Any) -> pd.DataFrame:
        """
        Hook implementing the strategy to load the reads fetched for the intervals
        :param config: configuration file containing the configuration object required by the load_fetched_reads
                       function
        :param extra_config: extra configuration that may be used in the hook implementation

        :return: pd.DataFrame containing the following columns:
                'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'reads_per_interval'. 
            The reads_per_interval columns contains for each interval the pysam.libcalignment.IterRowRegion object 
            from which the reads are fetched or a least of reads
        """
        pass

    @hookspec(firstresult=True)
    def transform_reads(self, reads_per_interval_container: pd.DataFrame, config: Any,
                        extra_config: Any):
        """
        :param pd.DataFrame containing the following columns:
                'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'reads_per_interval'. 
            The reads_per_interval columns contains for each interval the pysam.libcalignment.IterRowRegion object 
            from which the reads are fetched or a least of reads
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        pass

    @hookspec(firstresult=True)
    def transform_single_intervals(self, transformed_reads: pd.DataFrame, config: Any,
                                   extra_config: Any):
        """
        :param transformed_reads: pd.DataFrame containing the following columns:
                'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'reads_per_interval'. 
            The reads_per_interval columns contains for each interval the pysam.libcalignment.IterRowRegion object 
            from which the reads are fetched or a least of reads
        :param config: config specific to the function
        :param extra_config: config containing context information plus extra parameters
        """
        pass

    @hookspec(firstresult=True)
    def transform_all_intervals(self, single_intervals_transformed_reads: Signal, config: Any, extra_config: Any):
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        pass

    @hookspec(firstresult=True)
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
        :return: pd.DataFrame containing the following columns:
                'Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'reads_per_interval'. 
            The reads_per_interval columns contains for each interval the pysam.libcalignment.IterRowRegion object 
            from which the reads are fetched or a least of reads
        """
        pass

    @hookspec(firstresult=True)
    def plot_signal(self, signal: Signal, config: Any, extra_config: Any) -> matplotlib.figure.Figure:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        pass

    @hookspec(firstresult=True)
    def save_signal(self,
                    signal: Signal,
                    config: Any,
                    extra_config: Any) -> None:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        pass


class CliHookSpec:
    @hookspec_cli
    def get_command(self) -> Union[click.Command, List[click.Command]]:
        pass
