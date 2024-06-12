"""
cli_lib module
==============

In this module the commands for the feature extraction methods based on coverage are described. 
The following commands are included:

* ***extract_signal***:
    allows the combination of different packages, but requires the specification of the classes and
    configurations required by the different hooks.
* *** extract_coverage***:
    extract the coverage from all the genomic intervals contains in a BED file and report this summarized using
    a summarization method like mean or median at each position.
* *** extract_wps_coverage***:
    extract the windowed protection score over multiple genomic intervals, summarizes the signal and
    returns the normalized aggregated wps score per position
* *** extract_middle_point_coverage***:
    extracts midpoint coverage at each position of multiple genomic intervals contained in a
    BED file and summarizes results per position.
* *** extract_middle_n_points_coverage***:
    extract the n points from the left and right of the middle point of each fragment and
    calculate coverage for each genomic interval in a BED file taking only these position for each fragment into account
* *** extract_sliding_window_coverage***:
    extract a windowed average of the coverage calculating at each position the average of the
    coverage values of the n-th following positions.
* *** extract_peter_ulz_coverage (central-60 coverage)***:
    calculates central 60 base read coverage using the same coordinates utilized by Peter Ulz et al
    in its `Nature communication <https://www.nature.com/articles/s41467-019-12714-4>`_
* *** extract_coverage_dyads***:
    extract the coverage considering only the position of each fragment coming from where the dyads
    are probably located.
"""

from __future__ import annotations

import logging
import pathlib
import pickle
from typing import Optional, Union, List

import click
import numpy as np
import pysam
import scipy
import yaml
from scipy.signal import savgol_filter

import lbfextract.fextract
from lbfextract.core import App
from lbfextract.fextract.schemas import Config, AppExtraConfig, SignalSummarizer, SingleSignalTransformerConfig, \
    ReadFetcherConfig

logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def calculate_reference_distribution(path_to_sample, min_length, max_length, chr, start, end):
    alignment_file = pysam.AlignmentFile(path_to_sample, "rb")
    reads = alignment_file.fetch(chr, start, end)
    array_fragment_lengths = np.zeros(max_length - min_length)
    for i in reads:
        if i.tlen < min_length or i.tlen >= max_length:
            continue
        array_fragment_lengths[i.tlen - min_length] += 1

    return array_fragment_lengths / array_fragment_lengths.sum()


def get_peaks(distribution, height=0.2, distance=100):
    distribution = savgol_filter(distribution, 10, 3)
    peaks = scipy.signal.find_peaks(distribution, height=height, distance=distance)[0]
    return peaks


def open_yml(ctx, param, value) -> dict:
    with open(value, "r") as f:
        content = yaml.load(f, Loader=yaml.FullLoader)
    return content


def load_checks(ctx, param, value: tuple):
    checks = []
    for check_to_load in value:
        with open(check_to_load, "rb") as f:
            checks.append(pickle.load(f))
    return tuple(checks)


class CliHook:

    @lbfextract.hookimpl_cli
    def get_command(self) -> click.Command:
        @click.command()
        @click.option('--plugins_name', default=None, multiple=True, help='Names of the plugin to register.')
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
        @click.option("--skip_read_fetching", is_flag=True)
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option('--read_fetcher_config', default={}, callback=open_yml,
                      help='path to the yml file containing the configuration used in the fetch_reads hook')
        @click.option('--reads_transformer_config', default={}, callback=open_yml,
                      help='path to the yml file containing the configuration used in the transform_reads hook')
        @click.option('--single_signal_transformer_config', default={}, callback=open_yml,
                      help='path to the yml file containing the configuration used in the transform_single_intervals '
                           'hook')
        @click.option('--transform_all_intervals_config', default={}, callback=open_yml,
                      help='path to the yml file containing the configuration used in the transform_all_intervals hook')
        @click.option('--plot_signal_config', default={}, callback=open_yml,
                      help='path to the yml file containing the configuration used in the plot_signal hook')
        @click.option('--save_signal_config', default={}, callback=open_yml,
                      help='path to the yml file containing the configuration used in the save_signal hook')
        @click.option('--extra_config', default={}, callback=open_yml,
                      help='path to the yml file containing extra configuration')
        def extract_signal(plugins_name: list, path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                           output_path: pathlib.Path,
                           skip_read_fetching, reads_transformer_config,
                           read_fetcher_config: dict, single_signal_transformer_config: dict,
                           transform_all_intervals_config: dict,
                           plot_signal_config: dict, save_signal_config: dict, extra_config: dict,
                           exp_id: Optional[str]):
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)

            res = App(plugins_name=plugins_name,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=Config(single_signal_transformer_config),
                      transform_all_intervals_config=Config(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

        return extract_signal


class CliHookExtractCoverage:
    @lbfextract.hookimpl_cli
    def get_command(self) -> Union[click.Command, List[click.Command]]:
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
                      help='path to the BAM file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the BED file to be used')
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
                      help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest '
                           'timestamp of this run (identified by the id) is retrieved')
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option("--window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted around the middle point of an"
                           " interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper pairs, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a BAM file')
        def extract_coverage(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
            Given a set of genomic intervals having the same length w this feature extraction method extracts the 
            fragment coverage from each genomic interval and summarizes the resulting signals with one of the summarization
            methods provided: mean, median, max, min, skip. In case skip is selected then, no summary is generated and the 
            coverage value at each position for all genomic intervals is reported.
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
            res = App(plugins_name=None,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

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
                           " interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper pairs, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        @click.option('--window_size', type=int, default=120, show_default=True,
                      help='window size to be used to for the wps signal')
        @click.option('--min_fragment_length', type=int, default=121, show_default=True, help='minimum fragment length')
        @click.option('--max_fragment_length', type=int, default=180, show_default=True, help='maximum fragment length')
        def extract_wps_coverage(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                 gc_correction_tag: Optional[str],
                                 window_size: int,
                                 min_fragment_length: int,
                                 max_fragment_length: int
                                 ):
            """
            Given a set of genomic intervals having the same length w, extract_wps_coverage feature extraction method 
            extracts the windowed protection score at each position.
            """
            if min_fragment_length <= window_size:
                raise ValueError("min_fragment_length must be greater than window_size to avoid reads starting and "
                                 "ending in the same wps window")

            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "wps_coverage",
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag,
                "window_size": window_size,
                "min_fragment_length": min_fragment_length,
                "max_fragment_length": max_fragment_length,

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
            res = App(plugins_name=None,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

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
                      help='path to the BAM file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the BED file to be used')
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
                      help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest '
                           'timestamp of this run (identified by the id) is retrieved')
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option("--window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted around the middle point of an "
                           "interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper pairs, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--n", default=1, type=int, show_default=True,
                      help="number of bases to retain around the dyad")
        @click.option("--cores", type=int, show_default=True, default=1,
                      help="number of cores to be used for the computation")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        def extract_coverage_dyads(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                                   output_path: pathlib.Path,
                                   skip_read_fetching,
                                   window: int,
                                   flanking_window: int,
                                   extra_bases: int,
                                   n_binding_sites: int,
                                   summarization_method: str,
                                   n: int,
                                   cores: int,
                                   exp_id: Optional[str],
                                   flip_based_on_strand: bool,
                                   gc_correction_tag: Optional[str]):
            """
            Given a set of genomic intervals having the same length w this feature extraction method extracts the 
            fragment coverage around nucleosomes dyads. To achieve this, it models from which poly-nucleosomal structure
            a fragment might have originated from and infers the position of the dyad before degradation accordingly.
            """

            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "n": n,
                "signal_transformer": "coverage_dyads",
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag
            }
            distribution = calculate_reference_distribution(path_to_sample=path_to_bam,
                                                            min_length=0,
                                                            max_length=600,
                                                            chr="chr12",
                                                            start=34_300_000,
                                                            end=34_500_000
                                                            )
            peaks = get_peaks(distribution, height=0.01, distance=100)
            single_signal_transformer_config["peaks"] = [peaks[0]]
            transform_all_intervals_config = {
                "summarization_method": summarization_method,

            }
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=None,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

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
                      help='path to the BAM file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the BED file to be used')
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
                      help="Integer describing the number of bases to be extracted around the middle point of an "
                           "interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper paires, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min, skip }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a BAM file')
        def extract_middle_point_coverage(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                          gc_correction_tag: Optional[str]):
            """
            This command extracts the middle point coverage at each position for all provided genomic intervals having 
            the same length. This can be summarized with the following methods: mean, median, max, min, skip. In case skip 
            is selected no summary is generated and the profiles of each genomic interval are reported. The midpoint 
            coverage corresponds at each position to the number of fragments' midpoints overlapping that position. 
             """

            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "middle_point_coverage",
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag
            }
            transform_all_intervals_config = {
                "summarization_method": summarization_method,
            }
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=None,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

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
                      help='path to the BAM file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the BED file to be used')
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
                      help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest '
                           'timestamp of this run (identified by the id) is retrieved')
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option("--window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted around the middle point of an "
                           "interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper paires, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min, skip }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--n_middle_pos", default=10, type=int, show_default=True,
                      help="number of position around the middle point to be extracted")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        def extract_middle_n_points_coverage(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                                             output_path: pathlib.Path,
                                             skip_read_fetching,
                                             window: int,
                                             flanking_window: int,
                                             extra_bases: int,
                                             n_binding_sites: int,
                                             n_middle_pos: int,
                                             summarization_method: str,
                                             cores: int,
                                             exp_id: Optional[str],
                                             flip_based_on_strand: bool,
                                             gc_correction_tag: Optional[str]
                                             ):
            """
             This command extracts the middle n points coverage at each position for all provided genomic intervals having 
             the same length. This can be summarized with the following methods: mean, median, max, min, skip. In case skip
             is selected no summary is generated and the profiles of each genomic interval are reported. The middle n points 
             coverage corresponds at each position to the number of fragmetnts' positions lying within the ragne [m-n, m+n) 
             (in which m is the midpoint and n the number of position surrounding it from both sides) overlapping a genomic
             interval
             """

            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "middle_n_points_coverage",
                "n": n_middle_pos,
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag
            }
            transform_all_intervals_config = {
                "summarization_method": summarization_method,
            }
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=None,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

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
                      help='path to the BAM file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the BED file to be used')
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
                      help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest '
                           'timestamp of this run (identified by the id) is retrieved')
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option("--window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted around the middle point of an "
                           "interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper paires, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min, skip }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--window_size", default=5, type=int, show_default=True,
                      help="window size to be used to compute the sliding_window_coverage")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        def extract_sliding_window_coverage(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                                            output_path: pathlib.Path,
                                            skip_read_fetching,
                                            window: int,
                                            flanking_window: int,
                                            extra_bases: int,
                                            n_binding_sites: int,
                                            summarization_method: str,
                                            cores: int,
                                            exp_id: Optional[str],
                                            window_size: int,
                                            flip_based_on_strand: bool,
                                            gc_correction_tag: Optional[str]):
            """
            The sliding window coverage corresponds to a windowed average of the coverage signal
            """
            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "sliding_window_coverage",
                "window_size": window_size,
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag
            }
            transform_all_intervals_config = {
                "summarization_method": summarization_method,

            }
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=None,
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

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
                      help='path to the BAM file to be used')
        @click.option('--path_to_bed', type=click.Path(exists=False,
                                                       file_okay=True,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the BED file to be used')
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
                      help='Boolean flag. When it is set, the fetching of the reads is skipped and the latest '
                           'timestamp of this run (identified by the id) is retrieved')
        @click.option("--exp_id", default=None, type=str, show_default=True,
                      help="run id")
        @click.option("--window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted around the middle point of an"
                           " interval present in the BED file")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the BAM file when removing the "
                           "unused bases to be sure to get all the proper pairs, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min", "skip"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True,
                      show_default=True,
                      default=False,
                      help="flip the signal based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a BAM file')
        @click.option("--read_start", default=53, type=int, show_default=True,
                      help="start of the read to be used to extract coverage")
        @click.option("--read_end", default=113, type=int, show_default=True,
                      help="end of the read to be used to extract coverage")
        def extract_peter_ulz_coverage(
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
                gc_correction_tag: Optional[str],
                read_start: int,
                read_end: int,
        ):
            """
            Given a set of genomic intervals having the same length w, this feature extraction method extracts the read 
            coverage in which only the central position of a read are taken into acount. which part of the read is taken 
            into account is defined by the paramenters, which by default are 53-113
            """
            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "peter_ulz_coverage",
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag,
                "read_start": read_start,
                "read_end": read_end
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
            res = App(plugins_name=None,
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

        return [extract_coverage, extract_coverage_dyads, extract_middle_point_coverage,
                extract_middle_point_coverage, extract_middle_n_points_coverage, extract_sliding_window_coverage,
                extract_peter_ulz_coverage, extract_wps_coverage]
