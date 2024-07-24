from __future__ import annotations

import pathlib
from typing import Any, Optional

import click
import matplotlib
import numpy as np
import pandas as pd
import logging
import pysam
import scipy
from scipy.signal import savgol_filter

import lbfextract.fextract.signal_transformer
import lbfextract.fextract_fragment_length_distribution.signal_summarizers
from lbfextract.core import App
from lbfextract.fextract.schemas import AppExtraConfig, ReadFetcherConfig, Config
from lbfextract.utils import generate_time_stamp, sanitize_file_name
from lbfextract.utils_classes import Signal
from lbfextract.fextract_fragment_length_distribution.schemas import SingleSignalTransformerConfig
from lbfextract.plotting_lib.plotting_functions import plot_fragment_length_distribution

logger = logging.getLogger(__name__)


def calculate_reference_distribution(path_to_sample, min_length, max_length, chr, start, end):
    alignment_file = pysam.AlignmentFile(path_to_sample, "rb")
    reads = alignment_file.fetch(chr, start, end)
    array_fragment_lengths = np.zeros(max_length - min_length)
    for i in reads:
        if i.tlen < min_length or i.tlen >= max_length:
            continue
        array_fragment_lengths[i.tlen - min_length] += 1

    return array_fragment_lengths / array_fragment_lengths.sum()


def get_peaks(distribution, height=0.0001, distance=100):
    distribution = savgol_filter(distribution, 10, 3)
    peaks = scipy.signal.find_peaks(distribution, height=height, distance=distance)[0]
    return peaks


def subsample_fragment_lengths(x, n):
    new_x = np.zeros_like(x)
    subsampled_reads = np.random.choice(np.arange(0, x.shape[0]), size=n,
                                        p=np.where(x > 0, x / x.sum(), 0), replace=True)
    for i in range(x.shape[0]):
        new_x[i] = np.sum(subsampled_reads == i)
    return new_x


def get_position_coefficient(config, array):
    window = (config.ctx["read_fetcher_config"].window + config.ctx["read_fetcher_config"].flanking_region_window)
    return (window * 2) / array.shape[0]


class FextractHooks:

    @lbfextract.hookimpl
    def transform_single_intervals(self, transformed_reads: pd.DataFrame,
                                   config: SingleSignalTransformerConfig,
                                   extra_config: AppExtraConfig) -> Signal:
        """
        :param transformed_reads: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
            basically lists with information about start and end of the interval
        :param config: config specific to the function
        :param extra_config: config containing context information plus extra parameters
        """
        signal_transformers_dict = {
            "fld": {"class": "TfbsFragmentLengthDistribution",
                    "config": {"min_fragment_length": config.min_fragment_length,
                               "max_fragment_length": config.max_fragment_length,
                               "gc_correction": config.gc_correction,
                               "tag": config.tag
                               },
                    "tags": ("fld",)},
            "fld_middle": {
                "class": "TfbsFragmentLengthDistributionMiddlePoint",
                "config": {"min_fragment_length": config.min_fragment_length,
                           "max_fragment_length": config.max_fragment_length,
                           "gc_correction": config.gc_correction,
                           "tag": config.tag
                           },
                "tags": ("fld_middle",)
            },
            "fld_middle_n": {
                "class": "TfbsFragmentLengthDistributionMiddleNPoints",
                "config": {"min_fragment_length": config.min_fragment_length,
                           "max_fragment_length": config.max_fragment_length,
                           "gc_correction": config.gc_correction,
                           "tag": config.tag,
                           "n": config.w
                           },
                "tags": ("fld_middle_n",)
            },
            "fld_dyad": {
                "class": "TfbsFragmentLengthDistributionDyad",
                "config": {"min_fragment_length": config.min_fragment_length,
                           "max_fragment_length": config.max_fragment_length,
                           "gc_correction": config.gc_correction,
                           "tag": config.tag,
                           "n": config.w,
                           "peaks": config.peaks
                           },
                "tags": ("fld_dyad",)
            },
            "fld_peter_ulz": {
                "class": "PeterUlzFragmentLengthDistribution",
                "config": {"min_fragment_length": config.min_fragment_length,
                           "max_fragment_length": config.max_fragment_length,
                           "gc_correction": config.gc_correction,
                           "tag": config.tag,
                           "read_start": config.read_start,
                           "read_end": config.read_end,
                           },
                "tags": ("fld_peter_ulz",)
            }

        }
        fld_extractor_params = signal_transformers_dict[config.signal_transformer]["config"]
        fld_extractor = getattr(lbfextract.fextract_fragment_length_distribution.signal_summarizers,
                                signal_transformers_dict[config.signal_transformer]["class"])(
            **fld_extractor_params)
        tag = signal_transformers_dict[config.signal_transformer]["tags"]
        relative_fragment_len = config.max_fragment_length - config.min_fragment_length
        region_length = transformed_reads["End"][0] - transformed_reads["Start"][0]
        tensor = np.zeros((relative_fragment_len, region_length))
        for interval in transformed_reads.itertuples():
            if config.flip_based_on_strand:
                if interval.Strand == "+":
                    tensor += fld_extractor(interval)
                elif interval.Strand == "-":
                    tensor += np.fliplr(fld_extractor(interval))
                else:
                    logger.warning(f"Strand {interval.Strand} not recognized. Treating it as +.")
            else:
                tensor += fld_extractor(interval)

        if config.n_bins_pos:
            tensor = np.hstack(
                list(map(lambda x: x.sum(axis=1)[:, None], np.array_split(tensor, config.n_bins_pos, axis=1))))
        if config.n_bins_len:
            tensor = np.vstack(list(map(lambda x: x.sum(axis=0), np.array_split(tensor, config.n_bins_len, axis=0))))
        if config.subsample:
            n = config.n or int(tensor.sum(axis=0).min())
            tensor = np.apply_along_axis(
                lambda x: subsample_fragment_lengths(x, n),
                0,
                tensor
            )

        return Signal(array=tensor, tags=tag, metadata=None)

    @lbfextract.hookimpl
    def transform_all_intervals(self, single_intervals_transformed_reads: Signal, config: Any,
                                extra_config: Any) -> Signal:
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        array = single_intervals_transformed_reads.array
        col_sums = array.sum(axis=0)
        mask = col_sums > 0
        normalized_array = np.ones_like(array)
        normalized_array[:, mask] = array[:, mask] / col_sums[mask]
        return Signal(array=array, tags=single_intervals_transformed_reads.tags, metadata=None)

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal, extra_config: Any) -> matplotlib.figure.Figure:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """

        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags) if signal.tags else ""
        start_pos = extra_config.ctx["single_signal_transformer_config"].min_fragment_length
        end_pos = signal.array.shape[0] + start_pos
        sample_name = extra_config.ctx["path_to_bam"].stem
        interval_name = extra_config.ctx["path_to_bed"].stem.split(".", 1)[0]
        fig = plot_fragment_length_distribution(signal.array, start_pos, end_pos,
                                                title=f"{sample_name} {signal_type} {interval_name}")
        file_name = f"{time_stamp}__{run_id}__{signal_type}__heatmap.png"
        file_name_sanitized = sanitize_file_name(file_name)

        output_path = extra_config.ctx["output_path"] / file_name_sanitized
        fig.savefig(output_path, dpi=300)

        return fig

    @lbfextract.hookimpl
    def save_signal(self,
                    signal: Signal,
                    extra_config: Any) -> None:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """

        output_path = extra_config.ctx["output_path"]
        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags)
        file_name = f"{time_stamp}__{run_id}__{signal_type}__signal.csv"
        file_name_sanitized = sanitize_file_name(file_name)

        path_to_plot = output_path / file_name_sanitized
        df = pd.DataFrame(signal.array)
        df.to_csv(path_to_plot)
        return path_to_plot


class CliHook:
    r"""
        This CliHook implements the CLI interface for the extract_fragment_length_distribution feature extraction
        method.

        **extract_fragment_length_distribution**

        Given a set of genomic intervals having the same length w, extract_fragment_length_distribution calculates the
        fragment length distribution at each position, which can be represented as:

        .. math::
            \mathbf{d}_l = \left( \frac{1}{|F|} \sum_{\substack{f \in F \\ |f| = p \\  i \in f}} \mathbb{1}  \right)^{p_e}_{p_s}

        Where :math:`l` represents the genomic position, :math:`f` represents a fragment, :math:`p_e` represent
        the maximum fragment length
        and :math:`p_s` represents the minimum fragment length
    """

    @lbfextract.hookimpl_cli
    def get_command(self) -> click.Command:
        @click.command(
            short_help="It extracts the fragment length distribution signal from a BAM file for a BED file provided."
        )
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
                      help="Integer describing the number of bases to be extracted around the middle point of an "
                           "interval present in the bedfile")
        @click.option("--flanking_window", default=1000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted after the window")
        @click.option("--extra_bases", default=2000, type=int, show_default=True,
                      help="Integer describing the number of bases to be extracted from the bamfile when removing the "
                           "unused bases to be sure to get all the proper pairs, which may be mapping up to 2000 bs")
        @click.option("--n_binding_sites", default=1000, type=int, show_default=True,
                      help="number of intervals to be used to extract the signal, if it is higher then the provided"
                           "intervals, all the intervals will be used")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--min_fragment_length", default=100, type=int, show_default=True,
                      help="minimum fragment length to be considered")
        @click.option("--max_fragment_length", default=300, type=int, show_default=True,
                      help="maximum fragment length to be considered")
        @click.option("--n_reads", default=1000, type=int, show_default=True,
                      help="number of reads to be subsampled at each position")
        @click.option("--subsample", is_flag=True, show_default=False,
                      help="Boolean flag. When it is set, the reads are subsampled at each position")
        @click.option("--n_bins_len", type=int, show_default=True,
                      help="number of bins to be used to discretize the fragment length")
        @click.option("--n_bins_pos", type=int, show_default=True,
                      help="number of bins to be used to discretize the position")
        @click.option("--flip_based_on_strand", is_flag=True, show_default=False,
                      default=False,
                      help="Boolean flag. When it is set, the signal is flipped based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        @click.option("--w", default=5, type=int, show_default=True,
                      help="window used for the number of bases around either the middle point in the"
                           " fld_middle_around or the number of bases around the center of the dyad in fld_dyad")
        @click.option("--fld_type",
                      type=click.Choice(["fld", "fld_middle", "fld_middle_n", "fld_dyad", "fld_peter_ulz"],
                                        case_sensitive=False),
                      show_default=True, default="fld",
                      help="type of fragment length distribution to be extracted")
        @click.option("--read_start", default=53, type=int, show_default=True,
                      help="start of the read to be used to extract coverage")
        @click.option("--read_end", default=113, type=int, show_default=True,
                      help="end of the read to be used to extract coverage")
        def extract_fragment_length_distribution(
                path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                output_path: pathlib.Path,
                skip_read_fetching: bool,
                window: int,
                flanking_window: int,
                extra_bases: int,
                n_binding_sites: int,
                min_fragment_length: int,
                max_fragment_length: int,
                n_reads: int,
                subsample: bool,
                n_bins_pos: int,
                n_bins_len: int,
                cores: int,
                w: int,
                fld_type: str,
                exp_id: Optional[str],
                read_start: int,
                read_end: int,
                flip_based_on_strand: bool = False,
                gc_correction_tag: Optional[str] = None,

        ):
            """
            Given a set of genomic intervals having the same length w, the extract_fragment_length_distribution feature
            extraction method extracts the fragment length distribution at each position of the genomic intervals used.
            """
            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = Config({})
            single_signal_transformer_config = {
                "min_fragment_length": min_fragment_length,
                "max_fragment_length": max_fragment_length,
                "signal_transformer": fld_type,
                "n": n_reads,
                "subsample": subsample,
                "n_bins_pos": n_bins_pos,
                "n_bins_len": n_bins_len,
                "flip_based_on_strand": flip_based_on_strand,
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag,
                "read_start": read_start,
                "read_end": read_end

            }
            if w:
                single_signal_transformer_config["w"] = w
            if fld_type == "fld_dyad":
                distribution = calculate_reference_distribution(path_to_sample=path_to_bam,
                                                                min_length=min_fragment_length,
                                                                max_length=max_fragment_length,
                                                                chr="chr12",
                                                                start=34_300_000,
                                                                end=34_500_000
                                                                )
                peaks = get_peaks(distribution) + min_fragment_length
                single_signal_transformer_config["peaks"] = [peaks[0]]
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=["fragment_length_distribution", ],
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=Config({}),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

        return extract_fragment_length_distribution


hook = FextractHooks()
hook_cli = CliHook()
