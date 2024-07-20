import logging
import operator
import pathlib
from functools import reduce
from typing import Any, Optional

import click
import numpy as np
import pandas as pd
import pyranges
import scipy.stats
from matplotlib import pyplot as plt

import lbfextract.fextract.signal_transformer
import lbfextract.fextract_fragment_length_distribution
import lbfextract.fextract_fragment_length_distribution.signal_summarizers
from lbfextract.core import App
from lbfextract.fextract.schemas import Config, ReadFetcherConfig, AppExtraConfig
from lbfextract.fextract_batch_coverage.schemas import PlotConfig
from lbfextract.fextract_entropy_in_batch.schemas import SignalSummarizer
from lbfextract.fextract_fragment_length_distribution.schemas import SingleSignalTransformerConfig
from lbfextract.fextract_fragment_length_distribution_in_batch.plugin import IntervalIteratorFld
from lbfextract.plotting_lib.plotting_functions import plot_signal, plot_heatmap_kde_amplitude, plot_signal_batch
from lbfextract.utils import load_temporary_bed_file, get_tmp_fextract_file_name, filter_bam, generate_time_stamp, \
    check_input_bed, check_input_bam, filter_out_empty_bed_files, sanitize_file_name
from lbfextract.utils_classes import Signal

logger = logging.getLogger(__name__)


def extract_coverage(coverage_extractor, transformed_reads):
    array = np.vstack(transformed_reads.apply(lambda x: coverage_extractor(x), axis=1))
    return array


def identity(x):
    return x


class IntervalIteratorRelativeEntropyFlanking(IntervalIteratorFld):
    def __init__(self, df: pd.DataFrame, path_to_bam: pathlib.Path, multiple_iterators: bool, extra_bases: int,
                 flip_based_on_strand: bool | None = None, flanking_window: int = None):
        super().__init__(df, path_to_bam, multiple_iterators, extra_bases, flip_based_on_strand)
        self.flanking_window = flanking_window

    def __next__(self):
        fld_dict: dict = super().__next__()
        return {k: self.get_relative_entropy(v, self.flanking_window) for k, v in fld_dict.items()}

    @staticmethod
    def get_relative_entropy(array, flanking_window):
        array += 1e-10
        array /= array.sum(axis=0)
        default_flanking_region = array.shape[1] // 3
        flanking_window = flanking_window or default_flanking_region
        sum_left = array[:, :flanking_window].sum(axis=1)
        sum_right = array[:, -flanking_window:].sum(axis=1)
        flanking_array = (sum_left + sum_right) / (flanking_window * 2)
        flanking_distribution = flanking_array / flanking_array.sum()
        relative_entropy_to_flanking = np.apply_along_axis(
            lambda x: scipy.stats.entropy(x, flanking_distribution),
            0,
            array)
        return relative_entropy_to_flanking


class FextractHooks:
    @lbfextract.hookimpl
    def fetch_reads(self,
                    path_to_bam: pathlib.Path,
                    path_to_bed: pathlib.Path,
                    config: Any,
                    extra_config: Any) -> IntervalIteratorRelativeEntropyFlanking:
        """
        :param path_to_bam: path to the bam file
        :param path_to_bed: path to the bed file with the regions to be filtered
        :param config: configuration file containing the configuration object required by the fetch_reads function
        :param extra_config: extra configuration that may be used in the hook implementation
        :return: ReadsPerIntervalContainer object containing all the ReadsPerInterval objects in all the intervals
                 contained in the bed file
        """
        config_f = config.f or 2
        config_F = config.F or 3868

        check_input_bed(path_to_bed)
        check_input_bam(path_to_bam)

        bed_files_paths = filter_out_empty_bed_files(path_to_bed)

        temporary_bed_files_name = [
            load_temporary_bed_file(bed_file=bed,
                                    extra_bases=config.extra_bases,
                                    window=config.window,
                                    flanking_region_window=config.flanking_region_window,
                                    n_binding_sites=config.n_binding_sites,
                                    run_id=extra_config.ctx["run_id"])[1].as_df()
            for bed in bed_files_paths
        ]
        concat_bed = pd.concat(temporary_bed_files_name)
        path_to_tmp_file = get_tmp_fextract_file_name(extra_config.ctx["run_id"])
        pyranges.PyRanges(concat_bed).to_csv(path_to_tmp_file, sep="\t", header=None)

        tmp_bam_file = filter_bam(path_to_bam, path_to_tmp_file, cores=extra_config.cores,
                                  run_id=extra_config.ctx["run_id"],
                                  f=config_f, F=config_F)

        return IntervalIteratorRelativeEntropyFlanking(concat_bed,
                                                       tmp_bam_file,
                                                       multiple_iterators=True,
                                                       extra_bases=config.extra_bases,
                                                       flanking_window=config.flanking_region_window)

    @lbfextract.hookimpl
    def transform_single_intervals(self, transformed_reads: IntervalIteratorRelativeEntropyFlanking,
                                   config: SingleSignalTransformerConfig,
                                   extra_config: Any) -> IntervalIteratorRelativeEntropyFlanking:
        """
        :param transformed_reads: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
            basically lists with information about start and end of the interval
        :param config: config specific to the function
        :param extra_config: config containing context information plus extra parameters
        """
        transformed_reads.flip_based_on_strand = config.flip_based_on_strand
        transformed_reads.max_fragment_length = config.max_fragment_length
        transformed_reads.min_fragment_length = config.min_fragment_length
        transformed_reads.n_bins_pos = config.n_bins_pos
        transformed_reads.n_bins_len = config.n_bins_len
        transformed_reads.subsample = config.subsample
        transformed_reads.n = config.n

        signal_transformers_dict = {
            "fld": {"class": "TfbsFragmentLengthDistribution",
                    "config": {"min_fragment_length": config.min_fragment_length,
                               "max_fragment_length": config.max_fragment_length,
                               "gc_correction": config.gc_correction,
                               "tag": config.tag
                               },
                    "tags": ("fragment_length_distribution",)},

        }
        fld_extractor_params = signal_transformers_dict[config.signal_transformer]["config"]
        fld_extractor = getattr(lbfextract.fextract_fragment_length_distribution.signal_summarizers,
                                signal_transformers_dict[config.signal_transformer]["class"])(
            **fld_extractor_params)
        transformed_reads.set_signal_transformer(fld_extractor)
        return transformed_reads

    @lbfextract.hookimpl
    def transform_all_intervals(self, single_intervals_transformed_reads: IntervalIteratorFld, config: Any,
                                extra_config: Any) -> Signal:
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the transform_all_intervals hook
        :param extra_config: extra configuration that may be used in the hook implementation
        """

        summarized_signal_per_bed = [i for i in single_intervals_transformed_reads]
        summarized_signal_per_bed = pd.DataFrame(reduce(operator.ior, summarized_signal_per_bed, {})).T

        return Signal(
            array=summarized_signal_per_bed,
            metadata=summarized_signal_per_bed.index,
            tags=tuple(["relative_entropy_to_flanking_in_batch", ])
        )

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal, extra_config: Any, config: PlotConfig) -> dict[str, pathlib.Path]:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        fig_pths = {}
        df = pd.DataFrame(signal.array, index=signal.metadata)
        flanking = int((df.shape[1] // 5) * 2) if not config.flanking else config.flanking
        fig_plot_heatmap_kde_amplitude, ax = plot_heatmap_kde_amplitude(
            array=df,
            title=" ".join(signal.tags),
            title_font_size=config.title_font_size if config.title_font_size else 20,
            general_font_size=config.general_font_size if config.general_font_size else 15,
            tf_to_annotate=config.tf_to_annotate if config.tf_to_annotate else None,
            ylabel=config.ylabel if config.ylabel else None,
            flanking=flanking,
            annotation_center_line=config.annotation_center_line if config.annotation_center_line else "interval center",
            window_center=config.window_center if config.window_center else 50,
            top=config.top if config.top else 5,
            bottom=config.bottom if config.bottom else 5,

        )
        signal_type = "_".join(signal.tags)
        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        file_name = f"{time_stamp}__{run_id}__{signal_type}__heatmap_kde_amplitude_plot.png"
        file_name_sanitized = sanitize_file_name(file_name)
        output_path = extra_config.ctx["output_path"] / file_name_sanitized
        fig_plot_heatmap_kde_amplitude.savefig(output_path, dpi=300)
        fig_pths["plot_heatmap_kde_amplitude"] = output_path
        plt.close(fig_plot_heatmap_kde_amplitude)

        fig_plot_signal_batch, ax = plot_signal_batch(
            df,
            apply_savgol=config.apply_savgol if config.apply_savgol else False,
            savgol_window_length=config.savgol_window_length if config.apply_savgol else 11,
            savgol_polyorder=config.savgol_polyorder if config.apply_savgol else 3,
            signal=signal_type,
            title=f"{signal_type} all intervals",
            color=config.color if config.color else "blue",
            label=f"{signal_type} signal summary",
            flanking=flanking,
            xlabel=config.xlabel if config.xlabel else None,
            window_center=config.window_center if config.window_center else 50,
            top=config.top if config.top else 5,
            bottom=config.bottom if config.bottom else 5,
        )
        file_name = f"{time_stamp}__{run_id}__{signal_type}__batch_signals.png"
        file_name_sanitized = sanitize_file_name(file_name)
        output_path = extra_config.ctx["output_path"] / file_name_sanitized
        fig_plot_signal_batch.savefig(output_path, dpi=300)
        fig_pths["plot_signal_batch"] = output_path
        plt.close(fig_plot_signal_batch)

        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags) if signal.tags else ""

        for i in df.index.to_list():
            array = df.loc[i].values
            fig_plot_signal, ax = plot_signal(array, line_type="-")
            sample_name = extra_config.ctx["path_to_bam"].stem
            fig_plot_signal.suptitle(f"{sample_name} {signal_type} {i}")
            file_name = f"{time_stamp}__{run_id}__{signal_type}__{i}__signal.png"
            filename_sanitized = sanitize_file_name(file_name)
            output_path = extra_config.ctx["output_path"] / filename_sanitized
            fig_plot_signal.savefig(output_path, dpi=300)
            fig_pths[i] = output_path
            plt.close(fig_plot_signal)
        return fig_pths

    @lbfextract.hookimpl
    def save_signal(self,
                    signal: Signal,
                    config: Any,
                    extra_config: Any) -> pathlib.Path:
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
        file_path = output_path / file_name_sanitized
        logging.info(f"Saving signal to {file_path}")
        df = pd.DataFrame(signal.array, index=signal.metadata)
        pd.DataFrame(df).to_csv(file_path)

        return file_path


class CliHook:
    @lbfextract.hookimpl_cli
    def get_command(self) -> click.Command:
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
                                                       file_okay=False,
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
                           "unused bases to be sure to get all the proper paires, which may be mapping up to 2000 bs")
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
                      help="window used for the number of baseses around either the middle point in the fld_middle_around "
                           "or the number of bases around the center of the dyad in fld_dyad")
        @click.option("--fld_type",
                      type=click.Choice(["fld", "fld_middle", "fld_middle_n", "fld_dyad"],
                                        case_sensitive=False),
                      show_default=True, default="fld",
                      help="type of fragment length distribution to be extracted")
        def extract_relative_entropy_to_flanking_in_batch(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                                          flip_based_on_strand: bool = False,
                                                          gc_correction_tag: Optional[str] = None
                                                          ):
            """
            Extract the entropy signal for a TF from a bam file and a bed file.
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
                "w": w,
            }
            transform_all_intervals_config = {}
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            App(plugins_name=["relative_entropy_to_flanking_in_batch",
                              "fragment_length_distribution_in_batch",
                              "coverage_in_batch"],
                path_to_bam=path_to_bam,
                path_to_bed=path_to_bed,
                output_path=output_path_interval_spec or pathlib.Path.cwd(),
                skip_read_fetching=skip_read_fetching,
                read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                reads_transformer_config=reads_transformer_config,
                single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                transform_all_intervals_config=SignalSummarizer(transform_all_intervals_config),
                plot_signal_config=PlotConfig(plot_signal_config),
                save_signal_config=Config(save_signal_config),
                extra_config=AppExtraConfig(extra_config),
                id=exp_id).run()

        return extract_relative_entropy_to_flanking_in_batch


hook = FextractHooks()
hook_cli = CliHook()
