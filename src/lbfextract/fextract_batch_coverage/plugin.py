import logging
import operator
import pathlib
from functools import reduce
from typing import Any, Union, List, Optional

import click
import dill as pickle
import matplotlib.figure
import numpy as np
import pandas as pd
import pyranges
import pysam

import lbfextract.fextract
import lbfextract.fextract.signal_transformer
from lbfextract.core import App
from lbfextract.fextract.cli_lib import calculate_reference_distribution, get_peaks
from lbfextract.fextract.schemas import AppExtraConfig, Config, SingleSignalTransformerConfig, ReadFetcherConfig, \
    SignalSummarizer
from lbfextract.fextract_batch_coverage.schemas import PlotConfig
from lbfextract.plotting_lib.plotting_functions import plot_heatmap_kde_amplitude, correlation_map_plot, \
    plot_signal_batch, \
    plot_signal
from lbfextract.utils import load_temporary_bed_file, filter_bam, get_tmp_fextract_file_name, generate_time_stamp, \
    check_input_bed, check_input_bam, filter_out_empty_bed_files
from lbfextract.utils_classes import Signal

logger = logging.getLogger(__name__)


def extract_coverage(coverage_extractor, transformed_reads):
    array = np.vstack(transformed_reads.apply(lambda x: coverage_extractor(x), axis=1))
    return array


def identity(x, axis=None):
    return x


class IntervalIterator:
    def __init__(self, df: pd.DataFrame, path_to_bam: pathlib.Path, multiple_iterators: bool, extra_bases: int,
                 flanking_window: int, flip_based_on_strand: bool | None = None):
        self.df = df
        self.df_by_name = self.df.groupby("Name")
        self.sequence = list(self.df_by_name.groups.keys())
        self.path_to_bam = path_to_bam
        self._index = 0
        self.multiple_iterators = multiple_iterators
        self.extra_bases = extra_bases
        self.signal_transformer = identity
        self.signal_summarizer = identity
        self.flanking_window = flanking_window
        self.flip_based_on_strand = flip_based_on_strand

    def __iter__(self):
        return self

    def set_signal_transformer(self, signal_transformer):
        self.signal_transformer = signal_transformer

    def set_signal_summarizer(self, signal_summarizer):
        self.signal_summarizer = signal_summarizer

    def __next__(self):
        if self._index < len(self.sequence):
            key = self.sequence[self._index]
            bamfile = pysam.AlignmentFile(self.path_to_bam)
            df = self.df_by_name.get_group(key).copy()
            df["reads_per_interval"] = df.apply(
                lambda x: bamfile.fetch(x["Chromosome"], x["Start"], x["End"], multiple_iterators=True),
                axis=1
            )
            df.loc[:, "Start"] = df.Start + self.extra_bases
            df.loc[:, "End"] = df.End - self.extra_bases

            array = np.vstack(
                df.apply(lambda x: self.signal_transformer(x), axis=1)
            )

            if self.flip_based_on_strand:
                array[df["Strand"] == "-", :] = np.fliplr(array[df["Strand"] == "-"])

            self._index += 1

            index_flanking = np.logical_or(
                np.arange(array.shape[1]) < self.flanking_window,
                np.arange(array.shape[1]) > (array.shape[1] - self.flanking_window)
            )

            means_flanking = array[:, index_flanking].mean(axis=1)
            mask = means_flanking != 0
            normalized_array = np.zeros_like(array, dtype=float)
            normalized_array[mask] = array[mask] / means_flanking[mask, None]

            return {key: self.signal_summarizer(normalized_array, axis=0)}
        else:
            raise StopIteration


class FextractHooks:

    @lbfextract.hookimpl
    def fetch_reads(self,
                    path_to_bam: pathlib.Path,
                    path_to_bed: pathlib.Path,
                    config: Any,
                    extra_config: Any) -> IntervalIterator:
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
                                  run_id=extra_config.ctx["run_id"]
                                  , f=config_f, F=config_F)

        return IntervalIterator(concat_bed, tmp_bam_file,
                                multiple_iterators=True,
                                extra_bases=config.extra_bases,
                                flanking_window=config.flanking_region_window)

    @lbfextract.hookimpl
    def save_fetched_reads(self, reads_per_interval_container: IntervalIterator,
                           config: Any,
                           extra_config: Any
                           ) -> None:
        """
        Hook implementing the strategy to save the reads fetched from the intervals
        :param reads_per_interval_container: ReadsPerIntervalContainer containing information about the genomic region
            and the reads mapping to it
        :param extra_config: extra configuration that may be used in the hook implementation
        :return: None
        """

        output_path = extra_config.ctx["output_path"]
        run_id = extra_config.ctx["id"]
        signal_type = extra_config.ctx["single_signal_transformer_config"].signal_transformer
        with open(output_path / f"{run_id}__{signal_type}__reads_iterator.pkl", "wb") as f:
            pickle.dump(reads_per_interval_container, f)
        return output_path / f"{run_id}__{signal_type}__reads_iterator.pkl"

    @lbfextract.hookimpl
    def load_fetched_reads(self, config: Any, extra_config: AppExtraConfig) -> pd.DataFrame:
        """
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        output_path = extra_config.ctx["output_path"]
        run_id = extra_config.ctx["id"]
        signal_type = extra_config.ctx["single_signal_transformer_config"].signal_transformer
        with open(output_path / f"{run_id}__{signal_type}__reads_iterator.pkl", "rb") as f:
            return pickle.load(f)

    @lbfextract.hookimpl
    def transform_single_intervals(self, transformed_reads: IntervalIterator, config: Any,
                                   extra_config: Any) -> IntervalIterator:
        """
        :param transformed_reads: ReadsPerIntervalContainer containing a list of ReadsPerInterval which are
            basically lists with information about start and end of the interval
        :param config: config specific to the function
        :param extra_config: config containing context information plus extra parameters

        """
        transformed_reads.flip_based_on_strand = config.flip_based_on_strand

        signal_transformers_dict = {
            "coverage": {"class": "TFBSCoverage",
                         "config": {"gc_correction": config.gc_correction,
                                    "tag": config.tag}
                         },
            "coverage_dyads": {"class": "TFBSCoverageAroundDyads",
                               "config": {"n": config.n,
                                          "gc_correction": config.gc_correction,
                                          "tag": config.tag,
                                          "peaks": config.peaks}
                               },
            "middle_point_coverage": {"class": "TFBSMiddlePointCoverage",
                                      "config": {"gc_correction": config.gc_correction,
                                                 "tag": config.tag}
                                      },
            "middle_n_points_coverage": {"class": "TFBSNmiddlePointCoverage",
                                         "config": {"n": config.n,
                                                    "gc_correction": config.gc_correction,
                                                    "tag": config.tag}
                                         },
            "sliding_window_coverage": {"class": "TFBSSlidingWindowCoverage",
                                        "config": {"window_size": config.window_size,
                                                   "gc_correction": config.gc_correction,
                                                   "tag": config.tag}
                                        },
            "peter_ulz_coverage": {"class": "PeterUlzCoverage",
                                   "config": {"read_start": config.read_start,
                                              "read_end": config.read_end,
                                              "gc_correction": config.gc_correction,
                                              "tag": config.tag}
                                   },
            "wps_coverage": {"class": "WPSCoverage",
                             "config": {"window_size": config.window_size,
                                        "gc_correction": config.gc_correction,
                                        "min_fragment_length": config.min_fragment_length,
                                        "max_fragment_length": config.max_fragment_length,
                                        "tag": config.tag,
                                        }
                             },
        }
        coverage_extractor_params = signal_transformers_dict[config.signal_transformer]["config"]
        coverage_extractor = getattr(lbfextract.fextract.signal_transformer,
                                     signal_transformers_dict[config.signal_transformer]["class"])(
            **coverage_extractor_params)
        transformed_reads.set_signal_transformer(coverage_extractor)

        return transformed_reads

    @lbfextract.hookimpl
    def transform_all_intervals(self, single_intervals_transformed_reads: IntervalIterator, config: Any,
                                extra_config: Any) -> Signal:
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """

        summary_method = {
            "mean": np.mean,
            "median": np.median,
            "max": np.max,
            "min": np.min,
        }
        single_intervals_transformed_reads.set_signal_summarizer(summary_method[config.summarization_method])

        summarized_signal_per_bed = [i for i in single_intervals_transformed_reads]
        summarized_signal_per_bed = pd.DataFrame(reduce(operator.ior, summarized_signal_per_bed, {})).T

        return Signal(
            array=summarized_signal_per_bed.values,
            metadata=summarized_signal_per_bed.index,
            tags=tuple([extra_config.ctx["single_signal_transformer_config"].signal_transformer,
                        config.summarization_method, ])
        )

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal,
                    config: PlotConfig,
                    extra_config: Any) -> matplotlib.figure.Figure:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        df = pd.DataFrame(signal.array, index=signal.metadata)
        flanking = int((df.shape[1] // 5) * 2) if not config.flanking else config.flanking
        fig, ax = plot_heatmap_kde_amplitude(
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
        output_path = extra_config.ctx[
                          "output_path"] / f"{time_stamp}__{run_id}__{signal_type}__heatmap_kde_amplitude_plot.png"
        fig.savefig(output_path)
        fig = correlation_map_plot(df)
        output_path = extra_config.ctx["output_path"] / f"{time_stamp}__{run_id}__{signal_type}__corr_matrix.png"
        fig.savefig(output_path)
        fig, ax = plot_signal_batch(
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
        output_path = extra_config.ctx["output_path"] / f"{time_stamp}__{run_id}__{signal_type}__batch_signals.png"
        fig.savefig(output_path)
        for i in range(df.shape[0]):
            fig, ax = plot_signal(df.iloc[i, :].values,
                                  label="center",
                                  title=f"{signal_type} interval {df.index[i]}")
            output_path = extra_config.ctx[
                              "output_path"] / f"{time_stamp}__{run_id}__{signal_type}__{i}__signal_{df.index[i].replace('/', '_')}.png"
            fig.savefig(output_path)
        return fig

    @lbfextract.hookimpl
    def save_signal(self,
                    signal: Signal,
                    config: Any,
                    extra_config: Any) -> None:
        """
        :param signal: Signal object containing the signals per interval
        :param extra_config: extra configuration that may be used in the hook implementation
        """

        output_path = extra_config.ctx["output_path"]
        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags)
        file_path = output_path / f"{time_stamp}__{run_id}__{signal_type}__signal.csv"
        logging.info(f"Saving signal to {file_path}")
        df = pd.DataFrame(signal.array, index=signal.metadata)
        pd.DataFrame(df).to_csv(file_path)


class CliHook:
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
                      help='path to the directory containing the bed files to be used')
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
                      type=click.Choice(["mean", "median", "max", "min"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True, show_default=False,
                      default=False,
                      help="Boolean flag. When it is set, the signal is flipped based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        def extract_coverage_in_batch(path_to_bam: pathlib.Path,
                                      path_to_bed: pathlib.Path,
                                      output_path: pathlib.Path,
                                      skip_read_fetching,
                                      window: int,
                                      flanking_window: int,
                                      extra_bases: int,
                                      n_binding_sites: int,
                                      summarization_method: str,
                                      cores: int,
                                      exp_id: Optional[str],
                                      flip_based_on_strand: bool = True,
                                      gc_correction_tag: Optional[str] = None
                                      ):
            """
            extract_coverage_in_batch extracts the fragment coverage from multiple BED files at once and generates 
            a signal for each BED file provided.
            """
            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = {}
            single_signal_transformer_config = {
                "signal_transformer": "coverage",
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
            res = App(plugins_name=["coverage_in_batch", ],
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
                                                       file_okay=False,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the directory containing the bed files to be used')
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
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", type=int, show_default=True, default=1,
                      help="number of cores to be used for the computation")
        @click.option("--flip_based_on_strand", is_flag=True, show_default=False,
                      default=False,
                      help="Boolean flag. When it is set, the signal is flipped based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        @click.option("--n", default=1, type=int, show_default=True,
                      help="number of bases to retain around the dyad")
        def extract_coverage_around_dyads_in_batch(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                                   flip_based_on_strand: bool = True,
                                                   gc_correction_tag: Optional[str] = None
                                                   ):
            """
            extract_coverage_around_dyads_in_batch extracts the fragment coverage around dyads from multiple BED files 
            at once and generates a signal for each BED file provided.
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
            transform_all_intervals_config = {
                "summarization_method": summarization_method,
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

            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=["coverage_in_batch", ],
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
                                                       file_okay=False,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the directory containing the bed files to be used')
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
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--flip_based_on_strand", is_flag=True, show_default=False,
                      default=False,
                      help="Boolean flag. When it is set, the signal is flipped based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        def extract_middle_point_coverage_in_batch(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
                                                   output_path: pathlib.Path,
                                                   skip_read_fetching,
                                                   window: int,
                                                   flanking_window: int,
                                                   extra_bases: int,
                                                   n_binding_sites: int,
                                                   summarization_method: str,
                                                   cores: int,
                                                   exp_id: Optional[str],
                                                   flip_based_on_strand: bool = True,
                                                   gc_correction_tag: Optional[str] = None
                                                   ):
            """
            extract_middle_point_coverage_in_batch extracts the fragment coverage considering only the central position of
            each read from multiple BED files at once and generates a signal for each BED file provided.
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
            res = App(plugins_name=["coverage_in_batch", ],
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
                                                       file_okay=False,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the directory containing the bed files to be used')
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
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
        @click.option("--cores", default=1, type=int, show_default=True,
                      help="number of cores to be used")
        @click.option("--n_middle_pos", default=10, type=int, show_default=True,
                      help="number of position around the middle point to be extracted")
        @click.option("--flip_based_on_strand", is_flag=True, show_default=False,
                      default=False,
                      help="Boolean flag. When it is set, the signal is flipped based on the strand")
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        def extract_middle_n_points_coverage_in_batch(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                                      flip_based_on_strand: bool = True,
                                                      gc_correction_tag: Optional[str] = None):
            """
            extract_middle_n_points_coverage_in_batch extracts the fragment coverage considering only the central positions 
            of each read from multiple BED files at once and generates a signal for each BED file provided.
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
            res = App(plugins_name=["coverage_in_batch", ],
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
                                                       file_okay=False,
                                                       dir_okay=True,
                                                       writable=False,
                                                       readable=True,
                                                       resolve_path=False,
                                                       allow_dash=True,
                                                       path_type=pathlib.Path,
                                                       executable=False),
                      help='path to the directory containing the bed files to be used')
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
        @click.option("--summarization_method", default="mean",
                      type=click.Choice(["mean", "median", "max", "min"]),
                      show_default=True,
                      help=f"method to be used to summarize the signal: {{ mean, median, max, min }}")
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
        def extract_sliding_window_coverage_in_batch(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                                     flip_based_on_strand: bool = True,
                                                     gc_correction_tag: Optional[str] = None):
            """
            extract_sliding_window_coverage_in_batch extracts the fragment coverage using a sliding window aproach to 
            reduce the problems encounterd when using low coverage samples and generates a signal for each BED file provided.
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
            res = App(plugins_name=["coverage_in_batch", ],
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
            return res

        return [extract_coverage_in_batch, extract_coverage_around_dyads_in_batch,
                extract_middle_point_coverage_in_batch,
                extract_middle_point_coverage_in_batch, extract_middle_n_points_coverage_in_batch,
                extract_sliding_window_coverage_in_batch]


hook = FextractHooks()
hook_cli = CliHook()
