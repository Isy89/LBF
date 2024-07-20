import logging
import operator
import pathlib
from functools import reduce
from typing import Any, Optional

import click
import numpy as np
import pandas as pd
import pyranges
import pysam
from matplotlib import pyplot as plt

import lbfextract.fextract.signal_transformer
import lbfextract.fextract_fragment_length_distribution
import lbfextract.fextract_fragment_length_distribution.signal_summarizers
from lbfextract.core import App
from lbfextract.fextract.schemas import Config, ReadFetcherConfig, AppExtraConfig
from lbfextract.fextract_batch_coverage.schemas import PlotConfig
from lbfextract.fextract_entropy_in_batch.schemas import SignalSummarizer
from lbfextract.fextract_fragment_length_distribution.plugin import calculate_reference_distribution, get_peaks
from lbfextract.fextract_fragment_length_distribution.schemas import SingleSignalTransformerConfig
from lbfextract.plotting_lib.plotting_functions import plot_fragment_length_distribution
from lbfextract.utils import load_temporary_bed_file, get_tmp_fextract_file_name, filter_bam, generate_time_stamp, \
    check_input_bed, check_input_bam, filter_out_empty_bed_files, sanitize_file_name
from lbfextract.utils_classes import Signal

logger = logging.getLogger(__name__)


def subsample_fragment_lengths(x, n):
    mask = x > 0
    probabilities = np.zeros_like(x)
    probabilities[mask] = x[mask] / x[mask].sum()
    subsampled_reads = np.random.choice(len(x), size=n, p=probabilities, replace=True)
    new_x = np.bincount(subsampled_reads, minlength=len(x))
    return new_x


def optimize_tensor_subsampling(tensor, n):
    num_rows, num_columns = tensor.shape
    subsampled_tensor = np.zeros_like(tensor)
    for col in range(num_columns):
        subsampled_tensor[:, col] = subsample_fragment_lengths(tensor[:, col], n)
    return subsampled_tensor


def identity(x):
    return x


class IntervalIteratorFld:
    def __init__(self, df: pd.DataFrame, path_to_bam: pathlib.Path, multiple_iterators: bool, extra_bases: int,
                 flip_based_on_strand: bool | None = None):
        self.df = df
        self.df_by_name = self.df.groupby("Name")
        self.sequence = list(self.df_by_name.groups.keys())
        self.path_to_bam = path_to_bam
        self._index = 0
        self.multiple_iterators = multiple_iterators
        self.extra_bases = extra_bases
        self.signal_transformer = identity
        self.max_fragment_length = None
        self.min_fragment_length = None
        self.n_bins_pos = None
        self.n_bins_len = None
        self.subsample = None
        self.n = None
        self.flip_based_on_strand = flip_based_on_strand
        self.bamfile = pysam.AlignmentFile(self.path_to_bam)

    def __iter__(self):
        return self

    def set_signal_transformer(self, signal_transformer):
        self.signal_transformer = signal_transformer

    def __getstate__(self):
        state = self.__dict__.copy()
        state["bamfile"] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.path_to_bam = pysam.AlignmentFile(state['path_to_bam'])

    def __next__(self):
        if self._index < len(self.sequence):
            key = self.sequence[self._index]
            df = self.df_by_name.get_group(key).copy()
            df["reads_per_interval"] = [
                self.bamfile.fetch(row.Chromosome, row.Start, row.End, multiple_iterators=False)
                for row in df.itertuples()
            ]

            df["Start"] += self.extra_bases
            df["End"] -= self.extra_bases

            relative_fragment_len = self.max_fragment_length - self.min_fragment_length
            region_length = df["End"].iloc[0] - df["Start"].iloc[0]

            tensor = np.zeros((relative_fragment_len, region_length))
            for interval in df.itertuples():
                if self.flip_based_on_strand:
                    if interval.Strand == "+":
                        tensor += self.signal_transformer(interval)
                    if interval.Strand == "-":
                        tensor += np.fliplr(self.signal_transformer(interval))
                    else:
                        logger.warning(f"Strand {interval.Strand} not recognized. Trating as +")
                        tensor += self.signal_transformer(interval)
                else:
                    tensor += self.signal_transformer(interval)
            del df

            def bin_tensor(t, bins, axis):
                bin_edges = np.linspace(0, t.shape[axis], bins + 1, dtype=int)
                return np.add.reduceat(t, bin_edges[:-1], axis=axis)

            if self.n_bins_pos:
                tensor = bin_tensor(tensor, self.n_bins_pos, axis=1)

            if self.n_bins_len:
                tensor = bin_tensor(tensor, self.n_bins_len, axis=0)

            if self.subsample:
                n = self.n or int(tensor.sum(axis=0).min())
                tensor = optimize_tensor_subsampling(tensor, n)

            col_sums = tensor.sum(axis=0)
            non_zero_mask = col_sums > 0
            tensor[:, non_zero_mask] /= col_sums[non_zero_mask]

            self._index += 1
            return {key: tensor}
        else:
            raise StopIteration


class FextractHooks:
    @lbfextract.hookimpl
    def fetch_reads(self,
                    path_to_bam: pathlib.Path,
                    path_to_bed: pathlib.Path,
                    config: Any,
                    extra_config: Any) -> IntervalIteratorFld:
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

        return IntervalIteratorFld(concat_bed,
                                   tmp_bam_file,
                                   multiple_iterators=True,
                                   extra_bases=config.extra_bases)

    @lbfextract.hookimpl
    def transform_single_intervals(self, transformed_reads: IntervalIteratorFld, config: SingleSignalTransformerConfig,
                                   extra_config: Any) -> IntervalIteratorFld:
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
        summarized_signal_per_bed = reduce(operator.ior, summarized_signal_per_bed, {})

        return Signal(
            array=summarized_signal_per_bed,
            metadata=None,
            tags=tuple(["fragment_length_distribution_in_batch", ])
        )

    @lbfextract.hookimpl
    def save_signal(self,
                    signal: Signal,
                    config: Any,
                    extra_config: Any) -> pathlib.Path:
        """
        :param signal: Signal object containing the signals per interval
        :param config: config specific to the save signal hook
        :param extra_config: extra configuration that may be used in the hook implementation
        """

        output_path = extra_config.ctx["output_path"]
        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags)
        file_path = output_path / f"{time_stamp}__{run_id}__{signal_type}__signal"
        logger.info(f"Saving signal to {file_path}")
        np.savez_compressed(file_path, **signal.array)
        return file_path

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal, extra_config: Any, config: PlotConfig, ) -> dict[str, pathlib.Path]:
        """
        :param signal: Signal object containing the signals per interval
        :param config: object containing the configuration to create the plots
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags) if signal.tags else ""

        fig_pths = {}
        for i in signal.array:
            file_name = f"{time_stamp}__{run_id}__{signal_type}__{i}__heatmap.png"
            file_name_sanitized = sanitize_file_name(file_name)
            array = signal.array[i]
            start_pos = extra_config.ctx["single_signal_transformer_config"].min_fragment_length
            end_pos = array.shape[0] + start_pos
            fig = plot_fragment_length_distribution(array, start_pos, end_pos)
            sample_name = extra_config.ctx["path_to_bam"].stem
            fig.suptitle(f"{sample_name} {signal_type} {i}".capitalize(), fontsize=25)
            output_path = extra_config.ctx["output_path"] / file_name_sanitized
            fig.savefig(output_path, dpi=300)
            fig_pths[i] = output_path
            plt.close(fig)

        return fig_pths


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
        @click.option('--gc_correction_tag', type=str,
                      default=None, help='tag to be used to extract gc coefficient per read from a bam file')
        @click.option("--fld_type",
                      type=click.Choice(["fld", "fld_middle", "fld_middle_n", "fld_dyad"],
                                        case_sensitive=False),
                      show_default=True, default="fld",
                      help="type of fragment length distribution to be extracted")
        @click.option("--w", default=5, type=int, show_default=True,
                      help="window used for the number of bases around either the middle point in the "
                           "fld_middle_around or the number of bases around the center of the dyad in fld_dyad")
        def extract_fragment_length_distribution_in_batch(path_to_bam: pathlib.Path, path_to_bed: pathlib.Path,
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
                                                          fld_type: str,
                                                          w: int,
                                                          exp_id: Optional[str],
                                                          gc_correction_tag: Optional[str] = None):
            """
            Given a set of genomic intervals having the same length w, the extract_fragment_length_distribution feature
            extraction method extracts the fragment length distribution at each position of the genomic intervals used
            for multiple BED files at the same time.
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
                "gc_correction": True if gc_correction_tag else False,
                "tag": gc_correction_tag
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
            transform_all_intervals_config = {}
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            App(plugins_name=["fragment_length_distribution_in_batch", "coverage_in_batch"],
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

        return extract_fragment_length_distribution_in_batch


hook = FextractHooks()
hook_cli = CliHook()
