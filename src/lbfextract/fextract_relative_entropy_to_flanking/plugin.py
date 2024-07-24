import logging
import pathlib
from typing import Any, Optional
from typing import List

import click
import matplotlib
import numpy as np
import scipy.stats
from matplotlib import pyplot as plt

import lbfextract.fextract
from lbfextract.core import App
from lbfextract.fextract.schemas import Config, AppExtraConfig, ReadFetcherConfig
from lbfextract.fextract_fragment_length_distribution.plugin import calculate_reference_distribution, get_peaks
from lbfextract.fextract_fragment_length_distribution.schemas import SingleSignalTransformerConfig
from lbfextract.plotting_lib.plotting_functions import plot_signal
from lbfextract.utils import generate_time_stamp, sanitize_file_name
from lbfextract.utils_classes import Signal

logger = logging.getLogger(__name__)


class FextractHooks:

    @lbfextract.hookimpl
    def transform_all_intervals(self, single_intervals_transformed_reads: Signal, config: Any,
                                extra_config: Any) -> Signal:
        """
        :param single_intervals_transformed_reads: Signal object containing the signals per interval
        :param config: config specific to the function
        :param extra_config: extra configuration that may be used in the hook implementation
        """
        single_intervals_transformed_reads.array += 1e-10
        single_intervals_transformed_reads.array /= single_intervals_transformed_reads.array.sum(axis=0)
        default_flanking_region = single_intervals_transformed_reads.array.shape[1] // 3
        flanking_window = extra_config.ctx["read_fetcher_config"].flanking_region_window or default_flanking_region
        sum_left = single_intervals_transformed_reads.array[:, :flanking_window].sum(axis=1)
        sum_right = single_intervals_transformed_reads.array[:, -flanking_window:].sum(axis=1)
        flanking_array = (sum_left + sum_right) / (flanking_window * 2)
        flanking_distribution = flanking_array / flanking_array.sum()
        relative_entropy_to_flanking = np.apply_along_axis(lambda x: scipy.stats.entropy(x, flanking_distribution),
                                                           0,
                                                           single_intervals_transformed_reads.array)
        return Signal(array=relative_entropy_to_flanking, tags=("relative_entropy_to_flanking",), metadata=None)

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal,
                    config: Any,
                    extra_config: AppExtraConfig) -> matplotlib.figure.Figure:
        signal_type = "_".join(signal.tags) if signal.tags else ""
        with plt.style.context('seaborn-v0_8-whitegrid'):
            fig, ax = plt.subplots(1, figsize=(10, 10))
            ax.set_title(f"{signal_type}\n"
                         f"patient: {extra_config.ctx['path_to_bam'].stem} "
                         f"bed file: {extra_config.ctx['path_to_bed'].stem.split('.', 1)[0]}".capitalize(), fontsize=20)
            fig, _ = plot_signal(signal.array, apply_savgol=False, ax=ax, fig=fig, label=signal_type)
            ax.set_ylabel(signal_type)
            ax.set_xlabel("Position")
        file_name = f"{generate_time_stamp()}__{extra_config.ctx['id']}__{signal_type}_signal_plot.pdf"
        file_name_sanitized = sanitize_file_name(file_name)
        fig.savefig(
            extra_config.ctx["output_path"] / file_name_sanitized,
            dpi=300)
        return fig


class CliHook:
    r"""
        This CliHook implements the CLI interface for the extract_relative_entropy_to_flanking feature extraction
        method.

        **extract_relative_entropy_to_flanking**

        Given a set of genomic intervals having the same length w, extract_relative_entropy_to_flanking calculates the
        Relative Fragment Entropy at each position, which can be represented as:

        .. math::
            RFE_l = D_{KL}(Q_l | F) = \sum_{x \in X} Q_l(x) log(\frac{Q_l(x)}{F(x)})

        Where :math:`l` represents the genomic position, :math:`Q_l` represents the fragment length distribution at
        position :math:`l`, :math:`X` is the alphabet (fragment length range used) and :math:`F` is the fragment length
        distribution in the flanking region
    """

    @lbfextract.hookimpl_cli
    def get_command(self) -> click.Command | List[click.Command]:
        @click.command(short_help="It extracts the relative fragment entropy (RFE) (relative entropy to flanking)"
                                  " signal from a BAM file for a BED file provided.")
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
                      help="window used for the number of bases around either the middle point in "
                           "the fld_middle_around or the number of bases around the center of the dyad in fld_dyad")
        @click.option("--fld_type",
                      type=click.Choice(["fld", "fld_middle", "fld_middle_n", "fld_dyad"],
                                        case_sensitive=False),
                      show_default=True, default="fld",
                      help="type of fragment length distribution to be extracted")
        def extract_relative_entropy_to_flanking(
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
                gc_correction_tag: Optional[str] = None,
                flip_based_on_strand: bool = False,

        ):
            r"""
            Given a set of genomic intervals having the same length w, extract_relative_entropy_to_flanking calculates
            the Relative Fragment Entropy at each position, which can be represented as:

            .. math::
                RFE_l = D_{KL}(Q_l | F) = \sum_{x \in X} Q_l(x) log(\frac{Q_l(x)}{F(x)})

            Where l represents the genomic position, Q_l represents the fragment length distribution at position l and F
            is the fragment length distribution in the flanking region
            """
            read_fetcher_config = {
                "window": window,
                "flanking_region_window": flanking_window,
                "extra_bases": extra_bases,
                "n_binding_sites": n_binding_sites,
            }
            reads_transformer_config = Config({})
            transform_all_intervals_config = {}
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
            res = App(plugins_name=["relative_entropy_to_flanking",
                                    "fragment_length_distribution"],
                      # more plugin names can be added here to inherit hooks from them
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec,
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=Config(transform_all_intervals_config),
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

        return extract_relative_entropy_to_flanking


hook = FextractHooks()
hook_cli = CliHook()
