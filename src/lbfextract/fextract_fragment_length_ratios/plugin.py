from __future__ import annotations

import pathlib
from typing import Any, Optional
from typing import List

import click
import matplotlib
from matplotlib import pyplot as plt

import lbfextract.fextract
from lbfextract.core import App
from lbfextract.fextract.schemas import Config, AppExtraConfig, ReadFetcherConfig
from lbfextract.utils import generate_time_stamp, sanitize_file_name
from lbfextract.utils_classes import Signal
from lbfextract.fextract_fragment_length_distribution.schemas import SingleSignalTransformerConfig
from lbfextract.plotting_lib.plotting_functions import plot_signal


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
        nominator_upper_bound = config.nominator_upper_bound - config.min_fragment_length

        nominator_lower_bound = config.nominator_lower_bound - config.min_fragment_length
        nominator = single_intervals_transformed_reads.array[nominator_lower_bound:nominator_upper_bound].sum(axis=0)
        denominator_upper_bound = config.denominator_upper_bound - config.min_fragment_length
        denominator_lower_bound = config.denominator_lower_bound - config.min_fragment_length
        denominator = single_intervals_transformed_reads.array[denominator_lower_bound:denominator_upper_bound].sum(
            axis=0)
        return Signal(array=nominator / denominator, tags=("fl_ratios",), metadata=None)

    @lbfextract.hookimpl
    def plot_signal(self, signal: Signal,
                    config: Any,
                    extra_config: AppExtraConfig) -> matplotlib.figure.Figure:
        signal_type = "_".join(signal.tags) if signal.tags else ""
        with plt.style.context('seaborn-v0_8-whitegrid'):
            fig, ax = plt.subplots(1, figsize=(10, 10))
            ax.set_title(f"{signal_type}\n"
                         f"patient: {extra_config.ctx['path_to_bam'].stem} "
                         f"bed file: {extra_config.ctx['path_to_bed'].stem.split('.', 1)[0]}", fontsize=20)
            fig, _ = plot_signal(signal.array, apply_savgol=False, ax=ax, fig=fig, label=signal_type)
            ax.set_ylabel(signal_type)
            ax.set_xlabel("Position")

        file_name = f"{generate_time_stamp()}__{extra_config.ctx['id']}__{signal_type}_signal_plot.pdf"
        file_name_sanitized = sanitize_file_name(file_name)
        fig.savefig(
            extra_config.ctx["output_path"] / file_name_sanitized,
            dpi=300)
        return fig


def validate_nominator_denominator(ctx, param, value):
    nominator = [int(i) for i in value.split(",")]
    return nominator


class CliHook:
    r"""
        This CliHook implements the CLI interface for the extract_fragment_length_ratios feature extraction method.

        **extract_fragment_length_ratios**

        Given a set of genomic intervals having the same length w, the extract_fragment_length_ratios extracts the
        fragment length ratios of the proportion of reads contained in two different parts of the fragment length
        distribution at each position.

        .. math::
            \mathbf{c} = \left( \frac{\sum_{i \in [n,m)} \mathbf{d}^{l}_{i}}{\sum_{j \in [o, p)} \mathbf{d}^{l}_{j}} \right)^{w}_{l=0}
 
        in which :math:`\mathbf{d}` is the fragment length distribution at position :math:`l`, :math:`n,m` the start and
        the end of the range of fragment lengths used in the nominator and :math:`o,p` are the start and end of the
        range of fragment lengths used in the denominator.
    """

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
        @click.option("--max_fragment_length", default=600, type=int, show_default=True,
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
        @click.option("--w", default=None, type=int, show_default=True,
                      help="window used for the number of bases around either the middle point in the "
                           "fld_middle_around or the number of bases around the center of the dyad in fld_dyad")
        @click.option("--fld_type",
                      type=click.Choice(["fld", "fld_middle", "fld_middle_n", "fld_dyad"],
                                        case_sensitive=False),
                      show_default=True, default="fld",
                      help="type of fragment length distribution to be extracted")
        @click.option("--denominator", default="270,600", type=str, show_default=True,
                      callback=validate_nominator_denominator)
        @click.option("--nominator", default="100,270", type=str, show_default=True,
                      callback=validate_nominator_denominator)
        def extract_fragment_length_ratios(
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
                nominator: tuple,
                denominator: tuple,
                exp_id: Optional[str],
                flip_based_on_strand: bool = False,
                gc_correction_tag: Optional[str] = None,
        ):
            """
            given a set of genomic intervals having the same length w, the extract_fragment_length_ratios extracts the
            fragment length ratios of the reads contained in two different parts of the fragment length distribution at
            each position:

            .. math::
                \frac{|{ f: f \in F \land |f| > x }|}{|{ f: f \in F \land |f| < x }|}

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

            }
            if w:
                single_signal_transformer_config["w"] = w
            plot_signal_config = {}
            save_signal_config = {}
            extra_config = {
                "cores": cores
            }
            if max([*nominator, *denominator]) > max_fragment_length:
                raise ValueError(
                    "the maximum value of the denominator must not be greater than the max_fragment_length ")
            if min([*nominator, *denominator]) < min_fragment_length:
                raise ValueError(
                    "the min value of the denominator must not be lower than the min_fragment_length ")

            transform_all_intervals_config = Config({
                "denominator_upper_bound": nominator[1],
                "denominator_lower_bound": nominator[0],
                "nominator_upper_bound": denominator[1],
                "nominator_lower_bound": denominator[0],
                "min_fragment_length": min_fragment_length,
                "max_fragment_length": max_fragment_length,
            })
            output_path.mkdir(parents=True, exist_ok=True)
            output_path_interval_spec = output_path / f"{path_to_bam.stem}" / f"{path_to_bed.stem}"
            output_path_interval_spec.mkdir(parents=True, exist_ok=True)
            res = App(plugins_name=["fragment_length_ratios",
                                    "fragment_length_distribution"],
                      path_to_bam=path_to_bam,
                      path_to_bed=path_to_bed,
                      output_path=output_path_interval_spec or pathlib.Path.cwd(),
                      skip_read_fetching=skip_read_fetching,
                      read_fetcher_config=ReadFetcherConfig(read_fetcher_config),
                      reads_transformer_config=reads_transformer_config,
                      single_signal_transformer_config=SingleSignalTransformerConfig(single_signal_transformer_config),
                      transform_all_intervals_config=transform_all_intervals_config,
                      plot_signal_config=Config(plot_signal_config),
                      save_signal_config=Config(save_signal_config),
                      extra_config=AppExtraConfig(extra_config),
                      id=exp_id).run()
            return res

        return extract_fragment_length_ratios


hook = FextractHooks()
hook_cli = CliHook()
