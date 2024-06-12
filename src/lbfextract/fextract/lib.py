from __future__ import annotations

import logging
import os
import pathlib
import shutil
import tempfile
from typing import Any

import dill as pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges
import pysam

import lbfextract.fextract
import lbfextract.fextract.signal_transformer
from lbfextract.fextract.schemas import ReadFetcherConfig, Config, SingleSignalTransformerConfig, \
    SignalSummarizer, AppExtraConfig
from lbfextract.utils import filter_bam, load_temporary_bed_file, generate_time_stamp, get_tmp_bam_name, \
    write_yml, load_reads_from_dir
from lbfextract.utils_classes import Signal
from lbfextract.plotting_lib.plotting_functions import plot_signal

logger = logging.getLogger(__name__)


class FextractHooks:

    @lbfextract.hookimpl
    def fetch_reads(self, path_to_bam: pathlib.Path,
                    path_to_bed: pathlib.Path,
                    config: ReadFetcherConfig,
                    extra_config: AppExtraConfig) -> pd.DataFrame:
        # create a temporary file adding the flanking regions and the extra bases required to be sure the reads
        # to be fetched which partially overlapping with the regions will be included
        if not path_to_bam.exists():
            raise ValueError(f"The bam file ({path_to_bam}) does not exist")
        if not path_to_bed.exists():
            raise ValueError(f"The bed file ({path_to_bed}) does not exist")
        if path_to_bed.stat().st_size == 0:
            raise ValueError(f"The bed ({path_to_bed}) file is empty")
        if path_to_bam.stat().st_size == 0:
            raise ValueError(f"The bam file ({path_to_bam}) is empty")

        config_f = config.f or 2
        config_F = config.F or 3868
        temporary_bed_file_name, bed_file = load_temporary_bed_file(bed_file=path_to_bed,
                                                                    extra_bases=config.extra_bases,
                                                                    window=config.window,
                                                                    flanking_region_window=config.flanking_region_window,
                                                                    n_binding_sites=config.n_binding_sites,
                                                                    run_id=extra_config.ctx["run_id"])
        if config.window == 0:
            if any(bed_file.as_df()["Start"] == bed_file.as_df()["End"]):
                raise ValueError("The bed file contains intervals with the same start and end but window is set to 0."
                                 "Please either provide interval of size grater than 0 or set the window size to a value"
                                 "greater than 0")
        # filtering the bam file to avoid having to go through it each time while fatching
        if bed_file.empty:
            raise ValueError("The bed file is empty")

        tmp_bam_file = filter_bam(path_to_bam, temporary_bed_file_name, cores=extra_config.cores,
                                  run_id=extra_config.ctx["run_id"], f=config_f, F=config_F)

        bamfile = pysam.AlignmentFile(tmp_bam_file)
        list_of_reads = bed_file.as_df()
        list_of_reads["reads_per_interval"] = list_of_reads.apply(
            lambda x: bamfile.fetch(x["Chromosome"], x["Start"], x["End"], multiple_iterators=True),
            axis=1
        )
        return pyranges.PyRanges(list_of_reads).slack(-config.extra_bases).as_df()

    @lbfextract.hookimpl
    def save_fetched_reads(self,
                           reads_per_interval_container: pd.DataFrame,
                           config: Config,
                           extra_config: AppExtraConfig
                           ) -> pathlib.Path:
        """
        Hook implementing the strategy to save the reads fetched for the intervals
        :param reads_per_interval_container: ReadsPerIntervalContainer containing information about the genomic region
                                             and the reads mapping to it
        :param extra_config: AppExtraConfig containing the output path
        :return: None
        """

        sample = extra_config.ctx["path_to_bam"].stem
        temp_dir = pathlib.Path(os.environ.get("FRAGMENTOMICS_TMP") or tempfile.gettempdir())
        bam_file_name = get_tmp_bam_name(extra_config.ctx["path_to_bam"],
                                         run_id=extra_config.ctx["run_id"], )
        path_to_tmp_bam = (temp_dir / bam_file_name).with_suffix('.sorted.bam')
        path_to_tmp_bam_index = (temp_dir / bam_file_name).with_suffix('.sorted.bam.bai')
        output_dir = pathlib.Path(extra_config.ctx["output_path"]) / "fatched_reads"
        output_dir.mkdir(parents=True, exist_ok=True)
        (
            reads_per_interval_container[["Chromosome", "Start", "End"]]
            .to_csv(
                output_dir / f"{sample}_fetched_reads.bed",
                sep="\t",
                header=False,
                index=False
            )
        )
        write_yml(
            {
                "bed": f"{sample}_fetched_reads.bed",
                "bam": path_to_tmp_bam.name
            },
            output_dir / "metadata.yml"
        )
        if not (output_dir / path_to_tmp_bam.name).exists():
            logger.debug(f"Moving {path_to_tmp_bam} to {output_dir}")
            shutil.move(str(path_to_tmp_bam), str(output_dir))
        if not (output_dir / path_to_tmp_bam_index.name).exists():
            logger.debug(f"Moving {path_to_tmp_bam_index} to {output_dir}")
            shutil.move(str(path_to_tmp_bam_index), str(output_dir))
        return output_dir

    @lbfextract.hookimpl
    def load_fetched_reads(self, config: Config, extra_config: AppExtraConfig) -> pd.DataFrame:

        return load_reads_from_dir(
            pathlib.Path(extra_config.ctx["output_path"]) / "fatched_reads",
            extra_bases=extra_config.ctx["read_fetcher_config"].extra_bases
        )

    @lbfextract.hookimpl
    def transform_reads(self,
                        reads_per_interval_container: pd.DataFrame,
                        config: Config,
                        extra_config: AppExtraConfig) -> pd.DataFrame:
        return reads_per_interval_container

    @lbfextract.hookimpl
    def transform_single_intervals(self,
                                   transformed_reads: pd.DataFrame,
                                   config: SingleSignalTransformerConfig,
                                   extra_config: AppExtraConfig
                                   ) -> Signal:
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
        array = np.vstack(transformed_reads.apply(lambda x: coverage_extractor(x), axis=1))

        if config.flip_based_on_strand:
            array[transformed_reads["Strand"] == "-", :] = np.fliplr(array[transformed_reads["Strand"] == "-"])

        return Signal(
            array=array,
            metadata={"bed_file_df": transformed_reads[['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand']]},
            tags=(config.signal_transformer,)
        )

    @lbfextract.hookimpl
    def transform_all_intervals(self,
                                single_intervals_transformed_reads: Signal,
                                config: SignalSummarizer,
                                extra_config: AppExtraConfig) -> Signal:
        config.bed_file = extra_config.ctx["path_to_bed"]
        default_flanking_region = single_intervals_transformed_reads.array.shape[1] // 3
        flanking_window = extra_config.ctx["read_fetcher_config"].flanking_region_window or default_flanking_region
        summary_method = {
            "mean": np.mean,
            "median": np.median,
            "max": np.max,
            "min": np.min,
            "skip": lambda x, axis: x
        }
        index_flanking = np.logical_or(
            np.arange(single_intervals_transformed_reads.array.shape[1]) < flanking_window,
            np.arange(single_intervals_transformed_reads.array.shape[1]) > (
                    single_intervals_transformed_reads.array.shape[1] - flanking_window
            )
        )

        normalized_array = np.zeros_like(single_intervals_transformed_reads.array)
        means_flanking = single_intervals_transformed_reads.array[:, index_flanking].mean(axis=1)
        mask = means_flanking != 0
        normalized_array[mask, :] = single_intervals_transformed_reads.array[mask, :] / means_flanking[mask, None]
        array = summary_method[config.summarization_method](normalized_array, axis=0)

        return Signal(
            array=array,
            metadata=None,
            tags=tuple(list(single_intervals_transformed_reads.tags) + [config.summarization_method]))

    @lbfextract.hookimpl
    def save_signal(self,
                    signal: Signal,
                    config: Any,
                    extra_config: AppExtraConfig
                    ) -> pathlib.Path:
        output_path = extra_config.ctx["output_path"]
        time_stamp = generate_time_stamp()
        run_id = extra_config.ctx["id"]
        signal_type = "_".join(signal.tags)
        with open(output_path / f"{time_stamp}__{run_id}__{signal_type}__signal.pkl", "wb") as f:
            pickle.dump(signal, f)
        return output_path / f"{time_stamp}__{run_id}__{signal_type}__signal.pkl"

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
        fig.savefig(
            extra_config.ctx["output_path"] /
            f"{generate_time_stamp()}__{extra_config.ctx['id']}__{signal_type}_signal_plot.pdf",
            dpi=300)
        return fig
