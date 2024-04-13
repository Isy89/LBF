import hashlib
import json
import logging
import os
import pathlib
import uuid
from typing import List, Tuple

import matplotlib
import pandas as pd
import pluggy

import lbfextract.hookspecs
from lbfextract.fextract.lib import FextractHooks
from lbfextract.fextract.schemas import AppExtraConfig, Config
from lbfextract.utils import get_tmp_dir, clean_up_temporary_file, list_temporary_files
from lbfextract.utils_classes import Signal, TimerAndMemoryProfiler, Tracer

logger = logging.getLogger(__name__)


def get_plugin_manager(plugins_to_be_registered: List[str]) -> pluggy.PluginManager:
    pm = pluggy.PluginManager("lbfextract")
    pm.add_hookspecs(lbfextract.hookspecs.FextractHooksSpecs)
    # collecting all plugins
    pm.load_setuptools_entrypoints("lbfextract")
    list_valid_plugin_names = [pm.get_name(i) for i in pm.get_plugins()]
    # check if all plugins to be registered are valid
    to_be_registerd = []
    for plugin_name in plugins_to_be_registered:
        if plugin_name not in list_valid_plugin_names:
            logger.warning(f"{plugin_name} is not a valid plugin name")
            logger.warning(f"Valid plugin names are: {list_valid_plugin_names}")
        else:
            to_be_registerd.append(plugin_name)
    # unregister all plugins (load_setuptools_entrypoints register all plugins, but not in the correct order required
    # by the plugin implementation)
    for name in list_valid_plugin_names:
        pm.unregister(name=name)
    # register only the plugins to be registered in the correct order
    pm.register(FextractHooks(), name="base")

    for name in to_be_registerd[::-1]:
        pm.load_setuptools_entrypoints("lbfextract", name=name)

    logger.debug(f"Registered plugins are: {[pm.get_name(i) for i in pm.get_plugins()]}")
    return pm


class App:
    """
    This object is the main object of the LBFextract package. This object implements the signal extraction workflow.
    Each step of the workflow is represented by a hook. Hooks are implemneted by plugins, which are loaded by the
    plugin manager. Given a BAM file, a BED file and all the configuration needded for the different steps, it extracts
    the reads corresponding to the intervals specified in the BED file entries, save them, extracts the signal and save it
    together with possibly generated plots.
    """

    def __init__(self,
                 plugins_name: List[str] = None,
                 path_to_bam: pathlib.Path = None,
                 path_to_bed: pathlib.Path = None,
                 output_path: pathlib.Path = None,
                 skip_read_fetching: bool = False,
                 read_fetcher_config: Config = None,
                 save_fetched_reads_config: Config = None,
                 load_reads_config: Config = None,
                 reads_transformer_config: Config = None,
                 single_signal_transformer_config: Config = None,
                 transform_all_intervals_config: Config = None,
                 plot_signal_config: Config = None,
                 save_signal_config: Config = None,
                 extra_config: AppExtraConfig = None,
                 id: str = None,
                 debug: bool = True):

        self.id = id or str(uuid.uuid1())
        self.plugins_name = plugins_name
        self.path_to_bam = path_to_bam
        self.path_to_bed = path_to_bed
        self.output_path = output_path
        self.skip_read_fetching = skip_read_fetching
        self.read_fetcher_config = read_fetcher_config
        self.save_fetched_reads_config = save_fetched_reads_config
        self.load_reads_config = load_reads_config
        self.single_signal_transformer_config = single_signal_transformer_config
        self.reads_transformer_config = reads_transformer_config
        self.transform_all_intervals_config = transform_all_intervals_config
        self.plot_signal_config = plot_signal_config
        self.save_signal_config = save_signal_config
        self.extra_config = AppExtraConfig({**{"ctx": self.__dict__}, **extra_config.to_dict()})
        self.hook_manager = get_plugin_manager(
            plugins_to_be_registered=self.plugins_name or [])
        self.run_id = self._get_run_id()
        self.debug = debug

        os.environ['run_id'] = self.run_id

        if not self.path_to_bam.exists():
            raise FileNotFoundError(f"path {self.path_to_bam} not existing!")
        if not self.path_to_bed.exists():
            raise FileNotFoundError(f"path {self.path_to_bed} not existing!")
        if not self.output_path.exists():
            logger.info(
                f"output path: {self.output_path} does not exists. creating it with all missing parent folders.")
            self.output_path.mkdir(exist_ok=True, parents=True)

        if self.skip_read_fetching and self.id is None:
            raise ValueError("If you want to skip read fetching, you need to provide an id to load the reads from.")

        if debug:
            self.hook_manager.add_hookcall_monitoring(
                lambda hook_name, hook_impls, kwargs: logger.info(f"hook_name: {hook_name}, hook_impls: {hook_impls}"),
                lambda outcome, hook_name, hook_impls, kwargs: logger.info(f"hook_name: {hook_name},"
                                                                           f"hook_impls: {hook_impls},"
                                                                           f" outcome: {outcome}"),
            )

    def _get_run_id(self):
        h = hashlib.sha1()
        h.update(json.dumps(self.extra_config, sort_keys=True, default=str).encode())
        return h.hexdigest()

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def fetch_reads(self) -> pd.DataFrame:
        return self.hook_manager.hook.fetch_reads(
            path_to_bam=self.path_to_bam,
            path_to_bed=self.path_to_bed,
            config=self.read_fetcher_config,
            extra_config=self.extra_config)

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def load_reads(self) -> pd.DataFrame:
        reads = self.hook_manager.hook.load_fetched_reads(config=self.load_reads_config,
                                                          extra_config=self.extra_config)
        return reads

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def save_fetched_reads(self, reads_per_interval_container: pd.DataFrame):
        self.hook_manager.hook.save_fetched_reads(reads_per_interval_container=reads_per_interval_container,
                                                  config=self.save_fetched_reads_config,
                                                  extra_config=self.extra_config)

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def transform_reads(self, reads_per_interval_container: pd.DataFrame) -> pd.DataFrame:
        transformed_reads = self.hook_manager.hook.transform_reads(
            reads_per_interval_container=reads_per_interval_container,
            config=self.reads_transformer_config,
            extra_config=self.extra_config
        )
        return transformed_reads

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def transform_signal_single_interval(self, transformed_reads: pd.DataFrame) -> Signal:
        signal_single_interval = self.hook_manager.hook.transform_single_intervals(
            transformed_reads=transformed_reads,
            config=self.single_signal_transformer_config,
            extra_config=self.extra_config

        )
        return signal_single_interval

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def transform_all_single_interval_signals(self, single_interval_transformed_signals: Signal) -> Signal:
        all_single_interval_signals_transformed = self.hook_manager.hook.transform_all_intervals(
            single_intervals_transformed_reads=single_interval_transformed_signals,
            config=self.transform_all_intervals_config,
            extra_config=self.extra_config
        )
        return all_single_interval_signals_transformed

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def plot_signal(self, signal: Signal) -> matplotlib.figure.Figure:
        fig = self.hook_manager.hook.plot_signal(signal=signal,
                                                 config=self.plot_signal_config,
                                                 extra_config=self.extra_config)
        return fig

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def save_signal(self, signal: Signal) -> None:
        self.hook_manager.hook.save_signal(
            signal=signal,
            config=self.save_signal_config,
            extra_config=self.extra_config)
        logger.info("Signal was saved")

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def _clean_up_temporary_file(self):
        tmp_dir = get_tmp_dir()
        clean_up_temporary_file(tmp_dir, run_id=self.run_id)

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def _list_temporary_files(self):
        tmp_dir = get_tmp_dir()
        list_temporary_files(tmp_dir, run_id=self.run_id)

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def _log_fextract_exception(self, e: Exception | OSError):
        log_exc = logger.exception if isinstance(e, Exception) else logger.error
        logger.error(f"""An error occurred while extracting the signal. 
                 Please check the log file for more details.
                {f'Temporary files were not deleted. ' if self.debug else ''}
                {f'Temporary files can be found in ' if self.debug else ''}"""
                     )
        log_exc(e)

    @TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                            timer_logger=logger)
    @Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
    def run(self) -> Tuple[Signal, matplotlib.figure.Figure | list[matplotlib.figure.Figure] | None]:
        try:
            if self.skip_read_fetching:
                logger.info("Skipping read fetching. Loading reads from disk.")
                reads_per_interval_container = self.load_reads()
            else:
                logger.info("Fetching reads from bam file.")
                reads_per_interval_container = self.fetch_reads()

            if not self.skip_read_fetching:
                logger.info("Saving fetched reads to disk.")
                self.save_fetched_reads(reads_per_interval_container=reads_per_interval_container)
            logger.info("Transforming Reads.")
            transformed_reads = self.transform_reads(reads_per_interval_container=reads_per_interval_container)
            logger.info("Extracting signal per interval.")
            signal_single_interval = self.transform_signal_single_interval(transformed_reads)
            logger.info("Extracting the final signal from the single interval signals.")
            all_single_interval_signals_transformed = self.transform_all_single_interval_signals(signal_single_interval)
            logger.info("Saving signal.")
            self.save_signal(all_single_interval_signals_transformed)
            logger.info("Plotting signal.")
            fig = self.plot_signal(all_single_interval_signals_transformed)
            return all_single_interval_signals_transformed, fig

        except Exception as e:
            self._log_fextract_exception(e)
            raise e

        finally:
            if self.debug:
                logger.debug(f'TMP FOLDER = {os.environ.get("FRAGMENTOMICS_TMP")}')
                self._list_temporary_files()
                logger.debug(Tracer())
            else:
                self._clean_up_temporary_file()
