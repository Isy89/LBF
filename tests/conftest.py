import pathlib

import pytest
import yaml

from lbfextract.core import App
from lbfextract.fextract.schemas import Config, AppExtraConfig, SignalSummarizer, \
    SingleSignalTransformerConfig, ReadFetcherConfig


def load_yml(path):
    with open(path, "r") as f:
        content = yaml.load(f, Loader=yaml.FullLoader)
    return content


path_to_tests_folder = pathlib.Path(__file__).parent
path_to_bam = path_to_tests_folder  / "test_dataset" / "bam" / "fextract_anonymized_test.bam"
path_to_bed = path_to_tests_folder / "test_dataset" / "multi_bed_ar_ctcf_for_dyads" / "CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed"
output_path = path_to_tests_folder / "test_out" / "test_out"


@pytest.fixture
def app():
    app = App(plugins_name=None,
              path_to_bam=path_to_bam,
              path_to_bed=path_to_bed,
              output_path=output_path,
              read_fetcher_config=ReadFetcherConfig(
                  load_yml(path_to_tests_folder / "test_config" / "read_fetcher_config.yml")),
              reads_transformer_config=Config(
                  load_yml(path_to_tests_folder / "test_config" / "reads_transformer_config.yml")),
              single_signal_transformer_config=SingleSignalTransformerConfig(load_yml(
                  path_to_tests_folder / "test_config" / "single_signal_transformer_config.yml")),
              transform_all_intervals_config=SignalSummarizer(load_yml(
                  path_to_tests_folder / "test_config" / "signal_summarizer_transformer_config.yml")),
              plot_signal_config=Config(
                  load_yml(path_to_tests_folder / "test_config" / "plot_signal_config.yml")),
              save_signal_config=Config(),
              extra_config=AppExtraConfig(load_yml(path_to_tests_folder / "test_config" / "extra_config.yml")))
    return app


@pytest.fixture
def reads(app):
    return app.fetch_reads()


@pytest.fixture
def signal(app, reads):
    transformed_reads = app.transform_reads(reads)
    signal_single_interval = app.transform_signal_single_interval(transformed_reads)
    all_single_interval_signals_transformed = app.transform_all_single_interval_signals(signal_single_interval)
    return all_single_interval_signals_transformed

