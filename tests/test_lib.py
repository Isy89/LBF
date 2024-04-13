import pathlib

import pandas as pd
import pytest

from lbfextract.fextract.lib import FextractHooks
from lbfextract.fextract.schemas import Config, SingleSignalTransformerConfig, AppExtraConfig, ReadFetcherConfig
from lbfextract.utils import load_reads_from_dir
from lbfextract.utils_classes import Signal
from .conftest import load_yml

path_to_tests_folder = pathlib.Path(__file__).parent
path_to_bam = path_to_tests_folder / "test_dataset" / "bam" / "fextract_anonymized_test.bam"
path_to_bed = path_to_tests_folder / "test_dataset" / "multi_bed_ar_ctcf_for_dyads" / "CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed"
output_path = path_to_tests_folder / "test_out"


class TestFextractHooks:

    @pytest.mark.parametrize(
        "path_to_bam, path_to_bed, config, extra_config",
        [
            # testing error with different number of binding sites
            (path_to_bam,
             path_to_bed,
             {"window": 1000,
              "flanking_region_window": 1000,
              "extra_bases": 1000,
              "n_binding_sites": 1000},
             {"cores": 1}
             ),
            (path_to_bam,
             path_to_bed,
             {"window": 1000,
              "flanking_region_window": 1000,
              "extra_bases": 1000,
              "n_binding_sites": 10},
             {"cores": 1}),
            (path_to_bam,
             path_to_bed,
             {"window": 1000,
              "flanking_region_window": 1000,
              "extra_bases": 1000,
              "n_binding_sites": 1},
             {"cores": 1}),
            (path_to_bam,
             path_to_bed,
             {"window": 1000,
              "flanking_region_window": 1000,
              "extra_bases": 1000,
              "n_binding_sites": 0},
             {"cores": 1}),
            (path_to_bam,
             path_to_bed,
             {"window": 0,
              "flanking_region_window": 0,
              "extra_bases": 0,
              "n_binding_sites": 1000},
             {"cores": 1})
        ]
    )
    def test_fetch_reads(self, path_to_bam, path_to_bed, config, extra_config):
        bed_input = pd.read_csv(path_to_bed, sep="\t", header=None)
        config, extra_config = ReadFetcherConfig(config), AppExtraConfig(extra_config)
        extra_config.ctx = {
            "output_path": output_path,
            "path_to_bam": path_to_bam,
            "path_to_bed": path_to_bed,
            "run_id": "TEST_RUN_ID",
            "cores": 1
        }
        if bed_input.shape[0] <= 0 and config.n_binding_sites <= 0:
            with pytest.raises(ValueError):
                FextractHooks().fetch_reads(path_to_bam=path_to_bam,
                                            path_to_bed=path_to_bed,
                                            config=config,
                                            extra_config=extra_config)

        if bed_input.shape[0] <= 0 and config.n_binding_sites > 0:
            with pytest.raises(ValueError):
                FextractHooks().fetch_reads(path_to_bam=path_to_bam,
                                            path_to_bed=path_to_bed,
                                            config=config,
                                            extra_config=extra_config)
        if bed_input.shape[0] > 0 and config.n_binding_sites <= 0:
            with pytest.raises(ValueError):
                FextractHooks().fetch_reads(path_to_bam=path_to_bam,
                                            path_to_bed=path_to_bed,
                                            config=config,
                                            extra_config=extra_config)

        sum_window = config.window + config.flanking_region_window + config.extra_bases

        if bed_input.shape[0] > 0 and config.n_binding_sites > 0 and sum_window > 0:
            reads = FextractHooks().fetch_reads(path_to_bam=path_to_bam,
                                                path_to_bed=path_to_bed,
                                                config=config,
                                                extra_config=extra_config)

        if bed_input.shape[0] > 0 and config.n_binding_sites > 0 and sum_window == 0:
            with pytest.raises(ValueError) as exc_info:
                FextractHooks().fetch_reads(path_to_bam=path_to_bam,
                                                    path_to_bed=path_to_bed,
                                                    config=config,
                                                    extra_config=extra_config)
                assert exc_info.value == (
                    "The bed file contains intervals with the same start and end but window is set to 0."
                    "Please either provide interval of size grater than 0 or set the window size to a value"
                    "greater than 0")

    def test_save_fatched_reads(self, app, reads):
        path_to_saved_reads = FextractHooks().save_fetched_reads(
            reads_per_interval_container=reads,
            config=Config(),
            extra_config=AppExtraConfig({"ctx": app.extra_config.ctx })
        )
        loaded_reads = load_reads_from_dir(path_to_saved_reads, app.read_fetcher_config.extra_bases)
        assert len(reads) == len(loaded_reads)
        assert reads.Chromosome.eq(loaded_reads.Chromosome).all()
        assert reads.Start.eq(loaded_reads.Start).all()
        assert reads.End.eq(loaded_reads.End).all()
        reads.reads_per_interval = reads.reads_per_interval.apply(lambda x: len(list(x)))
        loaded_reads.reads_per_interval = loaded_reads.reads_per_interval.apply(lambda x: len(list(x)))
        assert reads.reads_per_interval.eq(loaded_reads.reads_per_interval).all()

    def test_transform_reads(self, app, reads):
        config = Config(load_yml(path_to_tests_folder / "test_config" / "reads_transformer_config.yml"))
        reads_per_interval_container = FextractHooks().transform_reads(reads_per_interval_container=reads,
                                                                       config=config,
                                                                       extra_config=app.extra_config)
        assert reads.equals(reads_per_interval_container)

    @pytest.mark.parametrize(
        "config, extra_config",
        [
            (path_to_tests_folder / "test_config" / "single_signal_transformer_config.yml", {"cores": 1}),
            (path_to_tests_folder / "test_config" / "single_signal_transformer_config.yml", {"cores": 10}),
        ]
    )
    def test_transform_single_intervals(self, app, reads, config: pathlib.Path, extra_config: dict):
        assert config
        FextractHooks().transform_single_intervals(
            transformed_reads=reads,
            config=SingleSignalTransformerConfig(load_yml(config)),
            extra_config=AppExtraConfig({**extra_config, **{"ctx": app.__dict__}})
        )

    @pytest.mark.parametrize(
        "config, extra_config",
        [
            (path_to_tests_folder / "test_config" / "signal_summarizer_transformer_config.yml", {"cores": 1}),
            (path_to_tests_folder / "test_config" / "signal_summarizer_transformer_config.yml", {"cores": 10})
        ]
    )
    def test_transform_all_intervals(self, app, reads, signal: Signal, config: pathlib.Path, extra_config: dict):
        app.extra_config = AppExtraConfig({**extra_config, **{"ctx": app.__dict__}})
        app.extract_signal(reads)

    def test_plot_signal(self, app):
        pass

    def test_save_signal(self, app):
        pass


if __name__ == "__main__":
    print(f"path: {path_to_tests_folder=}")
    print(load_yml(path_to_tests_folder / "test_config" / "read_fetcher_config.yml"))
