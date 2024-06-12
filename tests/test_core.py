import pytest


class TestCore:

    def test_App_fetch_reads(self, app):
        self.reeads = app.fetch_reads()

    @pytest.mark.dependency(depends=['TestCore::test_App_fetch_reads'])
    def test_App_save_fatched_reads(self, app, reads):
        app.save_fetched_reads(reads)

    @pytest.mark.dependency(depends=['TestCore::test_App_save_fatched_reads'])
    def test_App_load_reads(self, app, reads):
        loaded_reads = app.load_reads()
        assert len(reads) == len(loaded_reads)
        assert reads.Chromosome.eq(loaded_reads.Chromosome).all()
        assert reads.Start.eq(loaded_reads.Start).all()
        assert reads.End.eq(loaded_reads.End).all()
        reads.reads_per_interval = reads.reads_per_interval.apply(lambda x: len(list(x)))
        loaded_reads.reads_per_interval = loaded_reads.reads_per_interval.apply(lambda x: len(list(x)))
        assert reads.reads_per_interval.eq(loaded_reads.reads_per_interval).all()

    @pytest.mark.dependency(depends=['TestCore::test_App_fetch_reads'])
    def test_App_extract_signal(self, app, reads):
        transformed_reads = app.transform_reads(reads)
        signal_single_interval = app.transform_signal_single_interval(transformed_reads)
        all_single_interval_signals_transformed = app.transform_all_single_interval_signals(signal_single_interval)

    @pytest.mark.dependency(depends=['TestCore::test_App_extract_signal'])
    def test_App_save_signal(self, app, signal):
        app.save_signal(signal)

    @pytest.mark.dependency(depends=['TestCore::test_App_extract_signal'])
    def test_App_plot_signal(self, app, signal):
        app.plot_signal(signal)
