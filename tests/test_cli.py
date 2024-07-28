import logging
import pathlib

from click.testing import CliRunner

from lbfextract.cli import cli

path_to_tests_folder = pathlib.Path(__file__).parent
path_to_bam = path_to_tests_folder / "test_dataset" / "bam" / "fextract_anonymized_test.bam"
path_to_bed = path_to_tests_folder / "test_dataset" / "multi_bed_ar_ctcf_for_dyads" / "CTCF.sorted.gtrd_version_21_12.1000_sites.hg38.bed"
path_to_bed_dir = path_to_tests_folder / "test_dataset" / "multi_bed_ar_ctcf_for_dyads"
path_to_bam_for_in_batch_tests = path_to_tests_folder / "test_dataset" / "bam" / "fextract_anonymized_test.bam"
path_to_bam_for_dyads = path_to_tests_folder / "test_dataset" / "bam" / "fextract_anonymized_test.bam"
output_path = path_to_tests_folder / "test_out" / "test_cli"

logger = logging.getLogger(__name__)


class TestCli:

    def test_extract_coverage(self):
        cmd = [
            'feature_extraction_commands',
            'extract-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "mean",
            "--cores", "4",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_coverage_skip(self):
        cmd = [
            'feature_extraction_commands',
            'extract-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "skip",
            "--cores", "4",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage_skip")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_coverage_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            'extract-coverage-in-batch',
            "--path_to_bam", str(path_to_bam_for_in_batch_tests),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
            "--summarization_method", "mean",
            "--cores", "4",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)

        assert result.exit_code == 0

    def test_extract_coverage_around_dyads(self):
        cmd = [
            'feature_extraction_commands',
            'extract-coverage-dyads',
            "--path_to_bam", str(path_to_bam_for_dyads),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "mean",
            "--cores", "4",
            "--n", "30",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage_around_dyads")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_middle_point_coverage(self):
        cmd = [
            'feature_extraction_commands',
            'extract-middle-point-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "mean",
            "--cores", "4",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_middle_point_coverage")
        logger.info(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        print(" ".join(cmd))
        assert result.exit_code == 0

    def test_middle_n_points_coverage(self):
        cmd = [
            'feature_extraction_commands',
            'extract-middle-n-points-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "mean",
            "--cores", "4",
            "--n_middle_pos", "40",
        ]
        runner = CliRunner()
        logger.info("Testing test_middle_n_points_coverage")
        logger.info(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_sliding_window_coverage(self):
        cmd = [
            'feature_extraction_commands',
            'extract-sliding-window-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "mean",
            "--cores", "4",
            "--window_size", "5"
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_sliding_window_coverage")
        logger.info(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld")
        logger.info(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld_middle(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
            "--fld_type", "fld_middle"
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld_middle")
        logger.info(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld_middle_n(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
            "--max_fragment_length",
            "400",
            "--fld_type", "fld_middle_n"
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld_middle_n")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld_dyad(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution',
            "--path_to_bam", str(path_to_bam_for_dyads),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
            "--max_fragment_length",
            "400",
            "--fld_type", "fld_dyad"
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld_dyad")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld_peter_ulz(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
            "--fld_type", "fld_peter_ulz"
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld_peter_ulz")
        logger.info(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_entropy(self):
        cmd = ['feature_extraction_commands',
               'extract-entropy',
               "--path_to_bam", str(path_to_bam),
               "--path_to_bed", str(path_to_bed),
               "--output_path", str(output_path),
               ]
        runner = CliRunner()
        logger.info("Testing sync")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_entropy_in_batch(self):
        cmd = ['feature_extraction_commands',
               'extract-entropy-in-batch',
               "--path_to_bam", str(path_to_bam_for_in_batch_tests),
               "--path_to_bed", str(path_to_bed_dir),
               "--output_path", str(output_path),
               "--n_reads", str(25000),
               "--subsample"
               ]
        runner = CliRunner()
        logger.info("Testing test_extract_entropy_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_relative_entropy(self):
        cmd = [
            'feature_extraction_commands',
            'extract-relative-entropy-to-flanking',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_relative_entropy")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_relative_entropy_fld_middle_n(self):
        cmd = [
            'feature_extraction_commands',
            'extract-relative-entropy-to-flanking',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
            "--fld_type", "fld_middle_n"

        ]
        runner = CliRunner()
        logger.info("Testing test_extract_relative_entropy_fld_middle_n")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_relative_entropy_fld_middle(self):
        cmd = [
            'feature_extraction_commands',
            'extract-relative-entropy-to-flanking',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_relative_entropy_fld_middle")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_coverage_central60_bases(self):
        cmd = [
            'feature_extraction_commands',
            'extract-peter-ulz-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage_central60_bases")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_wps_coverage(self):
        cmd = [
            'feature_extraction_commands',
            'extract-wps-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_wps_coverage")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_ratios(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-ratios',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_ratios")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))

        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_middle_n_points_coverage_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            'extract-middle-n-points-coverage-in-batch',
            "--path_to_bam", str(path_to_bam_for_in_batch_tests),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_middle_n_points_coverage_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))

        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_middle_point_coverage_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            'extract-middle-point-coverage-in-batch',
            "--path_to_bam", str(path_to_bam_for_in_batch_tests),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_middle_point_coverage_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_sliding_window_coverage_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            ' extract-sliding-window-coverage-in-batch',
            "--path_to_bam", str(path_to_bam_for_in_batch_tests),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_sliding_window_coverage_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_coverage_around_dyads_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            ' extract-coverage-around-dyads-in-batch',
            "--path_to_bam", str(path_to_bam_for_in_batch_tests),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage_around_dyads_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution-in-batch',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))

        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_fragment_length_distribution_fld_dyads_in_batch(self):
        cmd = [
            'feature_extraction_commands',
            'extract-fragment-length-distribution-in-batch',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed_dir),
            "--output_path", str(output_path),
            "--subsample",
            "--fld_type", "fld_dyad",
            "--n_reads", "20000",
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_fragment_length_distribution_fld_dyads_in_batch")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))

        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_coverage__gc_corrected_skip(self):
        cmd = [
            'feature_extraction_commands',
            'extract-coverage',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--summarization_method", "skip",
            "--cores", "4",
            "--gc_correction_tag", "GC"
        ]
        runner = CliRunner()
        logger.info("Testing test_extract_coverage__gc_corrected_skip")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0

    def test_extract_relative_entropy_to_flanking_gc_corrected(self):
        cmd = [
            'feature_extraction_commands',
            'extract-relative-entropy-to-flanking',
            "--path_to_bam", str(path_to_bam),
            "--path_to_bed", str(path_to_bed),
            "--output_path", str(output_path),
            "--subsample",
            "--n_reads", "20000",
            "--gc_correction_tag", "GC"

        ]
        runner = CliRunner()
        logger.info("Testing test_extract_relative_entropy_to_flanking_gc_corrected")
        logger.info(" ".join(cmd))
        print(" ".join(cmd))
        result = runner.invoke(cli, cmd)
        assert result.exit_code == 0
