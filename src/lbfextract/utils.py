import logging
import multiprocessing
import os
import pathlib
import re
import shlex
import subprocess
import sys
import tempfile
from datetime import datetime
from typing import Tuple

import numpy as np
import pandas as pd
import pyranges
import pysam
import yaml
from rich.console import Console
from rich.console import group
from rich.panel import Panel
from rich.text import Text

import lbfextract.fextract
from lbfextract import PROFILER_DEBUG
from lbfextract.utils_classes import TimerAndMemoryProfiler, Tracer

logger = logging.getLogger(__name__)


@TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG, timer_logger=logger)
@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def split_array_in_regions(array, rows, cols):
    columns = np.array_split(array, cols, 1)
    matrix = [np.array_split(ar, rows, 0) for ar in columns]
    return matrix


@TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG, timer_logger=logger)
@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def join_splitted_arrays(splitted_matrix):
    columns = [np.vstack(col) for col in splitted_matrix]
    return np.hstack(columns)


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def adapt_indices(start, end, length_region):
    if end < 0:
        end = 0
    elif end > length_region:
        end = length_region if 0 <= start < length_region else 0
    start = start if 0 <= start < length_region else 0
    if start < end:
        return np.arange(start, end)
    else:
        return None


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def check_input_bed(path_to_bed: pathlib.Path):
    if not path_to_bed.exists():
        raise ValueError(f"The bed file ({path_to_bed}) does not exist")

    if not path_to_bed.is_dir():
        raise ValueError(f"The bed file ({path_to_bed}) should be a directory containing the bed files")


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def check_input_bam(path_to_bam: pathlib.Path):
    if not path_to_bam.exists():
        raise ValueError(f"The bam file ({path_to_bam}) does not exist")

    if path_to_bam.stat().st_size == 0:
        raise ValueError(f"The bam file ({path_to_bam}) is empty")


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def filter_out_empty_bed_files(path_to_bed: pathlib.Path):
    bed_files_paths = path_to_bed.glob("*.bed")
    non_empty_bed_files = []
    for bed_file in bed_files_paths:
        if bed_file.stat().st_size > 0:
            non_empty_bed_files.append(bed_file)
        else:
            # warn the user that the bed file is empty
            logger.warning(f"The bed file ({bed_file}) is empty")

    if len(non_empty_bed_files) == 0:
        raise ValueError(f"The bed directory ({path_to_bed}) does not contain any non-empty bed file")
    return non_empty_bed_files


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def clean_up_temporary_file(tmp, run_id: str):
    logger.info(f"Cleaning up temporary files in {tmp}")
    tmp_files = list(tmp.glob(f"fextract*{run_id}*"))
    for tmp_file in tmp_files:
        logger.debug(f"Deleting {tmp_file}")
        os.remove(tmp_file)


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def signal_handler(sig, frame):
    """
    Function to clean up tmp files when a signal is received.
    Signals may be: signal.SIGINT
    """

    tmp_dir = get_tmp_dir()
    try:
        logger.info('Ctrl+C pressed, deleting tmp files before exiting ... ')
        if not PROFILER_DEBUG:
            try:
                tmp_files = list(tmp_dir.glob(f"fextract*{os.environ.get('run_id')}*"))
                for tmp_file in tmp_files:
                    logger.info(f"Deleting {tmp_file}")
                    os.remove(tmp_file)

            except Exception as e:
                logger.exception(e)
        else:
            try:
                tmp_files = list(tmp_dir.glob(f"fextract*{os.environ.get('run_id')}*"))
                logger.info("The following file were not removed because debug in profiler_config.yml is set to True."
                            "To change this behaviour install the package with debug: False or through pip.")
                for tmp_file in tmp_files:
                    logger.info(f"  - {tmp_file}")

            except Exception as e:
                logger.exception(e)

    except Exception as e:
        logger.exception(e)
        logger.warning(f"Deleting files before exiting failed. There may still be file in the {tmp_dir} "
                       f"folder starting with fextract_{os.environ.get('run_id', 'run_id')} that need to be removed")
    finally:
        logger.info(" Exiting ... ")
        sys.exit(0)


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def list_temporary_files(tmp, run_id: str):
    logger.info(f"listing temporary files in {tmp}")
    tmp_files = list(tmp.glob(f"fextract*{run_id}*"))
    for tmp_file in tmp_files:
        logger.debug(f"Deleting {tmp_file}")


@TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG, timer_logger=logger)
@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def load_reads_from_dir(path_to_dir: pathlib.Path, extra_bases: int) -> pd.DataFrame:
    file_dict = load_yml(path_to_dir / "metadata.yml")
    bed_file = load_bed_file(path_to_dir / file_dict["bed"]).slack(extra_bases)
    samfile = pysam.AlignmentFile(path_to_dir / file_dict["bam"])
    list_of_reads = bed_file.as_df()
    list_of_reads["reads_per_interval"] = list_of_reads.apply(
        lambda x: samfile.fetch(x["Chromosome"], x["Start"], x["End"], multiple_iterators=True),
        axis=1)
    return pyranges.PyRanges(list_of_reads).slack(-extra_bases).as_df()


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def get_tmp_bam_name(path_to_bam: pathlib.Path, run_id: str) -> str:
    return path_to_bam.with_stem(f"fextract_{run_id}_" + path_to_bam.stem
                                 ).with_suffix(".tmp.filtered.bam").name


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def get_tmp_dir() -> pathlib.Path:
    temp_dir = os.environ.get("FRAGMENTOMICS_TMP") or tempfile.gettempdir()
    return pathlib.Path(temp_dir)


@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def get_tmp_fextract_file_name(run_id) -> pathlib.Path:
    temp_dir = os.environ.get("FRAGMENTOMICS_TMP") or tempfile.gettempdir()
    tmp_name = f"fextract_{run_id}" + ".tmp.bed"
    return pathlib.Path(temp_dir) / tmp_name


@TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG, timer_logger=logger)
@Tracer(debug=lbfextract.PROFILER_DEBUG, logger=logger)
def filter_bam(in_bam: pathlib.Path,
               bed: pathlib.Path,
               cores: int = None,
               overwrite=True,
               run_id: str | None = None,
               f: int = 2,
               F: int = 3868):
    cpu_counts = multiprocessing.cpu_count()
    view_cores = max(cores // 3, 1) if cores else max(cpu_counts // 3, 1)
    sort_cores = max(cores - view_cores, 1)
    index_cores = cores or 1
    temp_dir = pathlib.Path(os.environ.get("FRAGMENTOMICS_TMP") or tempfile.gettempdir())
    tmp = get_tmp_bam_name(in_bam, run_id)
    path_to_tmp_file = temp_dir / tmp

    logger.debug(f"path_to_tmp_file: {path_to_tmp_file}")

    samtools_view_cmd = (f"{lbfextract.PATH_TO_SAMTOOLS} view -b -h -F {F} -f {f} --region-file {bed} "
                         f"-@ {view_cores} {in_bam}")
    samtools_sort_command = (f"{lbfextract.PATH_TO_SAMTOOLS} sort -@ {sort_cores} -T {temp_dir} "
                             f"-o {path_to_tmp_file.with_suffix('.sorted.bam')} -")

    logger.debug(f"running \n {samtools_view_cmd + '|' + samtools_sort_command}")
    if not path_to_tmp_file.exists() or overwrite:
        try:
            with subprocess.Popen(shlex.split(samtools_view_cmd),
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE) as view_process:
                with subprocess.Popen(shlex.split(samtools_sort_command),
                                      stdin=view_process.stdout,
                                      stdout=subprocess.PIPE) as sort_process:
                    output, error = sort_process.communicate()

                    if sort_process.returncode != 0:
                        raise subprocess.CalledProcessError(sort_process.returncode, sort_process.args, output=output,
                                                            stderr=error)

        except subprocess.CalledProcessError as e:
            logger.debug(samtools_view_cmd + "|" + samtools_sort_command)
            logger.error(f"Error while filtering bam file {in_bam} with bed file {bed} and sorting it")
            logger.error(e)
            raise e

    cmd_index = f"{lbfextract.PATH_TO_SAMTOOLS} index -@ {index_cores} {path_to_tmp_file.with_suffix('.sorted.bam')}"

    logger.debug(f"running \n {cmd_index}")

    try:
        subprocess.run(shlex.split(cmd_index), check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error while indexing the bam with: {cmd_index}")
        logger.error(e)
        raise e

    return path_to_tmp_file.with_suffix('.sorted.bam')


@TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG, timer_logger=logger)
def load_temporary_bed_file(bed_file: pathlib.Path,
                            extra_bases: int = 2000,
                            window: int = 1000,
                            flanking_region_window: int = 0,
                            n_binding_sites=None,
                            run_id: str = None) -> Tuple[pathlib.Path, pyranges.PyRanges]:
    if bed_file is None or not bed_file.exists():
        raise ValueError("bed_file argument to the load_temporary_bed_file function is None or does not exist")
    if n_binding_sites is None or n_binding_sites <= 0:
        raise ValueError("n_binding_sites must be a positive integer")

    range_window = window + flanking_region_window + extra_bases
    temp_dir = os.environ.get("FRAGMENTOMICS_TMP") or tempfile.gettempdir()
    temporary_bed_file = load_bed_file(bed_file, n_binding_sites)
    middle = (temporary_bed_file.End + temporary_bed_file.Start) // 2
    temporary_bed_file.Start = middle - range_window
    temporary_bed_file.End = middle + range_window
    temporary_bed_file_name = get_tmp_fextract_file_name(run_id)
    path_to_tmp_file = pathlib.Path(temp_dir) / temporary_bed_file_name
    temporary_bed_file.to_csv(path_to_tmp_file, sep="\t", header=None)
    return path_to_tmp_file, temporary_bed_file


@TimerAndMemoryProfiler(debug=lbfextract.PROFILER_DEBUG,
                        timer_logger=logger)
def load_bed_file(bed_file: pathlib.Path, n_sites=None) -> pyranges.PyRanges:
    """Load a bed file and return a pyranges object filtering out
    the chromosomes that are not in {chr1, chr2, ..., chrX, chrY}"""

    if bed_file is None or not bed_file.exists():
        raise ValueError("bed_file argument to the load_bed_file function is None or does not exist")

    chromosomes = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    bed = pyranges.read_bed(str(bed_file))
    bed_chromosome_filtered = bed[bed.Chromosome.isin(chromosomes)]

    if n_sites is None:
        n_sites = bed_chromosome_filtered.length
    elif n_sites > bed_chromosome_filtered.length:
        n_sites = bed_chromosome_filtered.length
    elif n_sites <= 0:
        return pyranges.PyRanges()
    else:
        ValueError(f"n_sites must be a positive integer. This was found: {n_sites=} {type(n_sites)=}")

    try:
        bed_chromosome_filtered = bed_chromosome_filtered.as_df().sort_values(by=["Score"], ascending=False)[:n_sites]
    except KeyError:
        bed_chromosome_filtered = bed_chromosome_filtered.as_df()[:n_sites]

    return pyranges.PyRanges(bed_chromosome_filtered)


def load_yml(path_to_file_pkl):
    with open(path_to_file_pkl, "r") as f:
        loaded_file = yaml.load(f, Loader=yaml.FullLoader)
        return loaded_file


def write_yml(dictionary, file_name):
    with open(file_name, "w") as f:
        logger.info(file_name)
        logger.info(dictionary)
        yaml.dump(dictionary, f)
    return file_name


def generate_time_stamp():
    data_time_obj = datetime.now()
    time_stamp = data_time_obj.strftime("%Y%m%d__%H:%M:%S_%f")
    return time_stamp


def get_logo():
    msg = (
        r'''   
             __       _______  ________                     __                                  __     
            |  \     |       \|        \                   |  \                                |  \    
            | ▓▓     | ▓▓▓▓▓▓▓\ ▓▓▓▓▓▓▓▓ ______  __    __ _| ▓▓_    ______   ______   _______ _| ▓▓_   
            | ▓▓     | ▓▓__/ ▓▓ ▓▓__    /      \|  \  /  \   ▓▓ \  /      \ |      \ /       \   ▓▓ \  
            | ▓▓     | ▓▓    ▓▓ ▓▓  \  |  ▓▓▓▓▓▓\\▓▓\/  ▓▓\▓▓▓▓▓▓ |  ▓▓▓▓▓▓\ \▓▓▓▓▓▓\  ▓▓▓▓▓▓▓\▓▓▓▓▓▓  
            | ▓▓     | ▓▓▓▓▓▓▓\ ▓▓▓▓▓  | ▓▓    ▓▓ >▓▓  ▓▓  | ▓▓ __| ▓▓   \▓▓/      ▓▓ ▓▓       | ▓▓ __ 
            | ▓▓_____| ▓▓__/ ▓▓ ▓▓     | ▓▓▓▓▓▓▓▓/  ▓▓▓▓\  | ▓▓|  \ ▓▓     |  ▓▓▓▓▓▓▓ ▓▓_____  | ▓▓|  \
            | ▓▓     \ ▓▓    ▓▓ ▓▓      \▓▓     \  ▓▓ \▓▓\  \▓▓  ▓▓ ▓▓      \▓▓    ▓▓\▓▓     \  \▓▓  ▓▓
             \▓▓▓▓▓▓▓▓\▓▓▓▓▓▓▓ \▓▓       \▓▓▓▓▓▓▓\▓▓   \▓▓   \▓▓▓▓ \▓▓       \▓▓▓▓▓▓▓ \▓▓▓▓▓▓▓   \▓▓▓▓ 
        
            O       o O         o O       o
            | O   o | | O     o | | O   o |
            | | O | | | | LBF | | | | O | |
            | o   O | | o     O | | o   O |
            o       O o         O o       O
        
        '''
    )

    return msg


def get_starting_msg(logo):
    msg = (
        f'''
        {logo}


        Welcome to the LBFextract package version {lbfextract.__version__}

        For any questions please contact: 
        - <LBFextract@gmail.com>

        If you find any bugs please report them here:
        - <https://github.com/Isy89/LBF/issues>

        Copyright: Original work on LBFextract and accessory code Copyright (c) 2023 Isaac Lazzeri
        Licence: GNU General Public License v3.0


        ''')
    return msg


def start_msg():
    """Function to generate the start  message """
    console = Console()

    title = get_logo()
    with console.capture() as capture:
        console.print(get_starting_msg(title))

    title = capture.get()

    @group()
    def get_panels():
        yield Text.from_ansi(title)

    console.print(Panel(get_panels(), title="LBFextract", expand=True), style="green on black")


def sanitize_file_name(file_name):
    return re.sub(pattern=r'[<>:"/\\|?*\x00-\x1F\s]', repl='-', string=file_name)
