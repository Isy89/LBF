import logging
import os
import pathlib
import shlex
import subprocess
import sys

import click

logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def check_if_conda_env_is_base():
    return True if os.environ.get("CONDA_DEFAULT_ENV") == "base" else False


def check_in_conda_env():
    return True if os.environ.get("CONDA_DEFAULT_ENV", None) else False


def check_compatible_platform():
    plataform = sys.platform
    return True if any([i in plataform for i in ["darwin", "linux"]]) else False


@click.command()
def create_conda_env():
    if not check_compatible_platform():
        raise "The OS is incompatible with the package please use either mac or linux"
    path_to_conda_env_yml = str(pathlib.Path(__file__).parent / "config" / "conda_env.yml")
    in_conda_env = check_in_conda_env()
    if not in_conda_env:
        raise ValueError("Conda is not active. Be sure to have conda active and the conda init was run.")
    else:
        fextract_already_exists = pathlib.Path(os.environ.get("CONDA_PREFIX")) / "envs" / "@lbfextract"
        if fextract_already_exists:
            cmd = f"{os.environ.get('CONDA_EXE')} env update --name @lbfextract --file {path_to_conda_env_yml} --prune"
        else:
            cmd = f"{os.environ.get('CONDA_EXE')} env create -n @lbfextract --file={path_to_conda_env_yml}"

    try:
        log.info(cmd)
        subprocess.run(shlex.split(cmd), check=True)
    except subprocess.CalledProcessError:
        log.error(f"creation of conda env with command:\n{cmd}\nfailed")
        raise
