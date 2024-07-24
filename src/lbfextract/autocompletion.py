import os

import rich_click as click
from os.path import expanduser
import logging

logger = logging.getLogger(__name__)

home = expanduser("~")

lbfextract_completion = """
_lbfextract_completion() {
    local IFS=$'\n'
    local response

    response=$(env COMP_WORDS="${COMP_WORDS[*]}" COMP_CWORD=$COMP_CWORD _LBFEXTRACT_COMPLETE=bash_complete $1)

    for completion in $response; do
        IFS=',' read type value <<< "$completion"

        if [[ $type == 'dir' ]]; then
            COMPREPLY=()
            compopt -o dirnames
        elif [[ $type == 'file' ]]; then
            COMPREPLY=()
            compopt -o default
        elif [[ $type == 'plain' ]]; then
            COMPREPLY+=($value)
        fi
    done

    return 0
}

_lbfextract_completion_setup() {
    complete -F _lbfextract_completion lbfextract
}

_lbfextract_completion_setup;


"""


@click.command()
def enable_autocompletion():
    with open(home + '/.lbfextract-complete.bash', 'w') as f:
        f.write(lbfextract_completion)
    with open(home + '/.bashrc', 'a') as f:
        f.write('source ~/.lbfextract-complete.bash')
    logger.info("Autocompletion was enabled. To use it, open a new terminal or run `source ~/.bashrc`. "
                "In the latter case it may be needed to deactivate and reactivate your conda environment."
                "Run `conda deactivate && conda activate your_conda_env_where_LBF_is_installed`")


@click.command()
def disable_autocompletion():
    with open(home + '/.bashrc', 'r') as f:
        lines = f.readlines()

    with open(home + '/.bashrc', 'w') as f:
        for line in lines:
            if line.strip() != 'source ~/.lbfextract-complete.bash':
                f.write(line)

    os.remove(home + '/.lbfextract-complete.bash')
    logger.info("Autocompletion was disabled.")



if __name__ == "__main__":
    @click.group(name="autocompletion")
    def cli_autocompletion():
        pass


    cli_autocompletion.add_command(enable_autocompletion)
    cli_autocompletion.add_command(disable_autocompletion)
    cli_autocompletion()
