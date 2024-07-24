"""
In this module the general CLI object is created, which aggregates all the command and subcommands which can be used
in LBFextract.

Here the following clis are aggregated:
    - cli_tfa
    - cli_setup
    - cli_feature_extraction
    - start_tui

"""
import typing as t

import rich_click as click

from lbfextract.autocompletion import enable_autocompletion, disable_autocompletion
from lbfextract.generate_plugin_structure import new_plugin
from lbfextract.utils import start_msg
from lbfextract.transcription_factor_analysis.differential_signal_analysis import cli as cli_tfa
from lbfextract.setup_conda_env import create_conda_env
from lbfextract.fextract.tui.app import start_tui
from lbfextract.pluggin_manager import get_pluggin_manager


class FextractCli(click.RichGroup):
    def __init__(
            self,
            name: t.Optional[str] = None,
            commands: t.Optional[t.Union[t.Dict[str, click.Command], t.Sequence[click.Command]]] = None,
            **attrs: t.Any, ):
        start_msg()
        super().__init__(name, commands, **attrs)


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(cls=FextractCli, name="main", context_settings=CONTEXT_SETTINGS)
def cli():
    """
    LBFextract CLI
    
    This is the main command group for the lbfextract CLI and serves as the entry point for all lbfextract's 
    commands.

    \b
    For a list of all available commands and their short descriptions, use:
    \n
    \b
        - lbfextract --help
    \n
    
    \b
    For detailed help on a specific command, including available arguments and their default values, use:
    \n
    \b
        - lbfextract [command] --help
    \n
    
    """
    pass


def check_conda_env_is_present():
    from lbfextract import PATH_TO_SAMTOOLS
    if PATH_TO_SAMTOOLS is None:
        raise RuntimeError("To use lbf feature extraction methods the @lbfextract conda env must be present."
                           "Please be sure to have conda active and run `lbfextract setup create-conda-env`!")


@click.group(name="feature_extraction_commands")
def cli_feature_extraction(context_settings=CONTEXT_SETTINGS):
    """
        \b
        This command group contains all commands to extract features from a BAM file. 
        \n

        \b
        To display a list of available commands and their short descriptions, use:
        \n
        \b
            lbfextract feature_extraction_commands --help
        \n

        \b
        For detailed help on a specific command, including available arguments and their default values, use:
        \n
        \b
            lbfextract feature_extraction_commands [command] --help
        \n
    """
    print("here")
    check_conda_env_is_present()


@click.group(name="setup")
def cli_setup(context_settings=CONTEXT_SETTINGS):
    """
    This command group contains all setup-related commands.
    \n

    \b
    Use this group to initialize, configure, and manage the setup processes for your application.
    \n

    \b
    For a list of available setup commands and their short descriptions, use:
    \n
    \b
        lbfextract setup --help
    \n

    \b
    For detailed help on a specific setup command, including available arguments and their default values, use:
        lbfextract setup [command] --help
    """

    pass


cli_setup.add_command(new_plugin)
cli_setup.add_command(enable_autocompletion)
cli_setup.add_command(disable_autocompletion)
cli_setup.add_command(create_conda_env)

plugin_manager = get_pluggin_manager()
for command in plugin_manager.hook.get_command():
    if isinstance(command, list):
        for c in command:
            cli_feature_extraction.add_command(c)
    elif isinstance(command, click.Command):
        cli_feature_extraction.add_command(command)
    else:
        raise TypeError(f"The command returned by the hook.get_command() must be a list of click.Command"
                        f" or a click.Command, not {type(command)}")

cli.add_command(cli_tfa)
cli.add_command(cli_setup)
cli.add_command(cli_feature_extraction)
cli.add_command(start_tui)
