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

import click

from lbfextract.autocompletion import enable_autocompletion, disable_autocompletion
from lbfextract.generate_plugin_structure import new_plugin
from lbfextract.utils import start_msg
from lbfextract.transcription_factor_analysis.differential_signal_analysis import cli as cli_tfa
from lbfextract.setup_conda_env import create_conda_envs
from lbfextract.fextract.tui.app import start_tui
from lbfextract.pluggin_manager import get_pluggin_manager


class FextractCli(click.Group):
    def __init__(
            self,
            name: t.Optional[str] = None,
            commands: t.Optional[t.Union[t.Dict[str, click.Command], t.Sequence[click.Command]]] = None,
            **attrs: t.Any, ):
        start_msg()
        super().__init__(name, commands, **attrs)


@click.group(cls=FextractCli, name="main")
def cli():
    pass


@click.group(name="feature_extraction_commands")
def cli_feature_extraction():
    pass


@click.group(name="setup")
def cli_setup():
    pass


cli_setup.add_command(new_plugin)
cli_setup.add_command(enable_autocompletion)
cli_setup.add_command(disable_autocompletion)
cli_setup.add_command(create_conda_envs)

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
