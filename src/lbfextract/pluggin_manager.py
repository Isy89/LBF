import pluggy

import lbfextract.hookspecs
import lbfextract.hookspecs
from lbfextract.fextract.cli_lib import CliHook, CliHookExtractCoverage


def get_pluggin_manager():
    pm = pluggy.PluginManager("lbfextract_cli")
    pm.add_hookspecs(lbfextract.hookspecs.CliHookSpec)
    pm.register(CliHook())
    pm.register(CliHookExtractCoverage())
    pm.load_setuptools_entrypoints("lbfextract_cli")
    return pm
