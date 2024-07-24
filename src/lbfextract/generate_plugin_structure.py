import warnings

import rich_click as click
import jinja2
import shutil
import pathlib


def sanitize_name(ctx, param, value):
    """Sanitize the name of the signal and the cli command"""
    if "-" in value:
        warnings.warn("Found a dash in the name of the signal or the cli command. It  was replaced with an underscore.")
    return value.replace("-", "_")


@click.command()
@click.option("--name_of_the_signal", "-ns", type=str, help="name of signal", callback=sanitize_name)
@click.option("--name_of_the_cli_command", "-nc", type=str, help="name of cli command", callback=sanitize_name)
@click.option("--out_dir", "-od", type=click.Path(exists=False,
                                                  file_okay=False,
                                                  dir_okay=True,
                                                  writable=True,
                                                  readable=True,
                                                  resolve_path=False,
                                                  allow_dash=True,
                                                  path_type=pathlib.Path,
                                                  executable=False), default=pathlib.Path.cwd(),
              help="path to output dir where the plugin structure"
                   " will be created")
def new_plugin(name_of_the_signal: str, name_of_the_cli_command: str, out_dir: pathlib.Path):
    path_template_folder = pathlib.Path(__file__).parent / "plugin_template"
    out_dir_pkj = out_dir / f"fextract_{name_of_the_signal}"
    out_dir_pkj.mkdir(exist_ok=True, parents=True)

    template_loader = jinja2.FileSystemLoader(searchpath=path_template_folder)
    template_env = jinja2.Environment(loader=template_loader, lstrip_blocks=True, trim_blocks=True)
    setup_py_template = template_env.get_template("setup.template")
    plugin_py_template = template_env.get_template("plugin.template")
    rendered_plugin_template = plugin_py_template.render({
        "name_of_the_signal": name_of_the_signal,
        "name_of_the_cli_command": name_of_the_cli_command,
    })
    rendered_template = setup_py_template.render({
        "name_of_the_signal": name_of_the_signal,
        "name_of_the_cli_command": name_of_the_cli_command
    })

    shutil.copytree(path_template_folder / "src", out_dir_pkj / "src" / f"fextract_{name_of_the_signal}",
                    dirs_exist_ok=True)
    with open(out_dir_pkj / "setup.py", "w") as f:
        f.write(rendered_template)
    with open(out_dir_pkj / "src" / f"fextract_{name_of_the_signal}" / "plugin.py", "w") as f:
        f.write(rendered_plugin_template)
    shutil.copy2(path_template_folder / "requirements.txt", out_dir_pkj / "requirements.txt")


if __name__ == "__main__":
    new_plugin()
