from __future__ import annotations

import asyncio
import hashlib
import logging
import pathlib
import sys
import click
from lbfextract.feature_extractor import FeatureExtractor
from rich import print
from rich.markdown import Markdown
from textual.app import App, ComposeResult
from textual.containers import Horizontal, Vertical, Container
from textual.widgets import Input, Button, Static, Header, Footer, Checkbox, Welcome
from asyncio.subprocess import PIPE, STDOUT

import lbfextract.fextract
from lbfextract.fextract.tui.my_directory_tree import DirectoryTree


def get_element_id(element_name):
    h = hashlib.sha1()
    h.update(element_name.encode())
    return element_name + "_" + h.hexdigest()[:4]


logger = logging.getLogger(__name__)


def generate_path_element(label: str):
    logger.debug(f"Generating path element for {id} label {label}")
    return Horizontal(
        Button(label=f":file_folder: {label}", id=get_element_id(label), classes="button_tree"),
        Static(label, id=f"static_{get_element_id(label)}", classes="input_element"),
        classes="element"
    )


def generate_bool_element(label: str):
    logger.debug(f"Generating path element for {id} label {label}")
    return Horizontal(
        Static(f"{label}:", classes="label"),
        Checkbox(value=False, classes="checkbox", id=f"input_{get_element_id(label)}"),
        classes="element"
    )


def generate_static(label, help):
    return Horizontal(
        Static(f"{label}:", classes="label"),
        Input(placeholder=help,
              classes="input_element", id=f"input_{get_element_id(label)}"),
        classes="element"
    )


def generate_form_from_command(command):
    """
    Generate a form from a command.
    """
    renedered_commands = []
    for param in command.params:
        if param.type.name in ["path", "file", "directory"]:
            print(param.type.name)
            renedered_commands.append(generate_path_element(param.name))
        elif isinstance(param.type, click.types.BoolParamType):
            print(param.type.name)
            renedered_commands.append(generate_bool_element(param.name))
        else:
            print(param.type.name)
            renedered_commands.append(generate_static(param.name, help=param.help))

    form = Vertical(*renedered_commands, Button(label="Submit", id="submit"))

    return form


TITLE = (
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
    f'''Version {lbfextract.__version__}''')

WELCOME_MD = r"""
## A plugin implementation of feature extraction from BAM files and BED files.

LBFextract defines a series of hooks to carry out the feature extraction process from BAM files.
It extracts signals from the intervals defined in one or multiple BED files.

## Copyright 
Original work on LBFextract and accessory code Copyright (c) 2023 Isaac Lazzeri

## Licence
GNU General Public License v3.0

## Contact

For any questions please contact: 
* <LBFextract@gmail.com>

If you find any bugs please report them here:
* <https://github.com/Isy89/LBFextract/issues> 

"""


class Welcome(Static):
    DEFAULT_CSS = """
        Welcome {
            width: 100%;
            height: 100%;
            background: $surface;
        }

        Welcome Container {
            padding: 1;
            background: $panel;
            color: $text;
        }

        Welcome #text {
            margin:  0 1;
        }

        Welcome #close {
            dock: bottom;
            width: 100%;        
        }
    """

    def compose(self) -> ComposeResult:
        yield Container(
            Static(TITLE, id="text_1"),
            Static(Markdown(WELCOME_MD), id="text_2"),
            id="md")
        yield Button("OK", id="close", variant="success")


def Scrollview(param):
    pass


class FextractApp(App):
    CSS_PATH = "app.css"
    pressed_botton = None
    show_tree = False

    BINDINGS = [
        ("ctrl+t", "toggle_tree", "Toggle Files"),
        ("q", "quit", "Quit"),
        ("ctrl+p", "toggle_processes_pane", "Toggle Processes Pane"),
        ("ctrl+r", "check_task_status", "Check Task Status"),
    ]

    def __init__(self, path=str(pathlib.Path.cwd()), signal_type: str = "extract_coverage"):
        super().__init__()
        self.fe = FeatureExtractor()
        self.signal_type = signal_type
        self.command = self.fe.extractors.get(self.signal_type, None)
        self.show_process_pane = False
        if self.command is None:
            raise ValueError(f"Signal type {self.signal_type} not implemented")
        self.path = path
        self.task_number = 0
        self.running_processes = []

    def action_toggle_tree(self) -> None:
        """Called in response to key binding."""
        self.show_tree = not self.show_tree
        if not self.show_tree:
            self.query_one("#tree-view").remove_class("-active")
            self.query_one("#tree-view").add_class("hide")
            self.query_one("#tree_footer").add_class("hide")
        else:
            self.query_one("#tree-view").add_class("-active")
            self.query_one("#tree-view").remove_class("hide")
            self.query_one("#tree_footer").remove_class("hide")

    def action_toggle_processes_pane(self) -> None:
        """Called in response to key binding."""
        self.show_process_pane = not self.show_process_pane
        if not self.show_process_pane:
            self.query_one("#process_pane-view").remove_class("-active")
            self.query_one("#process_pane-view").add_class("hide")
        else:
            self.query_one("#process_pane-view").add_class("-active")
            self.query_one("#process_pane-view").remove_class("hide")

    def action_check_task_status(self) -> None:
        """Called in response to key binding."""
        for count, process in enumerate(self.running_processes):
            if process.returncode is None:
                if process.terminate():
                    print(f"Task {count} terminated")
                    task_static = self.query_one(f"#task_{count}")
                    task_static.renderable = f"Task {count} terminated"
                    task_static.remove_class("static_task")
                    task_static.add_class("terminated_task")
                else:
                    print(f"Task {count} still running")
            else:
                print(f"Task {count} terminated")
                task_static = self.query_one(f"#task_{count}")
                task_static.renderable = Markdown(f"Task {count} terminated")
                task_static_conatiner = self.query_one(f"#container_task_{count}")
                task_static_conatiner.remove_class("static_task")
                task_static_conatiner.add_class("terminated_task")

    def compose(self) -> ComposeResult:
        yield Header()
        yield Horizontal(
            Welcome(id="welkome"),
            Container(
                DirectoryTree(path=self.path),
                id="tree-view", classes="hide"
            ),
            Container(
                id="process_pane-view", classes="hide"
            ),
            generate_form_from_command(self.command),
        )
        yield Footer()
        yield Container(Static(id="static_tree_footer"),
                        id="tree_footer", classes="hide")

    async def run_command(self):

        dict_args_val = {}
        for key in self.command.params:
            if key.is_flag:
                dict_args_val[f"--{key.name}"] = "" if self.query_one(
                    f"#input_{get_element_id(key.name)}").value else None
            elif "path" in key.name:
                dict_args_val[f"--{key.name}"] = self.query_one(f"#static_{get_element_id(key.name)}").renderable
            else:
                dict_args_val[f"--{key.name}"] = self.query_one(f"#input_{get_element_id(key.name)}").value

        dict_args_val = {
            k: v for k, v in dict_args_val.items() if v
        }

        args_cmd = " ".join([f"lbfextract feature_extraction_commands {self.signal_type.replace('_', '-')}"] + [
            f"{k} {v}" for k, v in dict_args_val.items() if v or v != ""
        ])

        print(args_cmd)
        try:
            process = await asyncio.create_subprocess_shell(args_cmd, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
            self.running_processes.append(process)
            await self.query_one("#process_pane-view").mount(
                Container(
                    Static(id=f"task_{self.task_number}",
                           renderable=Markdown(
                               f"Task {self.task_number} running...")),
                    id=f"container_task_{self.task_number}",
                    classes="static_task"))
            self.task_number += 1
        except Exception as e:
            print(e)
        return

    async def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "submit":
            button = self.query_one(f"#{event.button.id}")
            button.remove_class("pressed")
            button.label = "Submitted :crossed_fingers:"
            await self.run_command()
            button.disabled = False

            self.query_one(f"#{event.button.id}").add_class("submit")
        elif event.button.id == "close":
            self.query_one(f"#welkome").remove()
        else:
            for button in self.query("Button"):
                button.remove_class("pressed")
            self.query_one(f"#{event.button.id}").add_class("pressed")
            self.pressed_botton = event.button.id
            self.query_one("#tree-view").add_class("-active")
            self.query_one("#tree-view").remove_class("hide")

            self.query_one("#tree_footer").remove_class("hide")
            self.query_one(f"#tree-view").focus()

    def on_input_changed(self, event: Input.Changed) -> None:
        self.query_one("#submit").disabled = False

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        self.query_one("#submit").disabled = False

    def on_directory_tree_file_click(self, event: DirectoryTree.FileClick) -> None:
        # self.query_one("#tree-view").remove_class("hide")
        print(f"botton: {self.pressed_botton}")
        button_id = f"#static_{self.pressed_botton}"
        print(f"button_id: {button_id}")
        static_path = self.query_one(button_id)
        static_path.renderable = str(event.path)
        static_path.refresh()
        tree_foter = self.query_one("#static_tree_footer")
        tree_foter.renderable = str(event.path)
        tree_foter.refresh()
        # self.query_one("#tree-view").add_class("hide")
        self.query_one(f"#{self.pressed_botton}").remove_class("pressed")


@click.command()
@click.option('--path_to_root_dir', type=click.Path(exists=False,
                                                    file_okay=True,
                                                    dir_okay=True,
                                                    writable=False,
                                                    readable=True,
                                                    resolve_path=False,
                                                    allow_dash=True,
                                                    executable=False),
              help='path defining the root of the directory tree '
                   'showed in the sidebar')
@click.option('--signal_type', type=str, default="extract_coverage",
              help='signal type to extract from the bam file')
def start_tui(path_to_root_dir: str, signal_type: str):
    app = FextractApp(path=path_to_root_dir, signal_type=signal_type)
    sys.exit(app.run())


if __name__ == "__main__":
    start_tui()
