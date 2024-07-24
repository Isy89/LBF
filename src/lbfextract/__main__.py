import signal
import logging

import rich_click as click
from lbfextract.utils import signal_handler
import collections.abc as cabc
import shutil


class CustomHelpFormatter(click.rich_help_formatter.RichHelpFormatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = shutil.get_terminal_size().columns
        self.max_width = 120

    def write_dl(self,
                 rows: cabc.Sequence[tuple[str, str]],
                 col_max: int = 45,
                 col_spacing: int = 2, ):
        super().write_dl(col_max=col_max, col_spacing=col_spacing, rows=rows)


logger = logging.getLogger(__name__)

signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)

from lbfextract.cli import cli
from rich_click import RichContext

RichContext.formatter_class = CustomHelpFormatter
cli.context_class = RichContext

if __name__ == "__main__":
    cli()
