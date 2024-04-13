import signal
import logging

from lbfextract.utils import signal_handler

logger = logging.getLogger(__name__)

signal.signal(signal.SIGINT, signal_handler)
signal.signal(signal.SIGTERM, signal_handler)

from lbfextract.cli import cli

if __name__ == "__main__":
    cli()
