#!/usr/bin/env python
"""
Various utilities for Larixite
"""

import os
import logging
from pathlib import Path
from packaging import version as pkg_version
from pyshortcuts import get_homedir, bytes2str, isotime

HAS_IPYTHON = False
try:
    from IPython.core.interactiveshell import InteractiveShell

    HAS_IPYTHON = True
except ImportError:
    pass

# Enable color output in Jupyter Notebook
if HAS_IPYTHON:
    InteractiveShell.ast_node_interactivity = "all"


home_dir = get_homedir()
user_folder = Path(home_dir, ".larch").absolute().as_posix()


def strict_ascii(s, replacement="_"):
    """for string to be truly ASCII with all characters below 128"""
    t = bytes(s, "UTF-8")
    return "".join([chr(a) if a < 128 else replacement for a in t])


def mkdir(name, mode=0o775):
    """create directory (and any intermediate subdirectories)

    Options:
    --------
      mode   permission mask to use for creating directory (default=0775)
    """
    path = Path(name)
    if path.exists():
        if path.is_dir():
            os.chmod(name, mode)
        else:
            raise FileExistsError(
                f"'{name}' is a file, cannot make folder with that name"
            )
    else:
        os.makedirs(name, mode=mode)


def version_ge(v1, v2):
    "returns whether version string 1 >= version_string2"
    return pkg_version.parse(v1) >= pkg_version.parse(v2)


def fcompact(val):
    """format a float value, removing extra trailing zeros"""
    val = f"{val:.6f}"
    while val.endswith("0"):
        val = val[:-1]
    if val.endswith("."):
        val = val + "0"
    return val


class ColorFormatter(logging.Formatter):
    """Colored logging formatter intended for the console output"""

    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    green = "\x1b[32;21m"
    reset = "\x1b[0m"
    format = "[%(name)-s | %(levelname)-8s] %(message)s" #: '%(asctime)s' redundant

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


log_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def get_color_logger(level="INFO"):
    """Return larixite logger with color output"""
    logger = logging.getLogger("larixite")
    logger.setLevel(log_levels[level])
    formatter = ColorFormatter()
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(handler)
    return logger
