#!/usr/bin/env python
"""
Various utilities for Larixite
"""

import os
import logging
import io
from typing import Union
from gzip import GzipFile
from pathlib import Path
from packaging import version as pkg_version
from charset_normalizer import from_bytes

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


def pprint(matrix):
    """Pretty print a list of lists
    from: https://stackoverflow.com/questions/13214809/pretty-print-2d-list
    an alternative could be: https://pypi.org/project/tabulate/
    """
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = "\t".join("{{:{}}}".format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print("\n".join(table))


class ColorFormatter(logging.Formatter):
    """Colored logging formatter intended for the console output"""

    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    green = "\x1b[32;21m"
    reset = "\x1b[0m"
    format = "[%(name)-15s | %(levelname)-8s] %(message)s"  #: '%(asctime)s' redundant

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


def get_logger(name="larixite", level="INFO"):
    """Get a custom logger for larixite"""
    logger = logging.getLogger(name)
    logger.setLevel(log_levels[level])
    formatter = ColorFormatter()
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(handler)
    return logger


def show_loggers(clear_handlers=False):
    """Show all existing loggers with their handlers (with the possibility to clear them)"""
    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    for logger in loggers:
        print(f"-> Logger: {logger.name}")
        if clear_handlers:
            logger.handlers.clear()
        for handler in logger.handlers:
            print(f"+-> Handler: {handler}")
            print(f"+-> Level: {handler.level} ({logging.getLevelName(handler.level)})")


def bytes2str(bytedata):
    "decode bytes using charset_normalizer.from_bytes"
    return str(from_bytes(bytedata).best())


def get_path(val: Union[Path, str, bytes]):
    """return best guess for a posix-style path name from an input string"""
    if isinstance(val, bytes):
        val = bytes2str(val)
    if isinstance(val, str):
        val = Path(val)
    return val.absolute()


def is_gzip(filename):
    "is a file gzipped?"
    with open(get_path(filename), "rb") as fh:
        return fh.read(3) == b"\x1f\x8b\x08"
    return False


def read_textfile(filename: Union[Path, io.IOBase, str, bytes], size=None) -> str:
    """read text from a file as string (decoding from bytes)

    Argument
    --------
    filename  (Path, str, bytes, or open File): file or file-like object to read
    size  (int or None): number of bytes to read

    Returns
    -------
    text of file as string.

    Notes
    ------
    1. the file encoding is detected with charset_normalizer.from_bytes
       which is then used to decode bytes read from file.
    2. line endings are normalized to be '\n', so that
       splitting on '\n' will give a list of lines.
    3. if filename is given, it can be a gzip-compressed file
    """
    text = ""

    if isinstance(filename, io.IOBase):
        text = filename.read(size)
        if filename.mode == "rb":
            text = bytes2str(text)
    else:
        fopen = GzipFile if is_gzip(filename) else open
        with fopen(get_path(filename), "rb") as fh:
            text = bytes2str(fh.read(size))
    return text.replace("\r\n", "\n").replace("\r", "\n")
