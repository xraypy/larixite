#!/usr/bin/env python
"""
site configuration
"""
import os
import sys
import time
from pathlib import Path
from packaging import version as pkg_version

from pyshortcuts import get_homedir

home_dir = get_homedir()
user_folder = Path(home_dir, '.larch').absolute().as_posix()

def bytes2str(s):
    if isinstance(s, str):
        return s
    elif isinstance(s, bytes):
        return s.decode(sys.stdout.encoding)
    return str(s, sys.stdout.encoding)


def strict_ascii(s, replacement='_'):
    """for string to be truly ASCII with all characters below 128"""
    t = bytes(s, 'UTF-8')
    return ''.join([chr(a) if a < 128 else replacement for a in t])


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
            raise FileExistsError(f"'{name}' is a file, cannot make folder with that name")
    else:
        os.makedirs(name, mode=mode)


def isotime(t=None):
    if t is None:
        t = time.time()
    sout = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t))
    return sout

def version_ge(v1, v2):
    "returns whether version string 1 >= version_string2"
    return pkg_version.parse(v1) >= pkg_version.parse(v2)
