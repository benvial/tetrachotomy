#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import importlib.metadata as metadata

data = metadata.metadata("tetrachotomy")
__version__ = metadata.version("tetrachotomy")
__author__ = metadata.metadata("tetrachotomy").get("Author-email")
__description__ = metadata.metadata("tetrachotomy").get("Summary")

try:
    import numdiff
    from numdiff import *
    from numdiff import _reload_package

    def set_backend(BACKEND):
        numdiff.set_backend(BACKEND)
        _reload_package("tetrachotomy")

    def get_backend():
        return numdiff.get_backend()

    BACKEND = numdiff.BACKEND

except:
    import numpy as backend

    def set_backend(BACKEND):
        pass

    def get_backend():
        return "numpy"

    BACKEND = "numpy"

from .log import *
from .tetrachotomy import *
