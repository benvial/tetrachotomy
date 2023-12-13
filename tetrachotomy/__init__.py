#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import importlib.metadata as metadata

data = metadata.metadata("tetrachotomy")
__version__ = metadata.version("tetrachotomy")
__author__ = metadata.metadata("tetrachotomy").get("Author-email")
__description__ = metadata.metadata("tetrachotomy").get("Summary")


from .tetrachotomy import *
