#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: GPLv3

from PIL import Image
import numpy as np
import glob
import os

def animate_field(filename="animation.gif", **kwargs):
    tmpdir = "."
    fp_in = f"{tmpdir}/animation_tmp_*.png"

    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    img.save(
        fp=filename,
        format="GIF",
        append_images=imgs,
        save_all=True,
        duration=100,
        loop=0,
    )
    os.system(f"rm -f {tmpdir}/animation_tmp_*.png")
