#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: GPLv3

import glob
import os

import numpy as np
from PIL import GifImagePlugin, Image

# GifImagePlugin.LOADING_STRATEGY = GifImagePlugin.LoadingStrategy.RGB_ALWAYS

# User editable values
method = Image.FASTOCTREE
colors = 250


def animate_field(filename="animation.gif", **kwargs):
    tmpdir = "."
    fp_in = f"{tmpdir}/animation_tmp_*.png"

    # img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    imgs = []
    for f in sorted(glob.glob(fp_in)):
        im = Image.open(f)
        # pImage = im.quantize(colors=colors, method=method, dither=0)
        pImage = im
        imgs.append(pImage)
    imgs = [imgs[0]] * 10 + imgs[1:-1] + [imgs[-1]] * 10
    imgs[0].save(
        fp=filename,
        format="GIF",
        append_images=imgs[1:],
        save_all=True,
        duration=100,
        loop=0,
        quality=100,
    )
    os.system(f"rm -f {tmpdir}/animation_tmp_*.jpg")
    os.system(f"rm -f {tmpdir}/animation_tmp_*.png")


if __name__ == "__main__":
    animate_field()
