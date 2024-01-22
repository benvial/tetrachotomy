#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: GPLv3

import sys

import matplotlib.pyplot as plt
import numpy as np

import tetrachotomy as tc

bk = tc.backend

# plt.ion()
# plt.close("all")

np.random.seed(int(sys.argv[1]))
z0 = -1 - 1 * 1j
z1 = 1 + 1 * 1j
x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag

Np = 4
order = 1
poles_test = (-0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np)) * 2
residues_test = -0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np)
poles_test = bk.array(poles_test)
residues_test = bk.array(residues_test)


def func(z):
    z = bk.array(z)
    N = len(poles_test)
    out = bk.exp(4 * z**2) + 0j
    out = 1
    for i in range(N):
        p = poles_test[i]
        r = residues_test[i]
        order_ = 1 if i == 0 else order
        out += r / (z - p) ** (order_)
    return out


npts = 100
xx = bk.linspace(x0, x1, npts)
yy = bk.linspace(y0, y1, npts)
X, Y = bk.meshgrid(xx, yy, indexing="xy")
Z = X + 1j * Y
mapf = func(Z)
delta = 0.01

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111)
ax.set_xlim((x0 - delta, x1 + delta))
ax.set_ylim((y0 - delta, y1 + delta))
plt.axis("scaled")
cm = ax.contourf(xx, yy, bk.log(bk.abs(mapf)), cmap="inferno", levels=9)
plt.axis("off")
plt.tight_layout()
options = tc.get_options()
options["plot"]["rectangle"] = True
options["plot"]["poles"] = False
options["plot"]["remove"] = False
options["plot"]["color"] = "#000000"

res = tc.find_poles(
    func,
    z0,
    z1,
    order=order,
    options=options,
)

poles = res.poles
residues = res.residues

plt.savefig("logo.png")
plt.savefig("logo.svg")
