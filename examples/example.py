import sys

import matplotlib.pyplot as plt
import numpy as np

import tetrachotomy as tc

# tc.set_backend(sys.argv[1])
tc.set_backend("numpy")
bk = tc.backend

plt.close("all")

np.random.seed(13)
plt.ion()

z0 = -1 - 1 * 1j
z1 = 1 + 1 * 1j
x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag

Np = 5
poles_test = (-0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np)) * 2
residues_test = -0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np)
poles_test = bk.array(poles_test)
residues_test = bk.array(residues_test)

isort = bk.argsort(poles_test.real)
poles_test = poles_test[isort]
residues_test = residues_test[isort]

order = 2


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


xx = bk.linspace(x0, x1, 300)
yy = bk.linspace(y0, y1, 300)
X, Y = bk.meshgrid(xx, yy, indexing="xy")
Z = X + 1j * Y
mapf = func(Z)
delta = 0.01

fig = plt.figure()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0 - delta, x1 + delta))
ax.set_ylim((y0 - delta, y1 + delta))
plt.xlabel(r"Re $z$")
plt.ylabel(r"Im $z$")
plt.axis("scaled")
cm = ax.imshow(
    bk.log(bk.abs(mapf)),
    extent=(z0.real, z1.real, z0.imag, z1.imag),
    origin="lower",
    cmap="inferno",
)
# cm = ax.pcolormesh(xx, yy, bk.log(bk.abs(mapf)),cmap="inferno")
# ax.plot(bk.real(poles_test), bk.imag(poles_test), "sk")
plt.colorbar(cm)
plt.tight_layout()

options = tc.get_options()
options["plot"]["rectangle"] = True
# options["plot"]["circle"] = True
options["plot"]["poles"] = True
options["plot"]["color"] = "#4d9aff"
options["plot"]["marker"]["value"] = "o"
options["plot"]["marker"]["size"] = 4

# q = (-1 + 5**0.5) / 2
q = 0.5
options["ratios"]["real"] = q
options["ratios"]["imag"] = 1 - q
options["ratios"]["circ"] = 0.5
options["refine"]["nref_max"] = 1
options["refine"]["method"] = "Nelder-Mead"
tolerances = options["tolerances"]
tol_moments = tolerances["moments"]
tol_moments["0"] = 1e-6
tol_moments["1"] = 1e-6
tol_moments["ratio"] = 1e-5

tolerances["pole"] = 1e-12
tolerances["residue"] = 1e-12
tolerances["integration"]["tol"] = 1e-8
tolerances["integration"]["divmax"] = 15


res = tc.find_poles(
    func,
    z0,
    z1,
    order=order,
    options=options,
)

poles = res.poles
residues = res.residues


print("error poles = ", bk.abs(poles - poles_test))
print("error residues = ", bk.abs(residues - residues_test))


assert bk.allclose(poles, poles_test)

assert bk.allclose(residues, residues_test)


# plt.plot(bk.real(poles), bk.imag(poles), "or")v


options["plot"]["color"] = "#41c76e"
options["plot"]["marker"]["value"] = "+"
options["plot"]["marker"]["size"] = 6
res = tc.find_zeros(
    func,
    z0,
    z1,
    order=order,
    options=options,
)

zeros = res.poles
residues_zeros = res.residues
print("zeros = ", zeros)
print("residues zeros= ", residues_zeros)
sys.exit(0)
# plt.plot(bk.real(zeros), bk.imag(zeros), "+b")


zeros_test = (-0.5 - 0.5 * 1j + bk.random.rand(Np) + 1j * bk.random.rand(Np)) * 2

isort = bk.argsort(zeros_test)
zeros_test = zeros_test[isort]


# def func2(z):
#     N = len(poles_test)
#     out = 1
#     for i in range(N):
#         zero = zeros_test[i]
#         pole = poles_test[i]
#         r = residues_test[i]
#         out *= (z - zero) ** order / (z - pole) ** order
#     return out


# mapf = func2(Z)

# fig = plt.figure()
# ax = fig.add_subplot(111)  # , aspect='equal')
# ax.set_xlim((x0, x1))
# ax.set_ylim((y0, y1))
# plt.xlabel(r"Re $z$")
# plt.ylabel(r"Im $z$")
# plt.axis("scaled")
# cm = ax.pcolormesh(xx, yy, bk.log(bk.abs(mapf)), cmap="inferno")
# ax.plot(bk.real(poles_test), bk.imag(poles_test), "sk")
# ax.plot(bk.real(zeros_test), bk.imag(zeros_test), "^g")
# plt.colorbar(cm)


# poles, residues, nb_cuts = tc.find_poles(
#     func2,
#     z0,
#     z1,
#     order=order,
#     options=options,
# )


# print("poles = ", poles)
# print("residues = ", residues)
# print("nb_cuts = ", nb_cuts)


# assert bk.allclose(poles, poles_test)


# plt.plot(bk.real(poles), bk.imag(poles), "or")

# zeros, residues_zeros, nb_cuts_zeros = tc.find_zeros(
#     func2,
#     z0,
#     z1,
#     order=order,
#     options=options,
# )

# print("zeros = ", zeros)
# print("residues zeros= ", residues_zeros)
# print("nb_cuts zeros = ", nb_cuts_zeros)
# plt.plot(bk.real(zeros), bk.imag(zeros), "+b")

# assert bk.allclose(zeros, zeros_test)


def func3(z):
    N = len(poles_test)
    out = 1e-12
    for i in range(N):
        zero = zeros_test[i]
        # pole = poles_test[i]
        # r = residues_test[i]
        out *= (z - zero) ** 2
    return out


def func3(z):
    N = len(poles_test)
    out = 0
    for i in range(N):
        zero = zeros_test[i]
        out += 1 / (z - zero) ** order
    return 1 / out


mapf = func3(Z)

fig = plt.figure()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
plt.xlabel(r"Re $z$")
plt.ylabel(r"Im $z$")
plt.axis("scaled")
cm = ax.pcolormesh(xx, yy, bk.log(bk.abs(mapf)), cmap="inferno")
# ax.plot(bk.real(poles_test), bk.imag(poles_test), "sk")
ax.plot(bk.real(zeros_test), bk.imag(zeros_test), "^g")
plt.colorbar(cm)

zeros, residues_zeros, nb_cuts_zeros = tc.find_zeros(
    func3,
    z0,
    z1,
    order=order,
    options=options,
)

print("zeros = ", zeros)
print("residues zeros= ", residues_zeros)
print("nb_cuts zeros = ", nb_cuts_zeros)
plt.plot(bk.real(zeros), bk.imag(zeros), "+b")

assert bk.allclose(zeros, zeros_test)
