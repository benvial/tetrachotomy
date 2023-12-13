import tetrachotomy as tc
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(13)
plt.ion()

z0 = -1 - 1 * 1j
z1 = 1 + 1 * 1j
x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag

Np = 10
poles_test = (-0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np)) * 2
residues_test = (
    (-0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np))
    * np.random.rand(1)
    * 100
)

isort = np.argsort(poles_test)
poles_test = poles_test[isort]
residues_test = residues_test[isort]



def func(z):
    N = len(poles_test)
    out = np.exp(10*z**2)
    for i in range(N):
        p = poles_test[i]
        r = residues_test[i]
        out += r / (z - p)
    return out


xx = np.linspace(x0, x1, 1000)
yy = np.linspace(y0, y1, 1000)
X, Y = np.meshgrid(xx, yy)
Z = X + 1j * Y
mapf = func(Z)

fig = plt.figure()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
plt.xlabel(r"Re $z$")
plt.ylabel(r"Im $z$")
plt.axis("scaled")
cm = ax.pcolormesh(xx, yy, np.log(np.abs(mapf)))
ax.plot(np.real(poles_test), np.imag(poles_test), "sk")
plt.colorbar(cm)


poles, residues, nb_cuts = tc.find_poles(
    func,
    z0,
    z1,
)


print("poles = ", poles)
print("residues = ", residues)
print("nb_cuts = ", nb_cuts)


plt.plot(np.real(poles), np.imag(poles), "or")

zeros, residues_zeros, nb_cuts_zeros = tc.find_zeros(
    func,
    z0,
    z1,
)

print("zeros = ", zeros)
print("residues zeros= ", residues_zeros)
print("nb_cuts zeros = ", nb_cuts_zeros)
plt.plot(np.real(zeros), np.imag(zeros), "+b")

assert np.allclose(poles,poles_test)
assert np.allclose(residues,residues_test)



zeros_test = (-0.5 - 0.5 * 1j + np.random.rand(Np) + 1j * np.random.rand(Np)) * 2

isort = np.argsort(zeros_test)
zeros_test = zeros_test[isort]

def func2(z):
    N = len(poles_test)
    out = 1
    for i in range(N):
        zero = zeros_test[i]
        pole = poles_test[i]
        r = residues_test[i]
        out *= (z - zero) / (z - pole)
    return out


mapf = func2(Z)

fig = plt.figure()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
plt.xlabel(r"Re $z$")
plt.ylabel(r"Im $z$")
plt.axis("scaled")
cm = ax.pcolormesh(xx, yy, np.log(np.abs(mapf)),cmap="inferno")
ax.plot(np.real(poles_test), np.imag(poles_test), "sk")
ax.plot(np.real(zeros_test), np.imag(zeros_test), "^g")
plt.colorbar(cm)


poles, residues, nb_cuts = tc.find_poles(
    func2,
    z0,
    z1,
)


print("poles = ", poles)
print("residues = ", residues)
print("nb_cuts = ", nb_cuts)


plt.plot(np.real(poles), np.imag(poles), "or")

zeros, residues_zeros, nb_cuts_zeros = tc.find_zeros(
    func2,
    z0,
    z1,
)

print("zeros = ", zeros)
print("residues zeros= ", residues_zeros)
print("nb_cuts zeros = ", nb_cuts_zeros)
plt.plot(np.real(zeros), np.imag(zeros), "+b")

assert np.allclose(poles,poles_test)
assert np.allclose(zeros,zeros_test)

