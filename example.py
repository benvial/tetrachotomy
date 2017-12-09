import numpy as np
import matplotlib.pyplot as plt
import tetrachotomy
# import importlib
# importlib.reload(tetrachotomy)


plt.ion()
pi = np.pi

tols = (1e-6 * (1 + 1j), 1e-6 * (1 + 1j), 1e-6 * (1 + 1j))
par_integ = (1e-6, 1e-6, 10)
tol_pol = 1e-6 * (1 + 1j)
tol_res = 1e-6 * (1 + 1j)
inv_golden_number = 2 / (1 + np.sqrt(5))
ratio = inv_golden_number
ratio_circ = 1-inv_golden_number
nref_max = 100
ratio_re, ratio_im = ratio, ratio


z0 = -1 - 1 * 1j
z1 = 1 + 1 * 1j
x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
#####################################

# poles_test = [0.2 - 0.02*1j, 0.4 - 0.1*1j]
# residues_test = [2 + 0.1*1j, -10 - 3*1j]
np.random.seed(11)
Np = 20
poles_test = (-0.5 - 0.5*1j + np.random.rand(Np) + 1j*np.random.rand(Np))*2
residues_test = (-0.5 - 0.5*1j + np.random.rand(Np) + 1j*np.random.rand(Np))* np.random.rand(1)*100

isort = np.argsort(poles_test)
poles_test = poles_test[isort]
residues_test = residues_test[isort]


def func(z):
    N = len(poles_test)
    out = 0
    for i in range(N):
        p = poles_test[i]
        r = residues_test[i]
        out += r/(z-p)
    return out



poles, residues, nb_cuts = tetrachotomy.pole_hunt(func, z0, z1, tols=tols, ratio_re=ratio_re, ratio_im=ratio_re,
                                        nref_max=nref_max, ratio_circ=ratio_circ, tol_pol=tol_pol, tol_res=tol_res,
                                        par_integ=par_integ, poles=[], residues=[], nb_cuts=0)
print('poles = ', poles)
print('residues = ', residues)
print('nb_cuts = ', nb_cuts)



fig = plt.figure()
fig = plt.gcf()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
plt.xlabel(r'Re $z$')
plt.ylabel(r'Im $z$')
plt.gca().plot(np.real(poles_test), np.imag(poles_test), 'sk')
plt.gca().plot(np.real(poles), np.imag(poles), 'or')

err_p = (poles - poles_test)
err_r = (residues - residues_test)

erel_p_re = np.abs(err_p.real/poles_test.real)
erel_p_im = np.abs(err_p.imag/poles_test.imag)
erel_r_re = np.abs(err_r.real/residues_test.real)
erel_r_im = np.abs(err_r.imag/residues_test.imag)
print("--- RELATIVE ERRORS ---")
print("  poles")
print("     real: ", erel_p_re)
print("     imag: ", erel_p_im)
print("  residues")
print("     real: ", erel_r_re)
print("     imag: ", erel_r_im)
