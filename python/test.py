import numpy as np
pi=np.pi
import matplotlib.pyplot as plt
plt.ion()
import seaborn as sns
sns.set_context("poster", font_scale = 1.)
sns.set_style("white")

import tetrachotomy
import importlib
importlib.reload(tetrachotomy)
import time

import tmm
import tetrachotomy as tc

z0 = 1e-5 - 0.01*1j
z1 = 0.05 - 0*1j
tols = (1e-2 *(1 + 1j), 1e-2 *(1 + 1j), 1e-1 *(1 + 1j))
par_integ = (1e-6, 1e-6, 15)

#####################################
pol = 'p'
h, index = 100, 2.2
# list of layer thicknesses in nm
d_list = [np.inf, h, np.inf]
# list of refractive indices
n_list = [1, index, 1]
th_0 = 0
# lam_vac = 500
# result = tmm.coh_tmm(pol, n_list, d_list, th_0, lam_vac)['t']


def func0(z):
    lam_vac = 2*pi/z /(2*pi/h)
    result = (tmm.coh_tmm(pol, n_list, d_list, th_0, lam_vac)['t'])
    return result
func = np.vectorize(func0)
Nmax = 4
alpha = (index+1)/(index-1)
poles_analytic = 1/(h*index)*(np.arange(1,Nmax+1)*pi - 1j*np.log(alpha))




z0 = 1e-3 - 0.2*1j
z1 = 1. + 0*1j


# ##################################
# z0 = 0
# z1 = 1 + 1*1j
#
#
# np.random.seed(5)
# N = 30
# x,y = np.random.rand(N), np.random.rand(N)
# zpole = x+1j*y
# xr,yr = np.random.rand(N), np.random.rand(N)
# zres = xr+1j*yr
#
# def func(z):
#     # time.sleep(1e-3)
#     out = 1
#     for zp, zr in zip(zpole, zres):
#         out += zr/(z-zp)
#     return out
#
# ##################################

x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
xp = np.linspace(x0, x1, 101)

# plt.close('all')
# fig = plt.figure()
plt.clf()
fig = plt.gcf()
ax = fig.add_subplot(111)#, aspect='equal')
ax.set_xlim((x0,x1))
ax.set_ylim((y0,y1))
plt.xlabel(r'Re $k_n$')
plt.ylabel(r'Im $k_n$')
plt.title(r'Spectrum in the complex $k$-plane')

plt.gca().plot(np.real(poles_analytic/(2*pi/h)), np.imag(poles_analytic/(2*pi/h)), 's')
# plt.gca().plot(np.real(zpole), np.imag(zpole), 's')





#
# fplot = func(xp)
# plt.clf()
# plt.plot(xp, fplot.real, label = 'Re')
# plt.plot(xp, fplot.imag, '--', label = 'Im')
# plt.legend(loc = 0)
# plt.xlabel(r'Re $k$')
# plt.ylabel(r'$f$')
# plt.title(r'Function $f$ on the real $k$ line')
# plt.xlim((x0, x1))



# tc.compute_integral(func, z0,z1,0, par_integ = par_integ)

inv_golden_number = 2/(1+np.sqrt(5))
# ratio_range = np.linspace(0.1,0.5,30)

ratio = inv_golden_number
# ncuts = []
# for ratio in ratio_range:
# plt.cla()
ratio_re, ratio_im = ratio, ratio

poles, residues, nb_cuts  = tc.pole_hunt(func, z0, z1, tols = tols,
ratio_re =ratio_re, ratio_im=ratio_im, par_integ= par_integ, poles = [], residues = [], nb_cuts = 0)
print('poles = ', poles)
print('residues = ', residues)
print('nb_cuts = ', nb_cuts)
# ncuts.append(nb_cuts)
#

#
# plt.figure()
# plt.plot(ratio_range, ncuts)

# plt.clf()
# plt.plot(np.real(poles), np.imag(poles), 'oy')
# plt.xlabel(r'Re $k_n$')
# plt.ylabel(r'Im $k_n$')
# plt.title(r'Spectrum in the complex $k$-plane')
#
# plt.xlim((x0, x1))
# plt.ylim((y0, y1))
#

#
# err_rel_re = 1-poles.real/poles_analytic.real
# err_rel_im = 1-poles.imag/poles_analytic.imag
# print('relative errors on real part: ', err_rel_re)
# print('relative errors on imaginary part: ', err_rel_im)
