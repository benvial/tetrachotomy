import numpy as np
pi=np.pi
import matplotlib.pyplot as plt
plt.ion()
import seaborn as sns
sns.set_context("poster", font_scale = 1.)
sns.set_style("white")
cmap = sns.diverging_palette(10, 220, sep=80, n=101, as_cmap = True)
import tetrachotomy
import importlib
importlib.reload(tetrachotomy)
import time
from matplotlib.ticker import MaxNLocator
import tmm
import tetrachotomy as tc


tols = (1e-6 *(1 + 1j), 1e-6 *(1 + 1j), 1e-6 *(1 + 1j))
par_integ = (1e-6, 1e-6, 11)
tol_pol = 1e-13 *(1 + 1j)
tol_res = 1e-13 *(1 + 1j)
inv_golden_number = 2/(1+np.sqrt(5))
ratio = inv_golden_number
ratio_circ = inv_golden_number
nref_max = 20
ratio_re, ratio_im = ratio, ratio
#####################################

z0 = -0.5 - 0.1*1j
z1 = 0.51 + 0.0*1j

pol = 'p'
h, index = 100, 2.2
d_list = [np.inf, h, np.inf]
n_list = [1, index, 1]

# d_list = [np.inf, 110, 40, np.inf]
# n_list = [1, 2,4, 1]

th_0 = 0

def func0(z):
    if z==0:
        z = 1e-16
    lam_vac = 2*pi/z /(2*pi/h)
    result = (tmm.coh_tmm(pol, n_list, d_list, th_0, lam_vac)['t'])
    return result

func = np.vectorize(func0)

Nmax = 2
alpha = (index+1)/(index-1)
poles_analytic = 1/(h*index)*(np.arange(-Nmax,Nmax+1)*pi - 1j*np.log(alpha))/(2*pi/h)


# ##################################
# z0 = 0
# z1 = 1 + 1*1j
#
# np.random.seed(5)
# N = 7
# x,y = np.random.rand(N), np.random.rand(N)
# zpole = x+1j*y
# xr,yr = np.random.rand(N), np.random.rand(N)
# zres = xr+1j*yr
# isort = zpole.real.argsort()
# zpole = zpole[isort]
# zres = zres[isort]
#
#
# def func(z):
#     # time.sleep(1e-15)
#     out = 0
#     for zp, zr in zip(zpole, zres):
#         out += zr/(z-zp)
#     return out
# poles_analytic = zpole
# ##################################

x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag


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

plt.gca().plot(np.real(poles_analytic), np.imag(poles_analytic), 's')
# plt.axis('equal')



poles, residues, nb_cuts  = tc.pole_hunt(func, z0, z1,tols = tols, ratio_re =ratio_re , ratio_im =ratio_re,
nref_max = nref_max, ratio_circ = ratio_circ, tol_pol = tol_pol, tol_res = tol_res,
par_integ = par_integ, poles = [], residues = [], nb_cuts = 0)
print('poles = ', poles)
print('residues = ', residues)
print('nb_cuts = ', nb_cuts)

err_re = np.abs(poles.real-poles_analytic.real)
err_im = np.abs(poles.imag-poles_analytic.imag)
err_abs = np.abs(poles_analytic-poles)
print('POLES')
print('  errors on real part: ', err_re)
print('  mean: ', np.mean(err_re))
print('  errors on imaginary part: ', err_im)
print('  mean: ', np.mean(err_im))
print('  errors module: ', err_abs)
print('  mean: ', np.mean(err_abs))


res_err_re = np.abs(residues.real-zres.real)
res_err_im = np.abs(residues.imag-zres.imag)
res_err_abs = np.abs(zres-residues)
print('RESIDUES')
print('  errors on real part: ', res_err_re)
print('  mean: ', np.mean(res_err_re))
print('  errors on imaginary part: ', res_err_im)
print('  mean: ', np.mean(res_err_im))
print('  errors module: ', res_err_abs)
print('  mean: ', np.mean(res_err_abs))


x = np.linspace(x0, x1, 101)
y = np.linspace(y0, y1, 101)
xx,yy = np.meshgrid(x, y)
zz = xx+1j*yy
fcplx = func(zz)
fplot = np.log10(np.abs(fcplx))
fplot= fcplx.imag


fig = plt.figure()
# plt.clf()
ax = fig.add_subplot(111)#, aspect='equal')
cmap = sns.diverging_palette(10, 220, sep=20, n=101, as_cmap = True)
levels = MaxNLocator(nbins=11).tick_values(fplot.min(), fplot.max()*0.1)


levels = MaxNLocator(nbins=11).tick_values(-11,11)
plt.contourf(xx, yy, fplot, cmap = cmap, levels=levels)
plt.colorbar()
plt.plot(np.real(poles), np.imag(poles), 'ok')
ax.set_xlim((x0,x1))
ax.set_ylim((y0,y1))
plt.xlabel(r'Re $k_n$')
plt.ylabel(r'Im $k_n$')
#
# xp = np.linspace(x0, x1, 101)
#
# plt.figure()
# plt.plot(ratio_range, ncuts)

# plt.clf()

# plt.xlabel(r'Re $k_n$')
# plt.ylabel(r'Im $k_n$')
# plt.title(r'Spectrum in the complex $k$-plane')
#
# plt.xlim((x0, x1))
# plt.ylim((y0, y1))
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
