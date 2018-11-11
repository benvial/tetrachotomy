
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

plt.ion()
sns.set_context("poster", font_scale=1.)
sns.set_style("white")
pi = np.pi
plot_rect = True
plot_circ = True
plot_poles = True

tols = (1e-6 * (1 + 1j), 1e-6 * (1 + 1j), 1e-4 * (1 + 1j))
par_integ = (1e-8, 1e-8, 15)
tol_pol = 1e-10 * (1 + 1j)
tol_res = 1e-10 * (1 + 1j)
tol = 1e3 * np.finfo(float).eps


def romberg(f, a, b, divmax=10, tol_re=1e-6, tol_im=1e-6, size=1, verbose=0):
    r = np.zeros((divmax + 2, divmax + 2, size), dtype=complex)
    h = b - a
    r[0, 0, :] = 0.5 * h * (f(a) + f(b))

    powerOf2 = 1
    for i in range(1, divmax + 2):
        h = 0.5 * h
        sum = 0.0
        powerOf2 = 2 * powerOf2
        for k in range(1, powerOf2, 2):
            sum = sum + f(a + k * h)
        r[i, 0, :] = 0.5 * r[i - 1, 0, :] + sum * h

        powerOf4 = 1
        for j in range(1, i + 1):
            powerOf4 = 4 * powerOf4
            r[i, j, :] = r[i, j - 1, :] + \
                (r[i, j - 1, :] - r[i - 1, j - 1, :]) / (powerOf4 - 1)

        current_tol_re = np.abs(np.real(r[i, j, :] - r[i - 1, j - 1, :]))
        current_tol_im = np.abs(np.imag(r[i, j, :] - r[i - 1, j - 1, :]))
        # current_rtol_re = np.abs(np.real(r[i, j,:] - r[i-1, j-1,:])/np.real(r[i, j,:]))
        # current_rtol_im = np.abs(np.imag(r[i, j,:] - r[i-1, j-1,:])/np.imag(r[i, j,:]))
        # print('current_tol = ', current_tol)
        # print('current_rtol = ', current_rtol)
        if i > 0 and (np.all(current_tol_re < tol_re) or np.all(current_tol_im < tol_im)):
            #  or np.all(current_rtol_re < rtol) or np.all(current_rtol_im < rtol)  ):
            return r[i, i, :]
    if verbose:
        print('Accuracy warning: divmax exceeded, tol_re = {}, tol_im = {}'.format(
            current_tol_re, current_tol_im))
    return r[i, i, :]


def trace_rect(ax, z0, z1, **kwargs):
    patch = patches.Rectangle(
        (z0.real, z0.imag),   # (x,y)
        (z1 - z0).real,          # width
        (z1 - z0).imag,          # height
        **kwargs
    )
    ax.add_patch(patch)
    plt.pause(0.0001)
    return patch


def trace_circ(ax, zp, r, **kwargs):
    patch = patches.Circle(
        (zp.real, zp.imag),   # center (x,y)
        r,          # radius
        **kwargs
    )
    ax.add_patch(patch)
    plt.pause(0.0001)
    return patch


def parametrization_line(z0, z1, loc):
    x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
    if loc == 0:
        a, b = 0, pi / 2
        A, B = 2 / pi * (x1 - x0), z0
    elif loc == 1:
        a, b = pi / 2, pi
        A, B = 1j * 2 / pi * (y1 - y0), x1 + 1j * (2 * y0 - y1)
    elif loc == 2:
        a, b = pi, 3 * pi / 2
        A, B = -2 / pi * (x1 - x0), 3 * x1 - 2 * x0 + 1j * y1
    elif loc == 3:
        a, b = 3 * pi / 2, 2 * pi
        A, B = -1j * 2 / pi * (y1 - y0), x0 + 1j * (4 * y1 - 3 * y0)
    else:
        raise TypeError('loc must be an integer betwwen 0 and 3')
    return a, b, A, B


def poly_period(x, a, b):
    P = ((a - x)**3 * (a**2 - 5 * a * b + 3 * a * x +
                       10 * b**2 - 15 * b * x + 6 * x**2)) / (a - b)**5
    dP = -(30 * (a - x)**2 * (b - x)**2) / (a - b)**5
    return P, dP


def get_integrands(a, b, A, B, x, func):
    P, dP = poly_period(x, a, b)
    Pt, dPt = (b - a) * P + a, (b - a) * dP
    PT = A * Pt + B
    fper = A * func(PT) * dPt
    return np.array([fper, PT * fper, PT**2 * fper])


def cut_in_four(z0, z1, ratio_re=0.5, ratio_im=0.5):
    x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
    xm, ym = x0 + ratio_re * (x1 - x0), y0 + ratio_im * (y1 - y0)
    zm = xm + 1j * ym
    zb = xm + 1j * y0
    zt = xm + 1j * y1
    zl = x0 + 1j * ym
    zr = x1 + 1j * ym
    return (z0, zm), (zb, zr), (zm, z1), (zl, zt)


def compute_integral(func, z0, z1, loc, par_integ=par_integ):
    a, b, A, B = parametrization_line(z0, z1, loc)
    Ire, Iim = np.zeros(3), np.zeros(3)
    tol_re, tol_im, divmax = par_integ

    def f(x): return get_integrands(a, b, A, B, x, func)
    I = romberg(f, a, b,  tol_re=tol_re, tol_im=tol_im, divmax=divmax, size=3)
    return I / (2 * 1j * pi)


def ispole(func, z0, z1, tols=tols, par_integ=par_integ):
    I = 0
    for loc in range(4):
        I += compute_integral(func, z0, z1, loc, par_integ=par_integ)
    R10 = I[1] / I[0]
    R21 = I[2] / I[1]
    # print('I = ', I)
    # print('R10 = ', R10)
    # print('R21 = ', R21)
    r0, r1, r2 = False, False, False
    pole, residue = [], []
    if abs(I[0]) < tols[0].real and abs(I[1]) < tols[1].real:
        # if abs(I[0].real) < tols[0].real and abs(I[0].imag) < tols[0].imag and abs(I[1].real) < tols[1].real and abs(I[1].imag) < tols[1].imag:
        message = 'No poles'
        r0 = True
    elif abs(R21 - R10) / abs(R10) < tols[2].real:
        # elif abs((R21-R10).real/R10.real) < tols[2].real and (abs((R21-R10).imag/R10.imag) < tols[2].real):
        message = 'One pole'
        pole, residue = R10, I[0]
        r1 = True
    else:
        message = 'Several poles'
        r2 = True
    return pole, residue, r0, r1, r2, message


def get_integrands_circ(t, zp, r, func):
    Z = zp + r * np.exp(1j * t)
    K = 1j * r * np.exp(1j * t)
    fr = K * func(Z)
    return np.array([fr, Z * fr])


def compute_integral_circ(func, zp, r):
    tol_re, tol_im, divmax = par_integ

    def f(t): return get_integrands_circ(t, zp, r, func)
    I = romberg(f, 0, 2 * pi,  tol_re=tol_re,
                tol_im=tol_im, divmax=divmax, size=2)
    return I / (2 * 1j * pi)


def min_dist(z0, z1, zp):
    x0, x1, y0, y1, xp, yp = z0.real, z1.real, z0.imag, z1.imag, zp.real, zp.imag
    return np.min(np.array([xp - x0, yp - y0, x1 - xp, y1 - yp]))


def refine_pole(func, z0, z1, zp, zr, tol_pol=tol_pol, tol_res=tol_res,
                nref_max=100, ratio_circ=0.5, verbose=0):
    r = min_dist(z0, z1, zp)
    if verbose:
        print('refining pole')
    conv_pol_re, conv_pol_im, conv_res_re, conv_res_im = 1, 1, 1, 1
    nref = 0
    while nref < nref_max and (conv_pol_re > tol_pol.real or
                               conv_pol_im > tol_pol.imag or conv_res_re > tol_res.real or conv_res_im > tol_res.imag):
        if plot_circ:
            circ_b = trace_circ(plt.gca(), zp, r, fill=False,
                                linewidth=3,  edgecolor="#ff884d")
            circ = trace_circ(plt.gca(), zp, r, fill=True,
                              linewidth=0, facecolor="#ff884d", alpha=0.2)
        I = compute_integral_circ(func, zp, r)
        pole_ref, residue_ref = I[1] / I[0], I[0]
        conv_pol_re = np.abs(np.real(zp - pole_ref))
        conv_pol_im = np.abs(np.imag(zp - pole_ref))
        conv_res_re = np.abs(np.real(zr - residue_ref))
        conv_res_im = np.abs(np.imag(zr - residue_ref))
        if verbose:
            print('Refined pole {}, with residue {}'.format(pole_ref, residue_ref))
            print('conv_pol_re = ', conv_pol_re)
            print('conv_pol_im = ', conv_pol_im)
            print('conv_res_re = ', conv_res_re)
            print('conv_res_im = ', conv_res_im)
        r = r * ratio_circ
        zp = pole_ref
        zr = residue_ref
        if plot_circ:
            circ.remove()
            circ_b.remove()
        nref += 1
        if nref >= nref_max:
            if verbose:
                print('Accuracy warning, nref_max exceeded')
    return pole_ref, residue_ref


def pole_hunt(func, z0, z1, tols=tols, ratio_re=0.5, ratio_im=0.5, nref_max=100, ratio_circ=0.5,
              tol_pol=tol_pol, tol_res=tol_res,
              par_integ=par_integ, poles=[], residues=[], nb_cuts=0, verbose=0):

    if plot_rect:
        trace_rect(plt.gca(), z0, z1, fill=True, linewidth=0,
                   facecolor="#ff884d", alpha=0.2)
        trace_rect(plt.gca(), z0, z1, fill=False,
                   linewidth=3,  edgecolor="#ff884d")
        trace_rect(plt.gca(), z0, z1, fill=True, facecolor="w",
                   linewidth=3, edgecolor="#000000")
    pole, residue, r0, r1, r2, message = ispole(
        func, z0, z1, tols=tols,  par_integ=par_integ)
    if verbose:
        print('nb_cuts = ', nb_cuts)
        print(message)
    if r2:
        if verbose:
            print('Cutting in four')
        nb_cuts += 1
        Z = cut_in_four(z0, z1, ratio_re=ratio_re, ratio_im=ratio_im)
        for z in Z:
            _, _, nb_cuts = pole_hunt(func, z[0], z[1], tols=tols, ratio_re=ratio_re, ratio_im=ratio_re,
                                      nref_max=nref_max, ratio_circ=ratio_circ, tol_pol=tol_pol, tol_res=tol_res,
                                      par_integ=par_integ, poles=poles, residues=residues, nb_cuts=nb_cuts)
    else:
        if r1:
            if verbose:
                print('Found a new pole {}, with residue {}'.format(pole, residue))
            pole, residue = refine_pole(func, z0, z1, pole, residue, tol_pol=tol_pol, tol_res=tol_res,
                                        nref_max=nref_max, ratio_circ=ratio_circ)
            poles.append(pole)
            residues.append(residue)
            if plot_poles:
                plt.gca().plot(np.real(pole), np.imag(pole), 'o', color="#ff884d", ms=7)
                plt.pause(0.0001)
    # trace_rect(plt.gca(), z0, z1, fill = False, linewidth=3, edgecolor="#000000")
    poles, residues = np.array(poles), np.array(residues)
    isort = poles.real.argsort()
    return poles[isort], residues[isort], nb_cuts


if __name__ == '__main__':
    print("This is the tetrachotomy module")
