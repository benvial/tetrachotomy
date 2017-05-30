
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import seaborn as sns
sns.set_context("poster", font_scale = 1.)
sns.set_style("white")
from scipy.integrate import romberg
pi = np.pi
import matplotlib.patches as patches

tols = (1e-6 *(1 + 1j), 1e-6 *(1 + 1j), 1e-4 *(1 + 1j))
par_integ = (1e-8, 1e-8, 15)

def trace_rect(ax, z0, z1, **kwargs):
    ax.add_patch(
    patches.Rectangle(
        (z0.real, z0.imag),   # (x,y)
        (z1-z0).real,          # width
        (z1-z0).imag,          # height
        **kwargs
        ),
    )
    plt.pause(0.01)


def parametrization_line(z0, z1, loc):
    x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
    if loc == 0:
        a, b = 0, pi/2
        A, B = 2/pi*(x1-x0), z0
    elif loc == 1:
        a, b = pi/2, pi
        A, B = 1j*2/pi*(y1-y0), x1 + 1j*(2*y0 - y1)
    elif loc == 2 :
        a, b = pi, 3*pi/2
        A, B = -2/pi*(x1-x0), 3*x1-2*x0+1j*y1
    elif loc == 3:
        a, b = 3*pi/2, 2*pi
        A, B = -1j*2/pi*(y1-y0), x0+1j*(4*y1-3*y0)
    else:
        raise TypeError('loc must be an integer betwwen 0 and 3')
    return a, b, A, B


def poly_period(x, a, b):
    P = ((a - x)**3*(a**2 - 5*a*b + 3*a*x + 10*b**2 - 15*b*x + 6*x**2))/(a - b)**5
    dP = -(30*(a - x)**2*(b - x)**2)/(a - b)**5
    return P, dP

def get_integrands(a, b, A, B, x, func):
    P, dP = poly_period(x, a, b)
    Pt, dPt = (b - a)*P + a, (b - a)*dP
    PT = A*Pt + B
    fper = A*func(PT)*dPt
    return fper, PT*fper, PT**2*fper

def cut_in_four(z0, z1, ratio_re = 0.5, ratio_im = 0.5):
    x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
    xm, ym = x0 + ratio_re*(x1 - x0), y0 + ratio_im*(y1- y0)
    zm = xm +1j*ym
    zb = xm +1j*y0
    zt = xm +1j*y1
    zl = x0 +1j*ym
    zr = x1 +1j*ym
    return (z0, zm), (zb, zr), (zm, z1), (zl, zt)

def compute_integral(func, z0, z1, loc, par_integ = par_integ):
    a, b, A, B = parametrization_line(z0, z1, loc)
    Ire, Iim = np.zeros(3), np.zeros(3)
    for ii in range(3):
        fre = lambda x: get_integrands(a, b, A, B, x, func)[ii].real
        fim = lambda x: get_integrands(a, b, A, B, x, func)[ii].imag
        tol, rtol, divmax = par_integ
        Ire[ii] = romberg(fre, a, b, vec_func = True, tol=tol, rtol=rtol, divmax=divmax)
        Iim[ii] = romberg(fim, a, b, vec_func = True, tol=tol, rtol=rtol, divmax=divmax)
    return ( Ire + 1j*Iim )/(2*1j*pi)

def ispole(func, z0, z1, tols = tols, par_integ = par_integ):
    I = 0
    for loc in range(4):
        I += compute_integral(func, z0, z1, loc, par_integ = par_integ)
    R10 = I[1]/I[0]
    R21 = I[2]/I[1]
    r0, r1, r2 = False, False, False
    pole, residue = [], []
    if abs(I[0].real) < tols[0].real and abs(I[0].imag) < tols[0].imag and abs(I[1].real) < tols[1].real and abs(I[1].imag) < tols[1].imag:
        message ='No poles'
        r0 = True
    elif abs((R21-R10).real/R10.real) < tols[2].real and (abs((R21-R10).imag/R10.imag) < tols[2].real):
        message ='One pole'
        pole, residue = R10, I[0]
        r1 = True
    else:
        message='Several poles'
        r2 = True
    return pole, residue, r0, r1, r2, message


def pole_hunt(func, z0, z1, tols = tols, par_integ = par_integ, poles = [], residues = [], nb_cuts = 0):
    print('nb_cuts = ', nb_cuts)
    trace_rect(plt.gca(), z0, z1, fill = False, linewidth=3, edgecolor="#ff884d")
    pole, residue, r0, r1, r2, message = ispole(func, z0, z1, tols = tols,  par_integ = par_integ)
    print(message)
    trace_rect(plt.gca(), z0, z1, fill = False, linewidth=3, edgecolor="#000000")
    if r2:
        print('Cutting in four')
        nb_cuts += 1
        inv_golden_number =  1/(0.5*(1+np.sqrt(5)))
        Z = cut_in_four(z0, z1, ratio_re =inv_golden_number , ratio_im =inv_golden_number)
        for z in Z:
            _ , _, nb_cuts = pole_hunt(func, z[0], z[1], tols = tols,
             par_integ = par_integ, poles = poles, residues= residues, nb_cuts = nb_cuts)
    else:
        if r1:
            print('Found a new pole {}, with residue {}'.format(pole, residue))
            poles.append(pole)
            residues.append(residue)
            plt.gca().plot(np.real(pole), np.imag(pole), 'o', color = "#ff884d")
    # trace_rect(plt.gca(), z0, z1, fill = False, linewidth=3, edgecolor="#000000")
    poles, residues = np.array(poles), np.array(residues)
    isort = poles.argsort()
    return poles[isort], residues[isort], nb_cuts

if __name__ == '__main__':
    print("This is the tetrachotomy module")
