import time
import numpy as np
import scipy as sc
from scipy.integrate import romberg
pi = np.pi


# rel_tol_I0 = 1e-4 *(1 + 1j),
# rel_tol_I1 = 1e-8 *(1 + 1j),
# rel_tol_R = 1e-4 *(1 + 1j),
# rel_tol_cv = 1e-3,


def poly_period(a, b, x):
    P = ((a - x)**3*(a**2 - 5*a*b + 3*a*x + 10*b**2 - 15*b*x + 6*x**2))/(a - b)**5
    dP = -(30*(a - x)**2*(b - x)**2)/(a - b)**5
    return P, dP

def get_integrands(a, b, A, B, x, func):
    P, dP = poly_period(a, b, x)
    Pt, dPt =(b - a)*P+a, (b - a)*dP
    PT = A*Pt + B;
    fper = A*func(PT)
    return fper, Pt*fper, Pt**2*fper

def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return sc.real(func(x))
    def imag_func(x):
        return sc.imag(func(x))
    real_integral = romberg(real_func, a, b, **kwargs)
    imag_integral = romberg(imag_func, a, b, **kwargs)
    return real_integral + 1j*imag_integral


class PoleHunt:
    """Tetrachotomy method to find poles in the complex plane.


    """
    def __init__(self,
        z0 = -1 - 1j, z1 = 1 + 1j,
        tols = (1e-4 *(1 + 1j), 1e-8 *(1 + 1j), 1e-4 *(1 + 1j), 1e-3)
        ):
        self.z0, self.z1 = z0, z1
        self. tols = tols


    def parametrization_line(self, loc):
        x0, x1, y0, y1 = self.z0.real, self.z1.real, self.z0.imag, self.z1.imag
        if loc == 0:
            a, b = 0, pi/2
            A, B = 2/pi*(x1-x0), self.z0
        elif loc == 1:
            a, b = pi/2, pi
            A, B = 1j*2/pi*(y1-y0), x1 + 1j*(2*y0 - y1)
        elif loc == 2 :
            a, b = pi, 3*pi/2
            A, B = -2/pi*(x1-x0), 3*x1-2*x0+i*y1
        elif loc == 3:
            a, b = pi, 3*pi/2
            A, B = -1j*2/pi*(y1-y0), x0+1j*(4*y1-3*y0)
        else:
            raise TypeError('loc must be an integer betwwen 0 and 3')
        return a, b, A, B

    def compute_integral(self, func, loc):
        a, b, A, B = self.parametrization_line(loc)
        it = 0
        t0 = time.time()
        while it < 9:
            print('-'*32)
            print('Iteration ', it)
            nx = 2**it + 1
            xint = np.linspace(a, b, nx)
            if it == 0:
                intgr = get_integrands(a, b, A, B, xint, func)[0]
            else:
                xint_new = xint[1::2]
                intgr_new = get_integrands(a, b, A, B, xint_new, func)[0]
                intgr = np.zeros((1, nx), dtype = complex)
                intgr[:, 0::2] = intgr_old
                intgr[:, 1::2] = intgr_new
            intgr_old = np.copy(intgr)
            integral = sc.trapz(intgr, xint)
            print('integral = ', integral)
            it += 1
        t1 = time.time() -t0
        print('t1 = ', t1)


        print('#'*55)
        print('Romberg')

        t0 = time.time()
        Ire, Iim = np.zeros(3), np.zeros(3)
        for ii in range(1):
            fre = lambda x: get_integrands(a, b, A, B, x, func)[ii].imag
            # fim = lambda x: get_integrands(a, b, A, B, x, func)[ii].imag
            Ire[ii] = romberg(fre, a, b, show =False,  vec_func = True)
            # Iim[ii] = romberg(fim, a, b, show =False,  vec_func = True)
        tr= time.time() -t0
        I = Ire+1j*Iim
        print('integral Romberg = ', I)
        print('tr = ', tr)
        return integral, I

        # t0 = time.time()
        # II = np.zeros(3, dtype = complex)
        # for ii in range(3):
        #     f = lambda x: get_integrands(a, b, A, B, x, func)[ii]
        #     II[ii] = complex_quadrature(f, a, b)
        # tr2= time.time() -t0
        # print('integral Romberg 2 = ', II)
        # print('tr2 = ', tr2)

        # complex_quadrature(func, a, b, **kwargs)






if __name__ == '__main__':
    print("This is the tetrachotomy module")
