


def complex_quadrature(func, a, b, **kwargs):
    def real_func(x):
        return sc.real(func(x))
    def imag_func(x):
        return sc.imag(func(x))
    real_integral = romberg(real_func, a, b, **kwargs)
    imag_integral = romberg(imag_func, a, b, **kwargs)
    return real_integral + 1j*imag_integral


a,b = 0,1
it = 0
t0 = time.time()
while it < 15:
    # print('-'*32)
    # print('Iteration ', it)
    nx = 2**it + 1
    xint = np.linspace(a, b, nx)
    if it == 0:

        intgr = fre(xint)
        # intgr = get_integrands(a, b, A, B, xint, func)[0]
    else:
        xint_new = xint[1::2]
        intgr_new = fre(xint_new)
        # intgr_new = get_integrands(a, b, A, B, xint_new, func)[0]
        intgr = np.zeros((1, nx), dtype = complex)
        intgr[:, 0::2] = intgr_old
        intgr[:, 1::2] = intgr_new
    intgr_old = np.copy(intgr)
    integral = sc.trapz(intgr, xint)
    print('integral = ', integral)
    it += 1
t1 = time.time() -t0
print('t1 = ', t1)

        # t0 = time.time()
        # II = np.zeros(3, dtype = complex)
        # for ii in range(3):
        #     f = lambda x: get_integrands(a, b, A, B, x, func)[ii]
        #     II[ii] = complex_quadrature(f, a, b)
        # tr2= time.time() -t0
        # print('integral Romberg 2 = ', II)
        # print('tr2 = ', tr2)

        # complex_quadrature(func, a, b, **kwargs)
