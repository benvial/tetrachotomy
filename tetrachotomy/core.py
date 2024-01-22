"""
Pole and zeros searching in the complex plane with a contour integral technique
"""

__all__ = ["get_options", "find_poles", "find_zeros"]

import copy
from collections import namedtuple
from math import factorial, pi

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from scipy.optimize import approx_fprime, minimize

from . import backend as bk
from . import logger

_epsilon = bk.sqrt(bk.finfo(float).eps)


def stencil5(f, x, h, epsilon):
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * epsilon)


def approx_fprime5(xk, f, epsilon=_epsilon):
    out = bk.zeros_like(xk)
    N = xk.shape[0]
    for i in range(N):
        h = bk.zeros(N)
        h[i] = epsilon
        q = stencil5(f, xk, h, epsilon)
        out[i] = q
    return out


def _tol_real2complex(tol):
    if bk.imag(tol) == 0:
        return tol * (1 + 1j)
    return tol


_ianim = 0


def _pause(save=False):
    global _ianim
    plt.pause(1e-4)
    if save:
        name = f"animation_tmp_{str(_ianim).zfill(4)}.png"
        plt.savefig(name)
        _ianim += 1


plot_options = dict(
    rectangle=False,
    circle=False,
    poles=False,
    remove=True,
    alpha=0.05,
    color="#ff884d",
    marker=dict(value="o", size=6),
)
tol_moments = {"0": 1e-5, "1": 1e-5, "ratio": 1e-4}
tolerances = dict(
    moments=tol_moments,
    integration=dict(tol=1e-7, divmax=13),
)
ratios = dict(real=0.5, imag=0.5)

ref_tol = dict(pole=1e-12, residue=1e-12)

refine_options = dict(
    nref_max=0,
    method="nelder-mead",
    circ_ratio=0.9,
    integration=dict(tol=1e-12, divmax=15),
    tolerances=ref_tol,
)

options = dict(
    tolerances=tolerances,
    plot=plot_options,
    ratios=ratios,
    ncuts_max=50,
    refine=refine_options,
)


def romberg(f, a, b, divmax=10, tol_re=1e-6, tol_im=1e-6, size=1):
    r = bk.zeros((divmax + 2, divmax + 2, size), dtype=bk.complex128)
    h = b - a
    r[0, 0] = 0.5 * h * (f(a) + f(b))
    neval = 2

    powerOf2 = 1
    for i in range(1, divmax + 2):
        h /= 2
        sigma = 0.0
        powerOf2 *= 2
        for k in range(1, powerOf2, 2):
            sigma += f(a + k * h)
            neval += 1
        r[i, 0] = 0.5 * r[i - 1, 0] + sigma * h

        powerOf4 = 1
        for j in range(1, i + 1):
            powerOf4 *= 4
            r[i, j] = r[i, j - 1] + (r[i, j - 1] - r[i - 1, j - 1]) / (powerOf4 - 1)

        current_tol_re = bk.abs(bk.real(r[i, j] - r[i - 1, j - 1]))
        current_tol_im = bk.abs(bk.imag(r[i, j] - r[i - 1, j - 1]))

        logger.debug(f"Integration iteration {i}")
        logger.debug(f"Integration current error (real) = {current_tol_re}")
        logger.debug(f"Integration current error (imag) = {current_tol_im}")
        if i > 0 and (
            bk.all(current_tol_re < tol_re) and bk.all(current_tol_im < tol_im)
        ):
            return r[i, i], neval
    logger.warning(f"Accuracy warning: divmax={divmax} exceeded")
    logger.warning(f"Integration current error (real) = {current_tol_re}")
    logger.warning(f"Integration current error (imag) = {current_tol_im}")
    return r[i, i], neval


def plot_rectangle(ax, z0, z1, **kwargs):
    patch = patches.Rectangle(
        (z0.real, z0.imag),  # (x,y)
        (z1 - z0).real,  # width
        (z1 - z0).imag,  # height
        **kwargs,
    )
    ax.add_patch(patch)
    _pause()
    return patch


def plot_circle(ax, zp, r, **kwargs):
    patch = patches.Circle((zp.real, zp.imag), r, **kwargs)  # center (x,y)  # radius
    ax.add_patch(patch)
    _pause()
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
        raise TypeError("loc must be an integer betwwen 0 and 3")
    return a, b, A, B


def poly_period(x, a, b):
    P = (
        (a - x) ** 3
        * (a**2 - 5 * a * b + 3 * a * x + 10 * b**2 - 15 * b * x + 6 * x**2)
    ) / (a - b) ** 5
    dP = -(30 * (a - x) ** 2 * (b - x) ** 2) / (a - b) ** 5
    return P, dP


def get_integrands(a, b, A, B, x, func, order):
    P, dP = poly_period(x, a, b)
    Pt, dPt = (b - a) * P + a, (b - a) * dP
    PT = A * Pt + B
    fper = A * func(PT) * dPt
    return bk.array(
        bk.stack(
            [PT ** (order - 1) * fper, PT ** (order) * fper, PT ** (order + 1) * fper]
        ),
        dtype=bk.complex128,
    )


def cut_in_four(z0, z1, ratio_re=0.5, ratio_im=0.5):
    logger.info("Cutting in four")
    x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
    xm, ym = x0 + ratio_re * (x1 - x0), y0 + ratio_im * (y1 - y0)
    zm = xm + 1j * ym
    zb = xm + 1j * y0
    zt = xm + 1j * y1
    zl = x0 + 1j * ym
    zr = x1 + 1j * ym
    return (z0, zm), (zb, zr), (zm, z1), (zl, zt)


def compute_integral(func, z0, z1, loc, par_integ, order):
    tol_re, tol_im, divmax = par_integ

    a, b, A, B = parametrization_line(z0, z1, loc)

    def f(x):
        return get_integrands(a, b, A, B, x, func, order)

    I, neval = romberg(f, a, b, tol_re=tol_re, tol_im=tol_im, divmax=divmax, size=3)
    return I / (2 * 1j * pi), neval


def compute_R10(I, order):
    R10 = factorial(order - 1) / factorial(order) * I[1] / I[0]
    return R10


def compute_R21(I, order):
    R21 = (2 * factorial(order)) / factorial(order + 1) * I[2] / I[1]
    return R21


def ispole(func, z0, z1, tols, par_integ, order):
    nfeval = 0
    I = 0
    for loc in range(4):
        val, neval = compute_integral(func, z0, z1, loc, par_integ, order)
        I += val
        nfeval += neval

    R10 = compute_R10(I, order)
    R21 = compute_R21(I, order)

    logger.info(f"moment 0  = {I[0]}")
    logger.info(f"moment 1  = {I[1]}")
    logger.info(f"moment 2  = {I[2]}")
    logger.info(f"ratio 1/0 = {R10}")
    logger.info(f"ratio 2/1 = {R21}")

    r0, r1, r2 = False, False, False
    pole, residue = [], []

    # mid = (z0 + z1) / 2
    # val_mid = func(mid)
    # peri = 2 * (z1.real - z0.real) + 2 * (z1.imag - z0.imag)
    # Iapprox = [1 / (2 * 1j * pi) * mid**n * val_mid * peri for n in range(2)]
    # cond0_re = abs(1 - I[0].real / Iapprox[0].real)
    # cond0_im = abs(1 - I[0].imag / Iapprox[0].imag)
    # cond1_re = abs(1 - I[1].real / Iapprox[1].real)
    # cond1_im = abs(1 - I[1].imag / Iapprox[1].imag)

    cond0_re = abs(I[0].real)
    cond0_im = abs(I[0].imag)
    cond1_re = abs(I[1].real)
    cond1_im = abs(I[1].imag)
    if (
        cond0_re < tols[0].real
        and cond0_im < tols[0].imag
        and cond1_re < tols[1].real
        and cond1_im < tols[1].imag
    ):
        # if (
        #     abs(I[0].real) < tols[0].real
        #     and abs(I[0].imag) < tols[0].imag
        #     and abs(I[1].real) < tols[1].real
        #     and abs(I[1].imag) < tols[1].imag
        # ):
        message = "No poles"
        r0 = True
    elif abs((R21 - R10).real / R10.real) < tols[2].real and (
        abs((R21 - R10).imag / R10.imag) < tols[2].imag
    ):
        # elif check_ratios < tols[2].real:
        message = "One pole"
        pole, residue = R10, I[0]
        r1 = True
    else:
        message = "Several poles"
        r2 = True
    return pole, residue, r0, r1, r2, message, nfeval


def get_integrands_circ(t, zp, r, func, order):
    Z = zp + r * bk.exp(1j * t)
    K = 1j * r * bk.exp(1j * t)
    fr = K * func(Z)

    return bk.array([Z ** (order - 1) * fr, Z ** (order) * fr], dtype=bk.complex128)


def compute_integral_circ(func, zp, r, par_integ, order):
    tol_re, tol_im, divmax = par_integ

    def f(t):
        return get_integrands_circ(t, zp, r, func, order)

    I, neval = romberg(
        f, 0, 2 * pi, tol_re=tol_re, tol_im=tol_im, divmax=divmax, size=2
    )
    return I / (2 * 1j * pi), neval


def min_dist(z0, z1, zp):
    x0, x1, y0, y1, xp, yp = z0.real, z1.real, z0.imag, z1.imag, zp.real, zp.imag
    return bk.min(bk.array([xp - x0, yp - y0, x1 - xp, y1 - yp]))


def refine_pole_optimize(func, zp, zr, options=options):
    x0 = [zp.real, zp.imag]
    tolerances = options["refine"]["tolerances"]
    tol_pol = tolerances["pole"]
    method = options["refine"]["method"].lower()

    def func_min(x):
        z = x[0] + x[1] * 1j
        return bk.abs(1 / func(z))

    if method in ("nelder-mead", "powell", "cobyla"):
        jac = None
    else:
        jac = lambda x: approx_fprime(x, func_min, 1.5e-8)

    opt = minimize(
        func_min,
        x0,
        method=method,
        tol=tol_pol.real,
        jac=jac,
        # options={"disp": True, "xtol": 1e-12, "eps": 1e-12},
    )
    pole_ref = opt.x[0] + opt.x[1] * 1j
    residue_ref = zr
    nfeval = opt.nfev

    logger.info(f"Refined pole {pole_ref}")
    logger.info(f"With residue {residue_ref}")

    return pole_ref, residue_ref, nfeval


def refine_pole_cauchy(func, z0, z1, zp, zr, order, options=options):
    tolerances = options["refine"]["tolerances"]
    tol_pol = _tol_real2complex(tolerances["pole"])
    tol_res = _tol_real2complex(tolerances["residue"])
    par_integ_ = options["refine"]["integration"]
    par_integ_tol = _tol_real2complex(par_integ_["tol"])
    par_integ = (
        par_integ_tol.real,
        par_integ_tol.imag,
        par_integ_["divmax"],
    )
    nref_max = options["refine"]["nref_max"]
    ratio_circ = options["refine"]["circ_ratio"]
    plot_options = options["plot"]
    color = plot_options["color"]

    r = min_dist(z0, z1, zp)
    conv_pol_re, conv_pol_im, conv_res_re, conv_res_im = 1, 1, 1, 1
    nref = 0
    nfeval = 0
    while nref < nref_max and (
        conv_pol_re > tol_pol.real
        or conv_pol_im > tol_pol.imag
        or conv_res_re > tol_res.real
        or conv_res_im > tol_res.imag
    ):
        if plot_options["circle"]:
            circ_b = plot_circle(
                plt.gca(), zp, r, fill=False, linewidth=0.1, edgecolor=color
            )
            circ = plot_circle(
                plt.gca(),
                zp,
                r,
                fill=True,
                linewidth=0,
                facecolor=color,
                alpha=plot_options["alpha"],
            )

        logger.info(f"Refinement iteration {nref}")
        I, neval = compute_integral_circ(func, zp, r, par_integ, order)
        nfeval += neval
        R10 = compute_R10(I, order)
        pole_ref, residue_ref = R10, I[0]
        conv_pol_re = bk.abs(bk.real(zp - pole_ref) / bk.real(zp))
        conv_pol_im = bk.abs(bk.imag(zp - pole_ref) / bk.imag(zp))
        conv_res_re = bk.abs(bk.real(zr - residue_ref) / bk.real(zr))
        conv_res_im = bk.abs(bk.imag(zr - residue_ref) / bk.imag(zr))

        logger.info(f"Refined pole {pole_ref}")
        logger.info(f"With residue {residue_ref}")
        logger.debug(f"Convergence pole (real) = {conv_pol_re},(imag) = {conv_pol_im}")
        logger.debug(
            f"Convergence residue (real) = {conv_res_re},(imag) = {conv_res_im}"
        )
        r *= ratio_circ
        zp = pole_ref
        zr = residue_ref
        if plot_options["circle"] and plot_options["remove"]:
            circ.remove()
            circ_b.remove()
        nref += 1
        if nref >= nref_max:
            logger.warning("Accuracy warning, nref_max exceeded")

    return pole_ref, residue_ref, nfeval


def refine_pole(func, z0, z1, zp, zr, order, options=options):
    method = options["refine"]["method"].lower()

    MINIMIZE_METHODS = [
        "cauchy",
        "nelder-mead",
        "powell",
        "cg",
        "bfgs",
        "newton-cg",
        "l-bfgs-b",
        "tnc",
        "cobyla",
        "slsqp",
        "trust-constr",
        "dogleg",
        "trust-ncg",
        "trust-exact",
        "trust-krylov",
    ]
    if method not in MINIMIZE_METHODS:
        raise ValueError(
            f"Unknown method {method}, should be one of {MINIMIZE_METHODS}"
        )
    logger.info(f"Refining pole using {method} method")
    if method == "cauchy":
        return refine_pole_cauchy(func, z0, z1, zp, zr, order, options=options)
    else:
        return refine_pole_optimize(func, zp, zr, options=options)


def pole_hunt(
    func,
    z0,
    z1,
    options=options,
    order=1,
    zeros=False,
    init=None,
):
    if init is None:
        init = namedtuple("Result", ["poles", "residues", "ncuts", "nfeval"])
        init.poles = []
        init.residues = []
        init.ncuts = 0
        init.nfeval = 0
    res = copy.copy(init)
    # poles = res.poles
    # residues = res.residues
    # nfeval = res.nfeval
    ratios = options["ratios"]
    ratio_re = ratios["real"]
    ratio_im = ratios["imag"]
    ratio_circ = options["refine"]["circ_ratio"]
    plot_options = options["plot"]
    tolerances = options["tolerances"]
    tols_ = tolerances["moments"]
    tols = (
        _tol_real2complex(tols_["0"]),
        _tol_real2complex(tols_["1"]),
        _tol_real2complex(tols_["ratio"]),
    )
    par_integ_ = tolerances["integration"]
    par_integ_tol = _tol_real2complex(par_integ_["tol"])
    par_integ = par_integ_tol.real, par_integ_tol.imag, par_integ_["divmax"]
    ncuts_max = options["ncuts_max"]

    nref_max = options["refine"]["nref_max"]

    plot_options = options["plot"]
    color = plot_options["color"]

    if res.ncuts > ncuts_max:
        return res
    if zeros:
        func1 = lambda z: 1 / func(z)
    else:
        func1 = func
    if plot_options["rectangle"]:
        # if ncuts == 0:
        #     plt.gca().set_xlim(z0.real, z1.real)
        #     plt.gca().set_ylim(z0.imag, z1.imag)
        patch1 = plot_rectangle(
            plt.gca(),
            z0,
            z1,
            fill=True,
            linewidth=1,
            facecolor=color,
            alpha=plot_options["alpha"],
        )
        patch2 = plot_rectangle(
            plt.gca(),
            z0,
            z1,
            fill=False,
            linewidth=3,
            edgecolor=color,
            alpha=1,
        )
    logger.info(f"Number of subdivisions = {res.ncuts}")
    pole, residue, r0, r1, r2, message, nfeval_ = ispole(
        func1,
        z0,
        z1,
        tols=tols,
        par_integ=par_integ,
        order=order,
    )
    res.nfeval += nfeval_
    logger.info(message)
    if r2:
        res.ncuts += 1
        Z = cut_in_four(z0, z1, ratio_re=ratio_re, ratio_im=ratio_im)
        for z in Z:
            res = pole_hunt(
                func1,
                z[0],
                z[1],
                options=options,
                order=order,
                init=res,
            )
    else:
        if r1:
            logger.info(f"Found a new pole {pole}")
            logger.info(f"Residue {residue}")
            if nref_max > 0:
                pole, residue, nfeval_ref = refine_pole(
                    func1,
                    z0,
                    z1,
                    pole,
                    residue,
                    order,
                    options=options,
                )
                res.nfeval += nfeval_ref
            res.poles.append(pole)
            res.residues.append(residue)
            if plot_options["poles"]:
                plt.gca().plot(
                    bk.real(pole),
                    bk.imag(pole),
                    plot_options["marker"]["value"],
                    ms=plot_options["marker"]["size"],
                    c=color,
                    linewidth=1,
                )
                _pause()
    if plot_options["rectangle"] and plot_options["remove"]:
        patch1.remove()
        patch2.remove()

    _pause()
    return res


def final_message(poles, residues, ncuts, nfeval):
    n = len(poles)
    message_result = f"Found {n} pole{'s' if n>1 else ''}"
    message_result += f" with {ncuts} contour subdivision{'s' if ncuts>1 else ''} and {nfeval} function evaluations"
    message_result += "\n" + " " * 10 + "poles" + " " * 50 + "residues"
    message_result += "\n" + "-" * 90

    for i in range(len(poles)):
        message_result += f"\n{poles[i]}    {residues[i]}"

    return message_result


def get_options():
    return copy.deepcopy(options)


def _wrapper(*args, **kwargs):
    res = pole_hunt(*args, **kwargs)
    logger.info(final_message(res.poles, res.residues, res.ncuts, res.nfeval))
    poles, residues = bk.array(res.poles, dtype=bk.complex128), bk.array(
        res.residues, dtype=bk.complex128
    )
    isort = poles.real.argsort()
    res.poles = poles[isort]
    res.residues = residues[isort]
    return res


def find_poles(
    func,
    z0,
    z1,
    order=1,
    options=options,
):
    return _wrapper(func, z0, z1, options=options, order=order, zeros=False, init=None)


def find_zeros(
    func,
    z0,
    z1,
    order=1,
    options=options,
):
    return _wrapper(func, z0, z1, options=options, order=order, zeros=True, init=None)
