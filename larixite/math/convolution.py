#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generic discrete (one-dimensional) convolution

Description
-----------

 This is a manual (not optimized!) implementation of discrete one-dimensional
 (1D) convolution intended for spectroscopy analysis. The difference with
 commonly used methods is the possibility to adapt the width of the convolution
 kernel for each convolution point, e.g. change the FWHM of the Gaussian kernel
 as a function of the energy scale (x-axis).

Resources
---------

.. [WPconv] <http://en.wikipedia.org/wiki/Convolution#Discrete_convolution>
.. [Fisher] <http://homepages.inf.ed.ac.uk/rbf/HIPR2/convolve.htm>
.. [GP1202] <http://glowingpython.blogspot.fr/2012/02/convolution-with-numpy.html>

"""

import numpy as np

s2pi = np.sqrt(2 * np.pi)
tiny = 1.0e-15

DEBUG = 0


def gaussian(x, amplitude=1.0, center=0.0, sigma=1.0):
    """Return a 1-dimensional Gaussian function."""
    return (amplitude / (max(tiny, s2pi * sigma))) * np.exp(
        -((1.0 * x - center) ** 2) / max(tiny, (2 * sigma**2))
    )


def lorentzian(x, amplitude=1.0, center=0.0, sigma=1.0):
    """Return a 1-dimensional Lorentzian function."""
    return (amplitude / (1 + ((1.0 * x - center) / max(tiny, sigma)) ** 2)) / max(
        tiny, (np.pi * sigma)
    )


def polyfit(x, y, deg=1, reverse=False):
    """
    simple emulation of deprecated numpy.polyfit,
    including its ordering of coefficients
    """
    pfit = np.polynomial.Polynomial.fit(x, y, deg=int(deg))
    coefs = pfit.convert().coef
    if reverse:
        coefs = list(reversed(coefs))
    return list(coefs)


def get_ene_index(ene, cen, hwhm):
    """returns the min/max indexes for array ene at (cen-hwhm) and (cen+hwhm)
    very similar to index_of in larch
    """
    try:
        if (cen - hwhm) <= min(ene):
            ene_imin = 0
        else:
            ene_imin = max(np.where(ene < (cen - hwhm))[0])
        if (cen + hwhm) >= max(ene):
            ene_imax = len(e) - 1
        else:
            ene_imax = min(np.where(ene > (cen + hwhm))[0])
        return ene_imin, ene_imax
    except Exception:
        print("index not found for {0} +/- {1}".format(cen, hwhm))
        return None, None


def lin_gamma(ene, gamma_hole=0.5, linbroad=None):
    """returns constant or linear energy-dependent broadening

    Parameters
    ----------
    ene : energy array in eV
    gamma_hole : initial full width at half maximum in eV
    linbroad   : list of 3-elements giving
                    'final full width at half maximum'
                    'start energy point of the linear increase'
                    'end energy point of the linear increase'
    """
    w = np.ones_like(ene)
    if linbroad is None:
        return w * gamma_hole
    else:
        try:
            fwhm2 = linbroad[0]
            e1 = linbroad[1]
            e2 = linbroad[2]
        except:
            raise ValueError("wrong format for linbroad")
        for en, ee in enumerate(ene):
            if ee < e1:
                w[en] *= gamma_hole
            elif ee <= e2:
                wlin = gamma_hole + (ee - e1) * (fwhm2 - gamma_hole) / (e2 - e1)
                w[en] *= wlin
            elif ee >= e2:
                w[en] *= fwhm2
        return w


def atan_gamma(e, e_cut=0, e_cent=30, e_larg=30, gamma_hole=0.5, gamma_max=15):
    """Arctangent energy dependent broadening as implemented in FDMNES (described in the "Convolution" section of the manual)"""

    f = (e - e_cut) / e_cent
    a = np.pi * gamma_max * (f - 1 / f**2) / (3 * e_larg)
    gammas = gamma_hole + gamma_max * (0.5 + 1 / np.pi * np.arctan(a))

    # Set gamma to the gamma_hole below the cutoff energy.
    mask = np.where(e < e_cut)
    gammas[mask] = gamma_hole

    return gammas


def conv_fast(x, y, gammas, e_cut=None, kernel="gaussian", num=501, step=0.1):
    """A significantly faster version of the `conv` function

    Parameters
    ----------
    x : x-axis (energy)
    y : f(x) to convolve with g(x) kernel, y(energy)
    kernel : convolution kernel, g(x)
             'gaussian'
             'lorentzian'
    gammas : the full width half maximum in eV for the kernel
            broadening. It is an array of size 'e' with constants or
            an energy-dependent values determined by a function as
            'lin_gamma()' or 'atan_gamma()'
    """
    assert e_cut is not None, "starting energy for the convolution not given"
    assert "lor" or "gauss" in kernel.lower(), (
        "the kernel should be either Lorentzian or Gaussian"
    )

    # Extend the X-axis array.
    start = x[-1]
    stop = start + num * step
    x_ext = np.append(x, np.arange(start + step, stop, step))

    ids = x_ext > e_cut

    x1 = x_ext[ids]
    x1[0] = e_cut

    # Extend the intensity array by coping the last value `num - 1` times.
    y = y[x > e_cut]
    y = np.append(y, np.ones(num - 1) * y[-1])

    y_conv = np.zeros_like(x)
    for i, (xi, gamma) in enumerate(zip(x, gammas)):
        gamma = gamma / 2.0
        if "gauss" in kernel.lower():
            ky = gaussian(x1, center=xi, sigma=gamma)
        elif "lor" in kernel.lower():
            ky = lorentzian(x1, center=xi, sigma=gamma)
        y_conv[i] = np.sum(ky * y) / np.pi

    return y_conv


def conv(x, y, gammas, e_cut=None, kernel="gaussian"):
    """linear broadening

    Parameters
    ----------
    x : x-axis (energy)
    y : f(x) to convolve with g(x) kernel, y(energy)
    kernel : convolution kernel, g(x)
             'gaussian'
             'lorentzian'
    gammas : the full width half maximum in eV for the kernel
            broadening. It is an array of size 'e' with constants or
            an energy-dependent values determined by a function as
            'lin_gamma()' or 'atan_gamma()'
    """
    assert e_cut is not None, "starting energy for the convolution not given"

    f = y[:] * 1.0
    z = np.zeros_like(f)
    # ief = index_nearest(x, e_cut)
    ief = np.argmin(np.abs(x - e_cut))
    f[0:ief] *= 0

    if x.shape != gammas.shape:
        print("Error: 'gammas' array does not have the same shape of 'x'")
        return 0

    # linear fit upper part of the spectrum to avoid border effects
    # polyfit => pf
    lpf = int(len(x) / 2.0)
    cpf = polyfit(x[-lpf:], f[-lpf:], 1, reverse=False)
    fpf = np.polynomial.Polynomial(cpf)

    # extend upper energy border to 3*fhwm_e[-1]
    xstep = x[-1] - x[-2]
    xup = np.append(x, np.arange(x[-1] + xstep, x[-1] + 3 * gammas[-1], xstep))

    for n in range(len(f)):
        # from now on I change e with eup
        eimin, eimax = get_ene_index(xup, xup[n], 1.5 * gammas[n])
        if (eimin is None) or (eimax is None):
            if DEBUG:
                raise IndexError("e[{0}]".format(n))
        if len(range(eimin, eimax)) % 2 == 0:
            kx = xup[eimin : eimax + 1]  # odd range centered at the convolution point
        else:
            kx = xup[eimin:eimax]
        ### kernel ###
        hwhm = gammas[n] / 2.0
        if "gauss" in kernel.lower():
            ky = gaussian(kx, center=xup[n], sigma=hwhm)
        elif "lor" in kernel.lower():
            ky = lorentzian(kx, center=xup[n], sigma=hwhm)
        else:
            raise ValueError("convolution kernel '{0}' not implemented".format(kernel))
        ky = ky / ky.sum()  # normalize
        zn = 0
        lk = len(kx)
        for mf, mg in zip(range(-int(lk / 2), int(lk / 2) + 1), range(lk)):
            if ((n + mf) >= 0) and ((n + mf) < len(f)):
                zn += f[n + mf] * ky[mg]
            elif (n + mf) >= 0:
                zn += fpf(xup[n + mf]) * ky[mg]
        z[n] = zn
    return z
