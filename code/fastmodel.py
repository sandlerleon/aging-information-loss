# -*- coding: utf-8 -*-
"""Numba-compiled batch integrator for the same equations as model.rhs (identical RK4 scheme and step), used for the large
ensembles. model.simulate dispatches here for final-state runs; the NumPy implementation remains the reference and the two are
compared in verify.py (check V2b)."""
import numpy as np
from numba import njit, prange

KEYS = ["r0", "rg", "rm", "csig", "aI", "aF", "aM", "au", "lam", "kM", "phi", "ku", "nH", "pm", "gF", "cF", "kappa", "G0", "rf", "pf", "Fref"]
HCODE = {"hill": 0, "mm": 1, "exp": 2}
PCODE = {"power": 0, "threshold": 1, "saturating": 2}


def pack(P):
    return np.ascontiguousarray(np.stack([np.asarray(P[k], float) for k in KEYS], axis=1))


@njit(cache=True, inline="always")
def _rhs(I, M, F, uG, uM, uF, p, hc, pc, out):
    r0, rg, rm, csig, aI, aF, aM, au, lam, kM, phi, ku, nH, pm, gF, cF, kappa, G0, rf, pf, Fref = (
        p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18], p[19], p[20])
    G = G0 * (1.0 + kappa * aI * I + aF * F) / (1.0 + aM * M + au * uG)
    sig = 1.0 + kappa * csig * I
    d = 1.0 - M
    if d < 0.0:
        d = 0.0
    if pc == 0:
        pen = d ** pm
    elif pc == 1:
        x = d - 0.15
        if x < 0.0:
            x = 0.0
        pen = (x / 0.85) ** pm
    else:
        c = 0.25
        dp = d ** pm
        pen = dp * (1.0 + c) / (c + dp)
    ex = F - Fref
    if ex < 0.0:
        ex = 0.0
    fat = ex ** pf
    bracket = rg * G * sig + rm * pen + rf * fat
    dI = (r0 + kappa * bracket) * (1.0 - I)
    if hc == 0:
        un = uM ** nH
        H = un / (ku ** nH + un)
    elif hc == 1:
        H = uM / (ku + uM)
    else:
        H = 1.0 - np.exp(-uM / ku)
    dM = kM * ((1.0 - lam * I) * (1.0 - phi * H) - M)
    dF = gF * (0.4 + 0.4 * kappa * I - F) - cF * uF * F
    out[0] = dI
    out[1] = dM
    out[2] = dF


@njit(parallel=True, cache=True)
def _final(par, u, t_start, t_end, dt, hc, pc, I0, M0, F0, age0):
    n = par.shape[0]
    res = np.empty((n, 3))
    nsteps = int(round((t_end - age0) / dt))
    for i in prange(n):
        p = par[i]
        I, M, F = I0, M0, F0
        k1 = np.empty(3)
        k2 = np.empty(3)
        k3 = np.empty(3)
        k4 = np.empty(3)
        for k in range(nsteps):
            t = age0 + dt * k
            if t >= t_start - 1e-9:
                uG, uM, uF = u[i, 0], u[i, 1], u[i, 2]
            else:
                uG, uM, uF = 0.0, 0.0, 0.0
            _rhs(I, M, F, uG, uM, uF, p, hc, pc, k1)
            _rhs(I + 0.5 * dt * k1[0], M + 0.5 * dt * k1[1], F + 0.5 * dt * k1[2], uG, uM, uF, p, hc, pc, k2)
            _rhs(I + 0.5 * dt * k2[0], M + 0.5 * dt * k2[1], F + 0.5 * dt * k2[2], uG, uM, uF, p, hc, pc, k3)
            _rhs(I + dt * k3[0], M + dt * k3[1], F + dt * k3[2], uG, uM, uF, p, hc, pc, k4)
            I = I + dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            M = M + dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            F = F + dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        res[i, 0] = I
        res[i, 1] = M
        res[i, 2] = F
    return res


def final_state(P, u, t_start, t_end, dt, hform, pform, I0, M0, F0, age0):
    n = len(P["r0"])
    u = np.asarray(u, float)
    if u.ndim < 2:
        u = np.repeat(np.broadcast_to(u, (n,))[:, None], 3, axis=1)
    return _final(pack(P), np.ascontiguousarray(u), float(t_start), float(t_end), float(dt), HCODE[hform], PCODE[pform], I0, M0, F0, age0)
