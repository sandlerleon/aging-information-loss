# -*- coding: utf-8 -*-
"""Glucose-muscle-information-loss model of metabolic aging under a glucose-lowering control input.

State: I (information loss, 0..1), M (muscle reserve, young value 1), F (energy-storage compartment, adiposity).
Algebraic: G (circulating glucose, quasi-steady state), sigma (susceptibility), SIM = G sigma / M (derived readout).

    G      = G0 (1 + kappa a_I I + a_F F) / (1 + a_M M + a_u u)                      (1)
    sigma  = 1 + kappa c_sig I                                                       (2)
    dI/dt  = ( r0 + kappa [ r_g G sigma + r_m P(1-M) + r_f Phi(F) ] ) (1 - I)         (3)
    dM/dt  = k_M [ (1 - lam I)(1 - phi H(u)) - M ]                                   (4)
    dF/dt  = g_F (F_set(I) - F) - c_F u F,   F_set = 0.4 + 0.4 kappa I               (5)

P(d) = d^p_m (default); H(u) = u^nH / (ku^nH + u^nH) (default); Phi(F) = max(F - F_ref, 0)^p_f (r_f = 0 by default).
Everything is vectorised over a batch of parameter sets with a fixed-step RK4 integrator (checked against LSODA).
"""
import numpy as np

NOMINAL = dict(r0=0.0020, rg=0.0130, rm=0.0400, csig=0.90, aI=0.80, aF=0.45, aM=1.00, au=0.55, lam=0.28, kM=0.45,
               phi=0.50, ku=1.00, nH=2.0, pm=2.0, gF=0.055, cF=0.16, kappa=1.0, G0=1.0,
               rf=0.0, pf=1.0, Fref=0.4)
# parameters sampled in the robustness ensemble (scale parameters) and the two functional-form exponents
SCALE = ["r0", "rg", "rm", "csig", "aI", "aF", "aM", "au", "lam", "kM", "phi", "ku", "gF", "cF"]
EXPON = ["nH", "pm"]
AGE0, AGE_START, AGE_END = 25.0, 45.0, 75.0
I0, M0, F0 = 0.10, 1.0, 0.40


def make(n=1, **over):
    P = {k: np.full(n, float(v)) for k, v in NOMINAL.items()}
    for k, v in over.items():
        P[k] = np.broadcast_to(np.asarray(v, float), (n,)).copy()
    return P


def tile(P, reps):
    return {k: np.repeat(v, reps) for k, v in P.items()}


def _H(u, P, form):
    if form == "hill":
        un = u ** P["nH"]
        return un / (P["ku"] ** P["nH"] + un)
    if form == "mm":                       # Michaelis-Menten (no threshold)
        return u / (P["ku"] + u)
    if form == "exp":
        return 1.0 - np.exp(-u / P["ku"])
    raise ValueError(form)


def _pen(d, P, form):
    d = np.maximum(d, 0.0)
    if form == "power":
        return d ** P["pm"]
    if form == "threshold":                # no cost until 15 % of reserve is lost
        return (np.maximum(d - 0.15, 0.0) / 0.85) ** P["pm"]
    if form == "saturating":
        c = 0.25
        return d ** P["pm"] * (1 + c) / (c + d ** P["pm"])
    raise ValueError(form)


def algebra(I, M, F, u, P):
    G = P["G0"] * (1 + P["kappa"] * P["aI"] * I + P["aF"] * F) / (1 + P["aM"] * M + P["au"] * u)
    sig = 1 + P["kappa"] * P["csig"] * I
    return G, sig


def rhs(y, u, P, hform="hill", pform="power"):
    """u: scalar / (n,) dose acting on all channels, or (n,3) = separate doses (glucose, muscle, storage) for the
    exact pathway decomposition of dI(T)/du."""
    I, M, F = y[:, 0], y[:, 1], y[:, 2]
    u = np.asarray(u, float)
    uG, uM, uF = (u[:, 0], u[:, 1], u[:, 2]) if u.ndim == 2 else (u, u, u)
    G, sig = algebra(I, M, F, uG, P)
    bracket = P["rg"] * G * sig + P["rm"] * _pen(1 - M, P, pform) + P["rf"] * np.maximum(F - P["Fref"], 0.0) ** P["pf"]
    dI = (P["r0"] + P["kappa"] * bracket) * (1 - I)
    dM = P["kM"] * ((1 - P["lam"] * I) * (1 - P["phi"] * _H(uM, P, hform)) - M)
    dF = P["gF"] * (0.4 + 0.4 * P["kappa"] * I - F) - P["cF"] * uF * F
    return np.stack([dI, dM, dF], axis=1)


USE_FAST = True


def simulate(P, u, t_start=AGE_START, t_end=AGE_END, dt=0.125, hform="hill", pform="power", traj=False):
    """Integrate from age 25 to t_end for a batch (n,) of parameter sets and doses u (n,) switched on at t_start.
    Returns the final state (n,3) or, with traj=True, (ages, states[nt, n, 3]) on the dt grid."""
    n = len(P["r0"])
    u = np.asarray(u, float)
    if USE_FAST and not traj:
        import fastmodel
        return fastmodel.final_state(P, u, t_start, t_end, dt, hform, pform, I0, M0, F0, AGE0)
    u = np.broadcast_to(u, (n,)) if u.ndim < 2 else u
    y = np.tile(np.array([I0, M0, F0]), (n, 1))
    nsteps = int(round((t_end - AGE0) / dt))
    ages = AGE0 + dt * np.arange(nsteps + 1)
    out = [y.copy()] if traj else None
    zero = np.zeros_like(u)
    for k in range(nsteps):
        uk = u if (ages[k] >= t_start - 1e-9) else zero      # t_start is a multiple of dt, so the switch is exact
        k1 = rhs(y, uk, P, hform, pform)
        k2 = rhs(y + 0.5 * dt * k1, uk, P, hform, pform)
        k3 = rhs(y + 0.5 * dt * k2, uk, P, hform, pform)
        k4 = rhs(y + dt * k3, uk, P, hform, pform)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if traj:
            out.append(y.copy())
    return (ages, np.array(out)) if traj else y


def sim_scipy(P1, u, t_start=AGE_START, t_end=AGE_END, hform="hill", pform="power"):
    """Reference integration with LSODA for ONE parameter set (arrays of length 1)."""
    from scipy.integrate import solve_ivp
    P = {k: np.array([float(np.ravel(v)[0])]) for k, v in P1.items()}
    pre = solve_ivp(lambda t, y: rhs(y[None, :], 0.0, P, hform, pform)[0], (AGE0, t_start), [I0, M0, F0],
                    method="LSODA", rtol=1e-10, atol=1e-12)
    post = solve_ivp(lambda t, y: rhs(y[None, :], u, P, hform, pform)[0], (t_start, t_end), pre.y[:, -1],
                     method="LSODA", rtol=1e-10, atol=1e-12)
    return post.y[:, -1]


def dose_grid(umax=2.0, n=41):
    return np.linspace(0.0, umax, n)


def dose_response(P, doses=None, **kw):
    """I(T) for every parameter set (rows) and dose (columns), and the full final state."""
    doses = dose_grid() if doses is None else np.asarray(doses, float)
    n = len(P["r0"])
    PP = tile(P, len(doses))
    uu = np.tile(doses, n)
    y = simulate(PP, uu, **kw)
    return y[:, 0].reshape(n, len(doses)), y.reshape(n, len(doses), 3)


def classify(J, doses, thr=1e-3):
    """Interior optimum: grid minimum of I(T) strictly inside the dose range, with benefit (I(0) - min) and
    backfire (I(umax) - min) each above thr. Returns (flag, u*, benefit, backfire)."""
    k = J.argmin(axis=1)
    nd = J.shape[1]
    ar = np.arange(len(J))
    jmin = J[ar, k]
    benefit = J[:, 0] - jmin
    backfire = J[:, -1] - jmin
    flag = (k > 0) & (k < nd - 1) & (benefit > thr) & (backfire > thr)
    return flag, doses[k], benefit, backfire


def refine_opt(J, doses):
    """Parabolic refinement of the grid minimum (for interior minima)."""
    k = J.argmin(axis=1)
    kk = np.clip(k, 1, J.shape[1] - 2)
    ar = np.arange(len(J))
    a, b, c = J[ar, kk - 1], J[ar, kk], J[ar, kk + 1]
    den = a - 2 * b + c
    h = doses[1] - doses[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        off = 0.5 * h * (a - c) / den
    return np.where((k > 0) & (k < J.shape[1] - 1) & (den > 0), doses[kk] + off, doses[k])


_UNTREATED = {}


def metabolic_age(Ival):
    """Age at which the untreated nominal model reaches information loss Ival (inversion of the untreated I(age) curve)."""
    if "c" not in _UNTREATED:
        ages, tr = simulate(make(1), 0.0, traj=True)
        _UNTREATED["c"] = (ages, tr[:, 0, 0])
    ages, I = _UNTREATED["c"]
    return np.interp(Ival, I, ages, left=AGE0, right=AGE_END)


def sim_index(y, u, P):
    G, sig = algebra(y[..., 0], y[..., 1], y[..., 2], u, P)
    return G * sig / y[..., 1]
