# -*- coding: utf-8 -*-
"""All numerical analyses of the manuscript. Every job writes results/<name>.json.

    python analysis.py nominal | decomposition | static | ensemble | oat | phase | sobol | identifiability | all

Seeds: numpy SeedSequence([20261115, job code]); every job can be repeated alone.
"""
import json
import os
import sys
import time

import numpy as np

import model as m

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
os.makedirs(RES, exist_ok=True)
MASTER = 20261115
DOSES = m.dose_grid()
UMAX = float(DOSES[-1])
CHUNK = 2500


def save(name, obj):
    with open(os.path.join(RES, name + ".json"), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    print("   wrote", name + ".json")


def rng_for(code):
    return np.random.default_rng(np.random.SeedSequence([MASTER, code]))


# ------------------------------------------------------------------------------------------------ sampling
def draw(n, prior="logu24", rng=None, u01=None):
    """Return (P, theta) for n parameter sets. theta holds the sampled scale factors / exponents (n, 16)."""
    D = len(m.SCALE) + len(m.EXPON)
    if u01 is None:
        u01 = rng.random((n, D))
    P = m.make(n)
    th = np.zeros((n, D))
    for j, k in enumerate(m.SCALE):
        x = u01[:, j]
        if prior == "logu24":
            f = np.exp(np.log(0.5) + x * (np.log(2.0) - np.log(0.5)))
        elif prior == "logu416":
            f = np.exp(np.log(0.25) + x * (np.log(4.0) - np.log(0.25)))
        elif prior == "lin":
            f = 0.5 + x
        elif prior == "lognorm":
            from scipy.stats import norm
            f = np.exp(0.25 * norm.ppf(np.clip(x, 1e-9, 1 - 1e-9)))
        else:
            raise ValueError(prior)
        P[k] = P[k] * f
        th[:, j] = f
    for j, k in enumerate(m.EXPON):
        x = u01[:, len(m.SCALE) + j]
        lo, hi = (1.0, 4.0) if prior == "logu416" else (1.0, 3.0)
        P[k] = lo + x * (hi - lo)
        th[:, len(m.SCALE) + j] = P[k]
    P["phi"] = np.minimum(P["phi"], 0.98)
    P["lam"] = np.minimum(P["lam"], 0.9)          # keeps the muscle set point (1 - lam I) positive for every I <= 1
    return P, th


def subset(P, idx):
    return {k: v[idx] for k, v in P.items()}


def with_fat(P, psi, pf=1.0):
    """Adiposity penalty r_f Phi(F), scaled so that at the untreated age-75 state it equals psi times the muscle penalty.
    The scale uses each case's own untreated end state."""
    out = {k: v.copy() for k, v in P.items()}
    if psi == 0:
        return out
    y = m.simulate(P, 0.0)
    musc = P["rm"] * np.maximum(1 - y[:, 1], 0) ** P["pm"]
    fat = np.maximum(y[:, 2] - P["Fref"], 1e-9) ** pf
    out["rf"] = psi * musc / fat
    out["pf"] = np.full(len(P["r0"]), pf)
    return out


def run_J(P, doses=DOSES, **kw):
    """Dose-response I(T) for all cases, processed in chunks."""
    n = len(P["r0"])
    J = np.zeros((n, len(doses)))
    for a in range(0, n, CHUNK):
        sl = slice(a, min(n, a + CHUNK))
        J[sl], _ = m.dose_response(subset(P, sl), doses, **kw)
    return J


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [(c - h) / d, (c + h) / d]


def summarize(J, thr=1e-3):
    flag, us, ben, bf = m.classify(J, DOSES, thr)
    uref = m.refine_opt(J, DOSES)
    n = len(J)
    out = dict(n=n, frac=float(flag.mean()), ci95=wilson(int(flag.sum()), n))
    if flag.any():
        out["u_star_median"] = float(np.median(uref[flag]))
        out["u_star_q"] = [float(x) for x in np.quantile(uref[flag], [0.1, 0.9])]
    return out, flag, uref, ben, bf


# ------------------------------------------------------------------------------------------------ nominal
def job_nominal():
    P = m.make(1)
    ages, tr0 = m.simulate(P, 0.0, traj=True)
    J, Y = m.dose_response(P, DOSES)
    flag, us, ben, bf = m.classify(J, DOSES)
    uref = float(m.refine_opt(J, DOSES)[0])
    # dose-response of three outcomes
    S0 = m.sim_index(tr0[:, 0, :], 0.0, P)
    S = m.sim_index(Y[0], DOSES, P)
    age_sim = np.interp(S, S0, ages)
    age_I = m.metabolic_age(J[0])
    G0, sig0 = m.algebra(tr0[:, 0, 0], tr0[:, 0, 1], tr0[:, 0, 2], 0.0, P)
    kopt = int(J[0].argmin())
    ages_t, trt = m.simulate(m.make(1), UMAX * 0 + DOSES[kopt], traj=True)
    Gt, _ = m.algebra(trt[:, 0, 0], trt[:, 0, 1], trt[:, 0, 2], np.where(ages_t >= m.AGE_START, DOSES[kopt], 0.0), P)
    St = m.sim_index(trt[:, 0, :], np.where(ages_t >= m.AGE_START, DOSES[kopt], 0.0), P)
    # SIM along the treatment course, untreated vs optimum vs supratherapeutic
    def course(u):
        a, tr = m.simulate(m.make(1), u, traj=True)
        uu = np.where(a >= m.AGE_START, u, 0.0)
        G, _ = m.algebra(tr[:, 0, 0], tr[:, 0, 1], tr[:, 0, 2], uu, P)
        return dict(I=tr[:, 0, 0], M=tr[:, 0, 1], F=tr[:, 0, 2], G=G, SIM=m.sim_index(tr[:, 0, :], uu, P))
    out = dict(ages=ages, doses=DOSES, I75=J[0], M75=Y[0, :, 1], F75=Y[0, :, 2], SIM75=S, age_from_I=age_I, age_from_SIM=age_sim,
               untreated=dict(I=tr0[:, 0, 0], M=tr0[:, 0, 1], F=tr0[:, 0, 2], G=G0, SIM=S0),
               u_star_grid=float(us[0]), u_star_refined=uref, benefit=float(ben[0]), backfire=float(bf[0]),
               flip=bool(flag[0]), course_opt=course(float(DOSES[kopt])), course_hi=course(2.0), course_u04=course(0.4),
               argmin_I=float(DOSES[J[0].argmin()]), argmin_SIM=float(DOSES[S.argmin()]), argmin_ageSIM=float(DOSES[age_sim.argmin()]),
               untreated_rise_G=float(G0[-1] / G0[0] - 1), untreated_M_loss=float(1 - tr0[-1, 0, 1]), untreated_I75=float(tr0[-1, 0, 0]))
    save("nominal", out)


# ------------------------------------------------------------------------------------------------ pathway decomposition
def dJ_channels(P, u, h=1e-4):
    """Exact decomposition dI(T)/du = d/du_G + d/du_M + d/du_S by central differences on separate channel doses."""
    base = np.array([u, u, u], float)
    res = []
    for c in range(3):
        up, dn = base.copy(), base.copy()
        up[c] += h
        dn[c] -= h
        a = m.simulate(P, np.tile(up, (1, 1)))[0, 0]
        b = m.simulate(P, np.tile(dn, (1, 1)))[0, 0]
        res.append((a - b) / (2 * h))
    return res


def job_decomposition():
    out = {}
    for label, over in (("nominal", {}), ("fat_psi1", "psi1")):
        P = m.make(1)
        if over == "psi1":
            P = with_fat(P, 1.0)
        us = np.linspace(0.025, 1.9, 76)
        ch = np.array([dJ_channels(P, u) for u in us])
        tot = ch.sum(axis=1)
        h = 2e-3
        J = np.array([m.simulate(P, np.array([[u, u, u]]))[0, 0] for u in us])
        Jp = np.array([m.simulate(P, np.array([[u + h] * 3]))[0, 0] for u in us])
        Jm = np.array([m.simulate(P, np.array([[u - h] * 3]))[0, 0] for u in us])
        d1 = (Jp - Jm) / (2 * h)
        d2 = (Jp - 2 * J + Jm) / h ** 2
        # sign changes of the total derivative
        sc = np.where((tot[:-1] < 0) & (tot[1:] >= 0))[0]
        roots = []
        for k in sc:
            u0 = us[k] - tot[k] * (us[k + 1] - us[k]) / (tot[k + 1] - tot[k])
            roots.append(float(u0))
        out[label] = dict(doses=us, glucose=ch[:, 0], muscle=ch[:, 1], storage=ch[:, 2], total=tot, total_check=d1, second=d2,
                          roots=roots, max_abs_check=float(np.abs(tot - d1).max()))
        if roots:
            u0 = roots[0]
            out[label]["second_at_root"] = float(np.interp(u0, us, d2))
    # small-dose behaviour: dJ/du at u -> 0
    P = m.make(1)
    out["small_dose"] = {str(u): dJ_channels(P, u, h=u * 0.1) for u in (0.01, 0.02, 0.05)}
    save("decomposition", out)


# ------------------------------------------------------------------------------------------------ static (quasi-static) criterion
def static_rate(P, I, u):
    """Instantaneous information-loss rate dI/dt / (1 - I) at frozen I when M and F have relaxed to their set points."""
    M = (1 - P["lam"] * I) * (1 - P["phi"] * m._H(u, P, "hill"))
    Fset = 0.4 + 0.4 * P["kappa"] * I
    F = P["gF"] * Fset / (P["gF"] + P["cF"] * u)
    G, sig = m.algebra(I, M, F, u, P)
    return P["r0"] + P["kappa"] * (P["rg"] * G * sig + P["rm"] * m._pen(1 - M, P, "power") + P["rf"] * np.maximum(F - P["Fref"], 0) ** P["pf"]), G, sig, M, F


def job_static():
    P = m.make(1)
    uu = np.linspace(0.0, 2.0, 801)
    out = dict(doses=uu)
    nominal_ratio = float(P["rm"][0] / P["rg"][0])
    out["nominal_ratio"] = nominal_ratio
    # critical ratio rm/rg at which the static rate has a stationary point at dose u
    for I in (0.2, 0.35, 0.5):
        Pn = m.make(len(uu))
        h = 1e-5
        G1, s1, M1, F1 = static_rate(Pn, I, uu + h)[1:]
        G0_, s0, M0_, F0_ = static_rate(Pn, I, np.maximum(uu - h, 0))[1:]
        dd = (uu + h) - np.maximum(uu - h, 0)
        dG = (G1 - G0_) / dd
        dM = (M1 - M0_) / dd
        sig = 1 + P["csig"] * I
        with np.errstate(divide="ignore", invalid="ignore"):
            rc = sig * (-dG) / (P["pm"] * (1 - 0.5 * (M1 + M0_)) ** (P["pm"] - 1) * (-dM))
        R = static_rate(Pn, I, uu)[0]
        out["I=%.2f" % I] = dict(crit_ratio=rc, rate=R, u_min_rate=float(uu[R.argmin()]),
                                 rate_at_0=float(R[0]), rate_min=float(R.min()), rate_at_umax=float(R[-1]))
    save("static", out)


# ------------------------------------------------------------------------------------------------ ensemble family
def job_ensemble():
    out = {}
    t0 = time.time()
    # (1) the original scheme: N = 6000, fixed seed, log-uniform [0.5,2], exponents U[1,3]
    P6, th6 = draw(6000, "logu24", rng_for(1))
    J6 = run_J(P6)
    s6, f6, u6, b6, bf6 = summarize(J6)
    out["original_scheme_N6000"] = s6
    # threshold on 'benefit' and 'backfire' used for the classification
    out["threshold"] = {str(t): summarize(J6, t)[0]["frac"] for t in (1e-4, 5e-4, 1e-3, 5e-3, 1e-2)}
    # (2) the same with a large ensemble, and the ratio dependence
    N = 40000
    P, th = draw(N, "logu24", rng_for(2))
    J = run_J(P)
    s, flag, uref, ben, bf = summarize(J)
    out["base_N%d" % N] = s
    ratio = P["rm"] / P["rg"]
    bins = np.exp(np.linspace(np.log(0.3), np.log(60), 16))
    idx = np.digitize(ratio, bins)
    rb = []
    for b in range(1, len(bins)):
        sel = idx == b
        if sel.sum() > 30:
            rb.append(dict(ratio=float(np.sqrt(bins[b - 1] * bins[b])), frac=float(flag[sel].mean()), n=int(sel.sum())))
    out["ratio_curve"] = rb
    out["u_star_hist"] = dict(edges=np.linspace(0, 2, 41), counts=np.histogram(uref[flag], bins=np.linspace(0, 2, 41))[0])
    out["frac_by_exponent"] = dict(nH=[dict(lo=a, frac=float(flag[(P["nH"] >= a) & (P["nH"] < a + 0.5)].mean())) for a in (1.0, 1.5, 2.0, 2.5)],
                                   pm=[dict(lo=a, frac=float(flag[(P["pm"] >= a) & (P["pm"] < a + 0.5)].mean())) for a in (1.0, 1.5, 2.0, 2.5)])
    # relation between the dynamic optimum and the static criterion is checked in job_static / figures
    out["sample_for_scatter"] = dict(ratio=ratio[:3000], u_star=np.where(flag[:3000], uref[:3000], np.nan), phi=P["phi"][:3000])
    # (3) alternative priors (common design where possible)
    pri = {}
    for j, pr in enumerate(("logu24", "logu416", "lin", "lognorm"), start=10):
        Pp, _ = draw(20000, pr, rng_for(j))
        pri[pr] = summarize(run_J(Pp))[0]
    from scipy.stats import qmc
    u01 = qmc.Sobol(d=len(m.SCALE) + len(m.EXPON), scramble=True, seed=MASTER).random(16384)
    Pq, _ = draw(16384, "logu24", None, u01=u01)
    pri["sobol_sequence_logu24"] = summarize(run_J(Pq))[0]
    out["priors"] = pri
    print("   priors done %.0fs" % (time.time() - t0))
    # (4) time horizon (years of treatment from age 45)
    hor = {}
    Ph, _ = draw(20000, "logu24", rng_for(20))
    for T in (5, 10, 20, 30):
        hor[str(T)] = summarize(run_J(Ph, t_end=m.AGE_START + T))[0]
    out["horizon"] = hor
    print("   horizon done %.0fs" % (time.time() - t0))
    # (5) adiposity penalty: psi = ratio of the fat to the muscle penalty at the untreated age-75 state
    fat = {}
    Pa, _ = draw(20000, "logu24", rng_for(30))
    nomP = m.make(1)
    for psi in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
        Pf = with_fat(Pa, psi)
        sf = summarize(run_J(Pf))[0]
        Pn = with_fat(nomP, psi)
        Jn, _ = m.dose_response(Pn, DOSES)
        fn, un, bn, bfn = m.classify(Jn, DOSES)
        sf["nominal"] = dict(flip=bool(fn[0]), u_star=float(m.refine_opt(Jn, DOSES)[0]), benefit=float(bn[0]), backfire=float(bfn[0]), I_dose=Jn[0])
        fat[str(psi)] = sf
    # convex fat penalty (pf = 2) at psi = 1
    Pf2 = with_fat(Pa, 1.0, pf=2.0)
    fat["psi1_pf2"] = summarize(run_J(Pf2))[0]
    out["adiposity"] = fat
    print("   adiposity done %.0fs" % (time.time() - t0))
    # (6) ablations (common random numbers: the same 20 000 draws)
    abl = {}
    Pb, _ = draw(20000, "logu24", rng_for(40))
    spec = {"complete model": {}, "no muscle-deficit cost (r_m = 0)": {"rm": 0.0}, "no drug-induced muscle loss (phi = 0)": {"phi": 0.0},
            "no direct glucose lowering (a_u = 0)": {"au": 0.0}, "no susceptibility feedback (c_sig = 0)": {"csig": 0.0},
            "no information-to-glucose feedback (a_I = 0)": {"aI": 0.0}, "no storage compartment effect (a_F = c_F = 0)": {"aF": 0.0, "cF": 0.0},
            "no muscle-to-glucose disposal (a_M = 0)": {"aM": 0.0}}
    for lab, ov in spec.items():
        Pk = {k: v.copy() for k, v in Pb.items()}
        for k, v in ov.items():
            Pk[k] = np.full(len(Pk[k]), v)
        sk = summarize(run_J(Pk))[0]
        Pn = m.make(1, **ov)
        Jn, _ = m.dose_response(Pn, DOSES)
        fn, un, bn, bfn = m.classify(Jn, DOSES)
        sk["nominal"] = dict(flip=bool(fn[0]), u_star=float(m.refine_opt(Jn, DOSES)[0]), benefit=float(bn[0]), backfire=float(bfn[0]),
                             I75_untreated=float(Jn[0, 0]), I75_min=float(Jn[0].min()), I75_umax=float(Jn[0, -1]))
        abl[lab] = sk
    out["ablation"] = abl
    print("   ablation done %.0fs" % (time.time() - t0))
    # (7) alternative functional forms (nominal and ensemble)
    forms = {}
    Pc, _ = draw(10000, "logu24", rng_for(50))
    for hf in ("hill", "mm", "exp"):
        for pf_ in ("power", "threshold", "saturating"):
            Jn, _ = m.dose_response(m.make(1), DOSES, hform=hf, pform=pf_)
            fn, un, bn, bfn = m.classify(Jn, DOSES)
            sk = summarize(run_J(Pc, hform=hf, pform=pf_))[0]
            sk["nominal"] = dict(flip=bool(fn[0]), u_star=float(m.refine_opt(Jn, DOSES)[0]), benefit=float(bn[0]), backfire=float(bfn[0]))
            forms["%s|%s" % (hf, pf_)] = sk
    out["forms"] = forms
    print("   forms done %.0fs" % (time.time() - t0))
    save("ensemble", out)


# ------------------------------------------------------------------------------------------------ one-at-a-time
def job_oat():
    out = {}
    scales = np.exp(np.linspace(np.log(0.5), np.log(2.0), 13))
    for k in m.SCALE:
        P = m.make(len(scales))
        P[k] = P[k] * scales
        if k == "phi":
            P[k] = np.minimum(P[k], 0.98)
        J = run_J(P)
        f, us, b, bf = m.classify(J, DOSES)
        out[k] = dict(scale=scales, flip=f, u_star=m.refine_opt(J, DOSES), backfire=bf, benefit=b)
    ex = np.linspace(1, 3, 9)
    for k in m.EXPON:
        P = m.make(len(ex))
        P[k] = ex
        J = run_J(P)
        f, us, b, bf = m.classify(J, DOSES)
        out[k] = dict(scale=ex, flip=f, u_star=m.refine_opt(J, DOSES), backfire=bf, benefit=b)
    save("oat", out)


# ------------------------------------------------------------------------------------------------ phase diagrams
def job_phase():
    out = {}
    n = 41
    a = np.exp(np.linspace(np.log(0.25), np.log(4.0), n))
    ph = np.linspace(0.1, 0.9, n)
    for psi in (0.0, 1.0):
        # (A) muscle-cost scale vs glucose-benefit (drug glucose-lowering potency) scale
        X, Y = np.meshgrid(a, a, indexing="ij")
        P = m.make(n * n)
        P["rm"] = P["rm"] * X.ravel()
        P["au"] = P["au"] * Y.ravel()
        P = with_fat(P, psi) if psi else P
        J = run_J(P)
        f, us, b, bf = m.classify(J, DOSES)
        out["rm_au|psi=%g" % psi] = dict(x=a, y=a, flip=f.reshape(n, n), u_star=m.refine_opt(J, DOSES).reshape(n, n), backfire=bf.reshape(n, n))
        # (B) muscle-cost scale vs maximum fractional muscle loss phi
        X, Y = np.meshgrid(a, ph, indexing="ij")
        P = m.make(n * n)
        P["rm"] = P["rm"] * X.ravel()
        P["phi"] = Y.ravel()
        P = with_fat(P, psi) if psi else P
        J = run_J(P)
        f, us, b, bf = m.classify(J, DOSES)
        out["rm_phi|psi=%g" % psi] = dict(x=a, y=ph, flip=f.reshape(n, n), u_star=m.refine_opt(J, DOSES).reshape(n, n), backfire=bf.reshape(n, n))
    # (C) the original Figure 5: optimal dose against phi, with and without the fat penalty
    P = m.make(len(ph))
    P["phi"] = ph
    J = run_J(P)
    out["phi_curve"] = dict(phi=ph, u_star=m.refine_opt(J, DOSES), flip=m.classify(J, DOSES)[0])
    save("phase", out)


# ------------------------------------------------------------------------------------------------ Sobol
def sobol_indices(YA, YB, YAB, nboot=300, rng=None):
    N, D = YAB.shape
    var = np.var(np.concatenate([YA, YB]), ddof=1)

    def est(ia, ib, iab):
        v = np.var(np.concatenate([YA[ia], YB[ib]]), ddof=1)
        s1 = np.mean(YB[ib][:, None] * (YAB[iab] - YA[ia][:, None]), axis=0) / v
        st = 0.5 * np.mean((YA[ia][:, None] - YAB[iab]) ** 2, axis=0) / v
        return s1, st
    s1, st = est(np.arange(N), np.arange(N), np.arange(N))
    bs1, bst = [], []
    for _ in range(nboot):
        ii = rng.integers(0, N, N)
        a, b = est(ii, ii, ii)
        bs1.append(a)
        bst.append(b)
    bs1, bst = np.array(bs1), np.array(bst)
    return dict(S1=s1, ST=st, S1_ci=np.quantile(bs1, [0.025, 0.975], axis=0), ST_ci=np.quantile(bst, [0.025, 0.975], axis=0), var=float(var))


def job_sobol():
    from scipy.stats import qmc
    D = len(m.SCALE) + len(m.EXPON)
    N = 4096
    sob = qmc.Sobol(d=2 * D, scramble=True, seed=np.random.SeedSequence([MASTER, 60]).generate_state(1)[0] % (2 ** 31)).random(N)
    A, B = sob[:, :D], sob[:, D:]
    names = m.SCALE + m.EXPON

    def outputs(u01):
        P, _ = draw(len(u01), "logu24", None, u01=u01)
        J = run_J(P)
        f, us, ben, bf = m.classify(J, DOSES)
        uref = m.refine_opt(J, DOSES)
        return dict(backfire=bf, u_star=uref, benefit=ben, flip=f.astype(float))
    t0 = time.time()
    oa, ob = outputs(A), outputs(B)
    oab = []
    for i in range(D):
        ABi = A.copy()
        ABi[:, i] = B[:, i]
        oab.append(outputs(ABi))
        print("   sobol column %d/%d %.0fs" % (i + 1, D, time.time() - t0))
    out = dict(names=names, N=N)
    rng = rng_for(61)
    for key in ("backfire", "u_star", "benefit", "flip"):
        YAB = np.stack([o[key] for o in oab], axis=1)
        out[key] = sobol_indices(oa[key], ob[key], YAB, rng=rng)
    save("sobol", out)


# ------------------------------------------------------------------------------------------------ identifiability
def observables(theta_log, names, scenario):
    """Cohort-mean observables (log scale) as functions of the log-parameters. A cohort is (dose from age 45, ages)."""
    P = m.make(1)
    for k, v in zip(names, theta_log):
        P[k] = np.array([np.exp(v)])
    out = []
    for dose, ages_obs in scenario["cohorts"]:
        ages, tr = m.simulate(P, dose, traj=True)
        for a in ages_obs:
            k = int(round((a - m.AGE0) / 0.125))
            uu = dose if a >= m.AGE_START else 0.0
            G, _ = m.algebra(tr[k, 0, 0], tr[k, 0, 1], tr[k, 0, 2], uu, P)
            out += [np.log(float(np.ravel(G)[0])), np.log(tr[k, 0, 1]), np.log(tr[k, 0, 2])]
            if "I" in scenario["obs"]:
                out.append(np.log(tr[k, 0, 0]))
    return np.array(out)


def job_identifiability():
    names = m.SCALE + m.EXPON
    nom = np.array([np.log(m.NOMINAL[k]) for k in names])
    prior_sd = 0.40            # log-uniform over a factor of four has sd ln4/sqrt(12) = 0.40
    U = (0.0, (35, 45, 55, 65, 75))
    T04 = (0.4, (55, 65, 75))
    T15 = (1.5, (55, 65, 75))
    scen = {"G, M, F; untreated cohort only; 5% noise": dict(obs="GMF", cohorts=[U], cv=0.05),
            "G, M, F; untreated + one treated cohort (u = 0.4); 5% noise": dict(obs="GMF", cohorts=[U, T04], cv=0.05),
            "G, M, F; untreated + one treated cohort (u = 0.4); 20% noise": dict(obs="GMF", cohorts=[U, T04], cv=0.20),
            "G, M, F and I; untreated + one treated cohort; 5% noise": dict(obs="GMFI", cohorts=[U, T04], cv=0.05),
            "G, M, F and I; dose-ranging (u = 0.4 and 1.5); 5% noise": dict(obs="GMFI", cohorts=[U, T04, T15], cv=0.05),
            "G, M, F; dose-ranging (u = 0.4 and 1.5); 5% noise": dict(obs="GMF", cohorts=[U, T04, T15], cv=0.05)}
    out = dict(names=names, prior_sd=prior_sd)
    for lab, sc in scen.items():
        h = 1e-3
        cols = []
        for j in range(len(names)):
            tp, tm = nom.copy(), nom.copy()
            tp[j] += h
            tm[j] -= h
            cols.append((observables(tp, names, sc) - observables(tm, names, sc)) / (2 * h))
        Jm = np.array(cols).T
        Fm = Jm.T @ Jm / sc["cv"] ** 2
        w = np.linalg.eigvalsh(Fm)[::-1]
        Fp = Fm + np.eye(len(names)) / prior_sd ** 2
        post_sd = np.sqrt(np.diag(np.linalg.inv(Fp)))
        Sig = np.linalg.inv(Fp)
        c = np.zeros(len(names))
        c[names.index("rm")], c[names.index("rg")] = 1.0, -1.0           # log of the ratio r_m / r_g
        sd_ratio = float(np.sqrt(c @ Sig @ c))
        out[lab] = dict(eigenvalues=w, n_directions_above_prior=int((w > 1 / prior_sd ** 2).sum()),
                        post_sd=post_sd, gain=1 - post_sd / prior_sd, n_obs=int(Jm.shape[0]),
                        ratio_gain=1 - sd_ratio / (prior_sd * np.sqrt(2.0)))
    save("identifiability", out)


JOBS = dict(nominal=job_nominal, decomposition=job_decomposition, static=job_static, ensemble=job_ensemble, oat=job_oat, phase=job_phase,
            sobol=job_sobol, identifiability=job_identifiability)

if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    for w in (list(JOBS) if what == ["all"] else what):
        t = time.time()
        print("== %s" % w)
        JOBS[w]()
        print("   %.1fs" % (time.time() - t))
