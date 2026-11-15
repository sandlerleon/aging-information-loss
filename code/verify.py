# -*- coding: utf-8 -*-
"""Executable checks of the implementation (Appendix B of the manuscript). Writes results/verification.json."""
import json
import os
import sys

import numpy as np

import model as m
import analysis as A

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
OUT = {}


def rec(key, passed, detail):
    OUT[key] = dict(passed=bool(passed), detail=detail)
    print("[%s] %s: %s" % ("PASS" if passed else "FAIL", key, detail))


# V1: reproduction of the quantities reported for the nominal case in the earlier manuscript
P = m.make(1)
ages, tr = m.simulate(P, 0.0, traj=True)
G, _ = m.algebra(tr[:, 0, 0], tr[:, 0, 1], tr[:, 0, 2], 0.0, P)
J, Y = m.dose_response(P)
flag, us, ben, bf = m.classify(J, m.dose_grid())
S0 = m.sim_index(tr[:, 0, :], 0.0, P)
S = m.sim_index(Y[0], m.dose_grid(), P)
age_sim = np.interp(S, S0, ages)
vals = dict(I75=float(tr[-1, 0, 0]), M_loss=float(1 - tr[-1, 0, 1]), G_rise=float(G[-1] / G[0] - 1), u_star=float(us[0]),
            age_SIM_min=float(age_sim.min()), I_min=float(J[0].min()), I_max=float(J[0].max()))
ok = (abs(vals["I75"] - 0.60) < 0.01 and abs(vals["M_loss"] - 0.16) < 0.01 and abs(vals["G_rise"] - 0.50) < 0.02 and abs(vals["u_star"] - 0.40) < 1e-9
      and abs(vals["age_SIM_min"] - 67.0) < 0.5 and abs(vals["I_min"] - 0.575) < 0.005 and abs(vals["I_max"] - 0.62) < 0.005)
rec("V1_reproduces_earlier_nominal_results", ok, "untreated I(75) = %.4f, muscle loss %.3f, glucose rise %.3f, optimum u = %.2f, SIM-based metabolic age at the optimum %.1f, I(75) range %.4f-%.4f (earlier manuscript: 0.60, ~0.16, ~0.50, 0.40, ~67, 0.575-0.62)" % (
    vals["I75"], vals["M_loss"], vals["G_rise"], vals["u_star"], vals["age_SIM_min"], vals["I_min"], vals["I_max"]))

# V2: RK4 against LSODA
d = []
for u in (0.0, 0.2, 0.4, 1.0, 2.0):
    d.append(np.abs(m.sim_scipy(P, u) - m.simulate(P, u)[0]).max())
rec("V2_RK4_against_LSODA", max(d) < 1e-7, "largest difference of the final state over five doses: %.1e (LSODA rtol 1e-10)" % max(d))

# V2b: compiled integrator against the NumPy reference (all functional forms, with and without the adiposity penalty)
Pt, _ = A.draw(1500, "logu416", A.rng_for(7))
ut = A.rng_for(8).uniform(0, 2, 1500)
dmax = 0.0
for hf in ("hill", "mm", "exp"):
    for pf in ("power", "threshold", "saturating"):
        m.USE_FAST = False
        a = m.simulate(Pt, ut, hform=hf, pform=pf)
        m.USE_FAST = True
        b = m.simulate(Pt, ut, hform=hf, pform=pf)
        dmax = max(dmax, float(np.abs(a - b).max()))
Pft = A.with_fat(Pt, 1.0)
m.USE_FAST = False
a = m.simulate(Pft, ut)
m.USE_FAST = True
b = m.simulate(Pft, ut)
dmax = max(dmax, float(np.abs(a - b).max()))
rec("V2b_compiled_integrator_equals_NumPy_reference", dmax < 1e-12, "1500 random parameter sets and doses, 9 functional-form combinations and the adiposity penalty: largest difference of the final state %.1e" % dmax)

# V3: step-size convergence
dts = [0.25, 0.125, 0.0625]
Jd = [m.simulate(P, 0.4, dt=dt)[0, 0] for dt in dts]
rec("V3_step_size_convergence", abs(Jd[1] - Jd[2]) < 1e-7, "I(75) at u = 0.4 for dt = 0.25, 0.125, 0.0625 yr: %s; adopted dt = 0.125, change from halving %.1e" % (
    ", ".join("%.9f" % x for x in Jd), abs(Jd[1] - Jd[2])))

# V4: invariants over random parameter sets
Pr, _ = A.draw(3000, "logu416", A.rng_for(90))
uu = A.rng_for(91).uniform(0, 2, 3000)
ag, trr = m.simulate(Pr, uu, traj=True)
I = trr[:, :, 0]
rec("V4a_information_loss_bounded_monotone", (I.min() >= 0) and (I.max() <= 1) and (np.diff(I, axis=0).min() >= -1e-12),
    "3000 random parameter sets (wide prior): I in [%.3f, %.3f]; smallest step change %.1e (non-negative)" % (I.min(), I.max(), np.diff(I, axis=0).min()))
rec("V4b_muscle_and_storage_positive", (trr[:, :, 1].min() > 0) and (trr[:, :, 1].max() <= 1.0 + 1e-9) and (trr[:, :, 2].min() > 0),
    "M in [%.3f, %.3f]; F minimum %.3f" % (trr[:, :, 1].min(), trr[:, :, 1].max(), trr[:, :, 2].min()))

# V5: exact channel decomposition
dc = json.load(open(os.path.join(RES, "decomposition.json"), encoding="utf-8"))
rec("V5_channel_decomposition_sums_to_total_derivative", dc["nominal"]["max_abs_check"] < 1e-5,
    "largest |sum of channel derivatives - derivative of I(T) along the common dose| over 76 doses: %.1e (nominal), %.1e (with adiposity penalty)" % (
        dc["nominal"]["max_abs_check"], dc["fat_psi1"]["max_abs_check"]))

# V6: the adiposity penalty construction
Pn = A.with_fat(m.make(1), 1.0)
y = m.simulate(m.make(1), 0.0)             # the scaling uses the untreated age-75 state of the corresponding model with r_f = 0
musc = Pn["rm"] * (1 - y[0, 1]) ** Pn["pm"]
fat = Pn["rf"] * (y[0, 2] - Pn["Fref"]) ** Pn["pf"]
P0 = A.with_fat(m.make(1), 0.0)
rec("V6_adiposity_penalty_scaling", abs(float(musc[0]) - float(fat[0])) < 1e-12 and float(P0["rf"][0]) == 0.0,
    "at the untreated age-75 state of the model with r_f = 0 the fat penalty (%.3e) equals psi times the muscle penalty (%.3e) for psi = 1; psi = 0 gives r_f = 0" % (float(fat[0]), float(musc[0])))

# V7: Sobol estimator on the Ishigami function (analytical indices known)
def ishigami(x, a=7.0, b=0.1):
    x = (x - 0.5) * 2 * np.pi
    return np.sin(x[:, 0]) + a * np.sin(x[:, 1]) ** 2 + b * x[:, 2] ** 4 * np.sin(x[:, 0])
from scipy.stats import qmc
N = 8192
sob = qmc.Sobol(d=6, scramble=True, seed=7).random(N)
Aa, Bb = sob[:, :3], sob[:, 3:]
YA, YB = ishigami(Aa), ishigami(Bb)
YAB = np.stack([ishigami(np.where(np.arange(3)[None, :] == i, Bb, Aa)) for i in range(3)], axis=1)
res = A.sobol_indices(YA, YB, YAB, nboot=50, rng=np.random.default_rng(1))
exact_S1 = np.array([0.3139, 0.4424, 0.0])
exact_ST = np.array([0.5576, 0.4424, 0.2437])
e1 = np.abs(np.array(res["S1"]) - exact_S1).max()
e2 = np.abs(np.array(res["ST"]) - exact_ST).max()
rec("V7_sobol_estimator_on_Ishigami", e1 < 0.03 and e2 < 0.03, "N = %d: largest error of first-order indices %.3f, of total indices %.3f (exact values 0.314, 0.442, 0 and 0.558, 0.442, 0.244)" % (N, e1, e2))

# V8: identifiability finite-difference step
names = m.SCALE + m.EXPON
nom = np.array([np.log(m.NOMINAL[k]) for k in names])
sc = dict(obs="GMF", cohorts=[(0.0, (35, 45, 55, 65, 75)), (0.4, (55, 65, 75))], cv=0.05)
def fim(h):
    cols = []
    for j in range(len(names)):
        tp, tm = nom.copy(), nom.copy()
        tp[j] += h
        tm[j] -= h
        cols.append((A.observables(tp, names, sc) - A.observables(tm, names, sc)) / (2 * h))
    Jm = np.array(cols).T
    return Jm.T @ Jm / sc["cv"] ** 2
F1, F2 = fim(1e-3), fim(1e-4)
rel = np.abs(F1 - F2).max() / np.abs(F1).max()
rec("V8_fisher_information_step_independence", rel < 1e-4, "largest relative difference between information matrices from steps 1e-3 and 1e-4: %.1e" % rel)

# V9: grid resolution of the interior-optimum classification
Pc, _ = A.draw(600, "logu24", A.rng_for(92))
J41 = A.run_J(Pc, m.dose_grid(2.0, 41))
J201 = A.run_J(Pc, m.dose_grid(2.0, 201))
f41 = m.classify(J41, m.dose_grid(2.0, 41))[0]
f201 = m.classify(J201, m.dose_grid(2.0, 201))[0]
agree = float((f41 == f201).mean())
rec("V9_dose_grid_resolution", agree > 0.97, "classification with 41 doses agrees with 201 doses in %.1f%% of 600 random parameter sets (fractions %.3f and %.3f)" % (100 * agree, f41.mean(), f201.mean()))

# V10: reproducibility of the stored ensemble statistic
Pq, _ = A.draw(6000, "logu24", A.rng_for(1))
fr = float(m.classify(A.run_J(Pq), m.dose_grid())[0].mean())
if os.path.exists(os.path.join(RES, "ensemble.json")):
    stored = json.load(open(os.path.join(RES, "ensemble.json"), encoding="utf-8"))["original_scheme_N6000"]["frac"]
    rec("V10_ensemble_reproducible", abs(fr - stored) < 1e-12, "re-running the N = 6000 ensemble with its seed gives %.6f, stored %.6f" % (fr, stored))

json.dump(OUT, open(os.path.join(RES, "verification.json"), "w", encoding="utf-8"), indent=1)
print("\n%s" % ("ALL CHECKS PASSED" if all(v["passed"] for v in OUT.values()) else "SOME CHECKS FAILED"))
sys.exit(0 if all(v["passed"] for v in OUT.values()) else 1)
