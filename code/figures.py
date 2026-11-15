# -*- coding: utf-8 -*-
"""Figures of the manuscript, drawn from results/*.json.   python figures.py [name ...]"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
os.makedirs(FIG, exist_ok=True)

BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#2b6a99", "#d9822b", "#3a8c5c", "#b5403a", "#6b4c9a", "#6b7785"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 8.5, "axes.titlepad": 7, "axes.titleweight": "bold",
                     "axes.labelsize": 9, "legend.frameon": False, "legend.fontsize": 8, "savefig.dpi": 300, "figure.dpi": 150,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6})


def L(name):
    return json.load(open(os.path.join(RES, name + ".json"), encoding="utf-8"))


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("   figure", name)


def tag(ax, t):
    ax.text(-0.14, 1.12, t, transform=ax.transAxes, fontsize=11, fontweight="bold")


# ------------------------------------------------------------------------------------------------ Fig 1: structure
def fig_structure():
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")

    def box(x, y, w, h, text, fc, ec, size=8.5, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12", fc=fc, ec=ec, lw=1.2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, fontweight="bold" if bold else "normal")

    def arrow(p, q, text=None, col="#333", rad=0.0, tx=0, ty=0.12):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=11, lw=1.2, color=col, connectionstyle="arc3,rad=%g" % rad))
        if text:
            ax.text((p[0] + q[0]) / 2 + tx, (p[1] + q[1]) / 2 + ty, text, ha="center", fontsize=7.5, color=col)

    box(0.2, 2.4, 2.0, 1.1, "Intervention\ndose  $u$\n(from age 45)", "#f3e3c8", ORANGE, bold=True)
    box(0.2, 4.6, 2.0, 1.0, "Storage  $F$\n(state)", "#f4e0de", RED)
    box(3.6, 4.6, 2.6, 1.0, "Circulating glucose  $G$\n(quasi-steady, derived)", "#dce9f4", BLUE)
    box(3.6, 2.4, 2.6, 1.1, "Information loss  $I$\n(posited latent state)", "#e6e0f0", PURPLE, bold=True)
    box(3.6, 0.3, 2.6, 1.0, "Muscle reserve  $M$\n(state)", "#dcefe3", GREEN)
    box(7.5, 2.4, 2.3, 1.1, "$\\mathrm{SIM} = G\\sigma(I)/M$\n(derived readout)", "#eeeeee", GREY)
    box(7.5, 0.3, 2.3, 1.0, "Metabolic age\n(inverted untreated\nSIM curve)", "#eeeeee", GREY, size=8)
    arrow((1.2, 3.5), (1.2, 4.6), "", RED)
    ax.text(1.35, 4.05, "removes fat", fontsize=7.5, color=RED)
    arrow((2.2, 3.4), (3.6, 4.75), "", BLUE)
    ax.text(2.0, 4.35, "lowers $G$\ndirectly", fontsize=7.5, color=BLUE, ha="right")
    arrow((2.2, 2.7), (3.6, 1.0), "", GREEN)
    ax.text(2.0, 1.5, "removes\nlean mass", fontsize=7.5, color=GREEN, ha="right")
    arrow((2.2, 5.1), (3.6, 5.1), "", RED)
    ax.text(2.9, 5.2, "enters $G$", fontsize=7.5, color=RED, ha="center")
    arrow((4.3, 4.6), (4.3, 3.5), "", BLUE)
    ax.text(4.2, 4.05, "rate $\\propto G\\,\\sigma$", fontsize=7.5, color=BLUE, ha="right")
    arrow((5.4, 3.5), (5.4, 4.6), "", PURPLE)
    ax.text(5.5, 4.05, "$I$ raises $G$ and $\\sigma$", fontsize=7.5, color=PURPLE, ha="left")
    arrow((4.9, 1.3), (4.9, 2.4), "", GREEN)
    ax.text(5.0, 1.85, "rate $\\propto (1-M)^{p_m}$", fontsize=7.5, color=GREEN, ha="left")
    arrow((6.2, 2.95), (7.5, 2.95), "", GREY)
    arrow((8.65, 2.4), (8.65, 1.3), "", GREY)
    ax.text(5.0, -0.15, "Arrows are the causal terms of Eqs. (1)–(5); the pair of arrows between $G$ and $I$ is the feedback loop.", ha="center", fontsize=7.5, color=GREY)
    save(fig, "Fig1_structure")


# ------------------------------------------------------------------------------------------------ Fig 2: untreated life course
def fig_lifecourse():
    N = L("nominal")
    a = np.array(N["ages"])
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for k, c, lab in (("I", PURPLE, "information loss $I$"), ("M", GREEN, "muscle reserve $M$"), ("F", RED, "storage $F$")):
        ax[0].plot(a, N["untreated"][k], color=c, lw=1.8, label=lab)
    ax[0].set_xlabel("age (years)"); ax[0].set_ylabel("state (normalised)"); ax[0].legend(loc="center right"); tag(ax[0], "A")
    ax[0].set_title("Untreated life course")
    G = np.array(N["untreated"]["G"])
    ax[1].plot(a, G, color=BLUE, lw=1.8, label="glucose $G$")
    ax[1].plot(a, np.array(N["untreated"]["SIM"]) / np.array(N["untreated"]["SIM"])[0] * G[0], color=GREY, lw=1.4, ls="--", label="SIM (scaled to $G$ at 25)")
    ax[1].set_xlabel("age (years)"); ax[1].set_ylabel("circulating glucose (units of $G_0$)"); ax[1].legend(loc="upper left"); tag(ax[1], "B")
    ax[1].set_title("Glucose rises as an output of the loop")
    fig.tight_layout()
    save(fig, "Fig2_lifecourse")


# ------------------------------------------------------------------------------------------------ Fig 3: interior optimum + channels
def fig_optimum():
    N, D = L("nominal"), L("decomposition")
    u = np.array(N["doses"])
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax[0].plot(u, N["I75"], color=PURPLE, lw=2)
    us = N["u_star_refined"]
    ax[0].axvline(us, color=ORANGE, ls="--", lw=1.2)
    ax[0].plot([us], [min(N["I75"])], "o", color=ORANGE)
    ax[0].text(us + 0.05, N["I75"][0] - 0.001, "optimum  $u^*\\approx%.2f$" % us, color=ORANGE, fontsize=8, va="top")
    ax[0].axvspan(0, us, color=GREEN, alpha=0.08); ax[0].axvspan(us, 2.0, color=RED, alpha=0.06)
    ax[0].text(0.05, max(N["I75"]) - 0.002, "benefit\nwins", color=GREEN, fontsize=8, va="top")
    ax[0].text(1.45, min(N["I75"]) + 0.004, "muscle cost\nwins", color=RED, fontsize=8)
    ax[0].set_xlabel("dose $u$ (from age 45)"); ax[0].set_ylabel("information loss at age 75, $I(75)$"); tag(ax[0], "A")
    ax[0].set_title("Information loss is U-shaped in dose")
    x = np.array(D["nominal"]["doses"])
    ax[1].plot(x, D["nominal"]["glucose"], color=BLUE, lw=1.6, label="glucose channel")
    ax[1].plot(x, D["nominal"]["storage"], color=RED, lw=1.6, label="storage channel")
    ax[1].plot(x, D["nominal"]["muscle"], color=GREEN, lw=1.6, label="muscle channel")
    ax[1].plot(x, D["nominal"]["total"], color="k", lw=2.2, label="total")
    ax[1].axhline(0, color="k", lw=0.8)
    if D["nominal"]["roots"]:
        ax[1].axvline(D["nominal"]["roots"][0], color=ORANGE, ls="--", lw=1.2)
    ax[1].set_xlabel("dose $u$"); ax[1].set_ylabel("$\\partial I(75)/\\partial u$"); ax[1].legend(loc="lower right", fontsize=7.5)
    ax[1].set_title("Pathway decomposition of the slope"); tag(ax[1], "B")
    fig.tight_layout()
    save(fig, "Fig3_optimum")


# ------------------------------------------------------------------------------------------------ Fig 4: three outcomes
def fig_outcomes():
    N = L("nominal")
    u = np.array(N["doses"])
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))
    ax[0].plot(u, N["age_from_I"], color=PURPLE, lw=1.8, label="from $I$ (direct)")
    ax[0].plot(u, N["age_from_SIM"], color=ORANGE, lw=1.8, label="from SIM (proxy)")
    ax[0].set_xlabel("dose $u$"); ax[0].set_ylabel("metabolic age at 75 (years)"); ax[0].legend(); ax[0].set_title("Two readings of age"); tag(ax[0], "A")
    ax[0].set_ylim(64, 76)
    a = np.array(N["ages"])
    for key, c, lab in (("untreated", GREY, "untreated"), ("course_u04", ORANGE, "$u = 0.4$"), ("course_hi", RED, "$u = 2.0$")):
        d = N[key]
        ax[1].plot(a, np.array(d["SIM"]) if key != "untreated" else N["untreated"]["SIM"], color=c, lw=1.6, label=lab)
    ax[1].axvline(45, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel("age (years)"); ax[1].set_ylabel("SIM"); ax[1].legend(); ax[1].set_title("SIM along the course"); tag(ax[1], "B")
    ax[2].plot(a, np.array(N["untreated"]["M"]), color=GREY, lw=1.6, label="untreated")
    ax[2].plot(a, N["course_u04"]["M"], color=ORANGE, lw=1.6, label="$u = 0.4$")
    ax[2].plot(a, N["course_hi"]["M"], color=RED, lw=1.6, label="$u = 2.0$")
    ax[2].axvline(45, color="k", ls=":", lw=0.8)
    ax[2].set_xlabel("age (years)"); ax[2].set_ylabel("muscle reserve $M$"); ax[2].legend(); ax[2].set_title("Muscle along the course"); tag(ax[2], "C")
    fig.tight_layout()
    save(fig, "Fig4_outcomes")


# ------------------------------------------------------------------------------------------------ Fig 5: static criterion, ratio, phi
def fig_criterion():
    S, E, P = L("static"), L("ensemble"), L("phase")
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.8))
    u = np.array(S["doses"])
    for key, c in (("I=0.20", BLUE), ("I=0.35", PURPLE), ("I=0.50", RED)):
        rc = np.array(S[key]["crit_ratio"], float)
        ax[0].plot(u, rc, color=c, lw=1.6, label="$I = %s$" % key[2:])
    ax[0].axhline(S["nominal_ratio"], color="k", ls="--", lw=1.0)
    ax[0].text(1.2, S["nominal_ratio"] * 1.08, "nominal $r_m/r_g$", fontsize=7.5)
    ax[0].set_yscale("log"); ax[0].set_ylim(0.3, 300); ax[0].set_xlim(0, 2)
    ax[0].set_xlabel("dose $u$"); ax[0].set_ylabel("critical ratio $r_m/r_g$"); ax[0].legend(); tag(ax[0], "A")
    ax[0].set_title("Static condition, Eq. (11)")
    rc = E["ratio_curve"]
    ax[1].plot([r["ratio"] for r in rc], [100 * r["frac"] for r in rc], "o-", color=PURPLE, ms=3.5)
    ax[1].axvline(S["nominal_ratio"], color="k", ls="--", lw=1.0)
    ax[1].set_xscale("log"); ax[1].set_xlabel("sampled $r_m/r_g$"); ax[1].set_ylabel("interior optimum present (%)"); tag(ax[1], "B")
    ax[1].set_title("Existence vs. coupling ratio")
    pc = P["phi_curve"]
    ph, us, fl = np.array(pc["phi"]), np.array(pc["u_star"]), np.array(pc["flip"])
    ax[2].plot(ph[fl], us[fl], "o-", color=GREEN, ms=3.5)
    ax[2].plot(ph[~fl], us[~fl], "s", color=GREY, ms=3.5, label="no interior optimum\n(minimum at $u_{max}$)")
    ax[2].set_ylim(0, 2.5)
    ax[2].set_xlabel("maximum fractional muscle loss $\\varphi$"); ax[2].set_ylabel("optimal dose $u^*$"); ax[2].legend(loc="center right", fontsize=7); tag(ax[2], "C")
    ax[2].set_title("Effect of muscle sparing")
    fig.tight_layout()
    save(fig, "Fig5_criterion")


# ------------------------------------------------------------------------------------------------ Fig 6: ensemble robustness
def fig_ensemble():
    E = L("ensemble")
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.8))
    h = E["u_star_hist"]
    ed = np.array(h["edges"])
    ax[0].bar(0.5 * (ed[:-1] + ed[1:]), h["counts"], width=ed[1] - ed[0], color=PURPLE, alpha=0.85)
    ax[0].set_xlabel("optimal dose $u^*$"); ax[0].set_ylabel("parameter sets"); tag(ax[0], "A")
    ax[0].set_title("Optimum present in %.0f%%" % (100 * E["base_N40000"]["frac"]))
    lab = {"logu24": "log-u.\n0.5–2", "logu416": "log-u.\n0.25–4", "lin": "unif.\n0.5–1.5", "lognorm": "log-n.\n0.25", "sobol_sequence_logu24": "Sobol\nseq."}
    keys = list(lab)
    fr = [100 * E["priors"][k]["frac"] for k in keys]
    ci = np.array([[100 * E["priors"][k]["ci95"][0], 100 * E["priors"][k]["ci95"][1]] for k in keys])
    ax[1].bar(range(len(keys)), fr, yerr=[np.array(fr) - ci[:, 0], ci[:, 1] - np.array(fr)], color=BLUE, capsize=2)
    ax[1].set_xticks(range(len(keys))); ax[1].set_xticklabels([lab[k] for k in keys], fontsize=6.5)
    ax[1].set_ylabel("interior optimum present (%)"); tag(ax[1], "B"); ax[1].set_title("Sampling prior")
    Ts = ["5", "10", "20", "30"]
    ax[2].errorbar([int(t) for t in Ts], [100 * E["horizon"][t]["frac"] for t in Ts],
                   yerr=[[100 * (E["horizon"][t]["frac"] - E["horizon"][t]["ci95"][0]) for t in Ts], [100 * (E["horizon"][t]["ci95"][1] - E["horizon"][t]["frac"]) for t in Ts]],
                   marker="o", color=ORANGE, capsize=2)
    ax[2].set_xlabel("treatment duration (years)"); ax[2].set_ylabel("interior optimum present (%)"); tag(ax[2], "C"); ax[2].set_title("Treatment horizon")
    fig.tight_layout()
    save(fig, "Fig6_ensemble")


# ------------------------------------------------------------------------------------------------ Fig 7: ablation + adiposity
def fig_ablation():
    E = L("ensemble")
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 3.0), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))
    ab = E["ablation"]
    labs = list(ab)
    short = {"complete model": "complete model", "no muscle-deficit cost (r_m = 0)": "no muscle cost $r_m$", "no drug-induced muscle loss (phi = 0)": "no muscle loss $\\varphi$",
             "no direct glucose lowering (a_u = 0)": "no direct drug effect $a_u$", "no susceptibility feedback (c_sig = 0)": "no susceptibility $c_\\sigma$",
             "no information-to-glucose feedback (a_I = 0)": "no $I\\to G$ feedback $a_I$", "no storage compartment effect (a_F = c_F = 0)": "no storage channel",
             "no muscle-to-glucose disposal (a_M = 0)": "no muscle disposal $a_M$"}
    y = np.arange(len(labs))[::-1]
    ax[0].barh(y, [100 * ab[k]["frac"] for k in labs], color=[PURPLE if k == "complete model" else BLUE for k in labs])
    ax[0].set_yticks(y); ax[0].set_yticklabels([short[k] for k in labs], fontsize=7)
    ax[0].set_xlabel("interior optimum present (%)"); tag(ax[0], "A"); ax[0].set_title("Ablation of model elements")
    ad = E["adiposity"]
    psis = ["0.0", "0.25", "0.5", "1.0", "2.0", "4.0"]
    xs = [float(p) for p in psis]
    ax[1].errorbar(xs, [100 * ad[p]["frac"] for p in psis], yerr=[[100 * (ad[p]["frac"] - ad[p]["ci95"][0]) for p in psis], [100 * (ad[p]["ci95"][1] - ad[p]["frac"]) for p in psis]],
                   marker="o", color=RED, capsize=2)
    ax[1].set_xlabel("fat-to-muscle penalty ratio $\\psi$"); ax[1].set_ylabel("interior optimum present (%)"); tag(ax[1], "B"); ax[1].set_title("Adiposity penalty")
    u = np.linspace(0, 2, 41)
    cm = plt.cm.viridis(np.linspace(0, 0.9, len(psis)))
    for p, c in zip(psis, cm):
        ax[2].plot(u, ad[p]["nominal"]["I_dose"], color=c, lw=1.4, label="$\\psi$ = %s" % p)
    ax[2].set_xlabel("dose $u$"); ax[2].set_ylabel("$I(75)$, nominal case"); ax[2].legend(fontsize=6.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.28), columnspacing=0.8, handlelength=1.2); tag(ax[2], "C"); ax[2].set_title("Nominal dose–response")
    fig.tight_layout()
    save(fig, "Fig7_ablation")


# ------------------------------------------------------------------------------------------------ Fig 8: Sobol
def fig_sobol():
    S = L("sobol")
    names = S["names"]
    pretty = {"r0": "$r_0$", "rg": "$r_g$", "rm": "$r_m$", "csig": "$c_\\sigma$", "aI": "$a_I$", "aF": "$a_F$", "aM": "$a_M$", "au": "$a_u$", "lam": "$\\lambda$", "kM": "$k_M$",
              "phi": "$\\varphi$", "ku": "$k_u$", "gF": "$g_F$", "cF": "$c_F$", "nH": "$n_H$", "pm": "$p_m$"}
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.0), sharey=True)
    for a_, key, title, t in ((ax[0], "backfire", "Backfire depth $I(u_{max}) - \\min I$", "A"), (ax[1], "u_star", "Optimal dose $u^*$", "B")):
        S1, ST = np.array(S[key]["S1"]), np.array(S[key]["ST"])
        order = np.argsort(-ST)[:10][::-1]
        y = np.arange(len(order))
        a_.barh(y + 0.2, ST[order], height=0.38, color=BLUE, label="total $S_T$")
        a_.barh(y - 0.2, S1[order], height=0.38, color=ORANGE, label="first-order $S_1$")
        a_.set_yticks(y); a_.set_yticklabels([pretty[names[i]] for i in order])
        a_.set_xlabel("Sobol index"); a_.set_title(title); tag(a_, t)
    ax[0].legend(loc="lower right")
    fig.tight_layout()
    save(fig, "Fig8_sobol")


# ------------------------------------------------------------------------------------------------ Fig 9: phase diagrams
def fig_phase():
    P = L("phase")
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.2), gridspec_kw=dict(hspace=0.5, wspace=0.32))
    for r, psi in enumerate(("0", "1")):
        for c, key in enumerate(("rm_au|psi=%s.0" % psi if False else "rm_au|psi=%g" % float(psi), "rm_phi|psi=%g" % float(psi))):
            d = P[key]
            x, y = np.array(d["x"]), np.array(d["y"])
            fl = np.array(d["flip"], float)
            us = np.array(d["u_star"], float)
            Z = np.where(fl > 0, us, np.nan)
            a_ = ax[r, c]
            a_.set_facecolor("#e9ecef")
            im = a_.pcolormesh(x, y, Z.T, cmap="viridis", vmin=0, vmax=2.0, shading="auto")
            a_.set_xscale("log")
            if c == 0:
                a_.set_yscale("log")
            a_.contour(x, y, fl.T, levels=[0.5], colors="k", linewidths=1.2)
            a_.plot([1], [1] if c == 0 else [0.5], "*", color="white", mec="k", ms=11)
            a_.set_xlabel("muscle-cost scale  $r_m / r_m^{nom}$")
            a_.set_ylabel("direct glucose-lowering scale  $a_u / a_u^{nom}$" if c == 0 else "maximum muscle loss  $\\varphi$")
            a_.set_title("%s adiposity penalty  ($\\psi$ = %s)" % ("No" if psi == "0" else "With", psi), fontsize=8.5)
            tag(a_, "ABCD"[r * 2 + c])
            a_.grid(False)
            from matplotlib.ticker import NullFormatter
            a_.xaxis.set_minor_formatter(NullFormatter())
            a_.set_xticks([0.25, 0.5, 1, 2, 4]); a_.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
            if c == 0:
                a_.yaxis.set_minor_formatter(NullFormatter())
                a_.set_yticks([0.25, 0.5, 1, 2, 4]); a_.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
    cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label("optimal dose $u^*$ (grey: no interior optimum)")
    save(fig, "Fig9_phase")


# ------------------------------------------------------------------------------------------------ Fig 10: identifiability
def fig_ident():
    I = L("identifiability")
    names = I["names"]
    scen = [k for k in I if isinstance(I[k], dict)]
    pretty = {"r0": "$r_0$", "rg": "$r_g$", "rm": "$r_m$", "csig": "$c_\\sigma$", "aI": "$a_I$", "aF": "$a_F$", "aM": "$a_M$", "au": "$a_u$", "lam": "$\\lambda$", "kM": "$k_M$",
              "phi": "$\\varphi$", "ku": "$k_u$", "gF": "$g_F$", "cF": "$c_F$", "nH": "$n_H$", "pm": "$p_m$"}
    M = np.array([I[k]["gain"] for k in scen])
    fig, ax = plt.subplots(figsize=(7.2, 3.1))
    im = ax.imshow(M, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([pretty[n] for n in names])
    short = ["G, M, F\nuntreated only", "G, M, F\n+ one dose, 5%", "G, M, F\n+ one dose, 20%", "G, M, F, I\n+ one dose", "G, M, F, I\ndose-ranging", "G, M, F\ndose-ranging"]
    ax.set_yticks(range(len(scen))); ax.set_yticklabels(short, fontsize=7)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=6.3, color="white" if M[i, j] > 0.5 else "k")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("information gain  $1 - \\mathrm{SD}_{post}/\\mathrm{SD}_{prior}$")
    save(fig, "Fig10_identifiability")


JOBS = dict(structure=fig_structure, lifecourse=fig_lifecourse, optimum=fig_optimum, outcomes=fig_outcomes, criterion=fig_criterion, ensemble=fig_ensemble,
            ablation=fig_ablation, sobol=fig_sobol, phase=fig_phase, ident=fig_ident)

if __name__ == "__main__":
    what = sys.argv[1:] or list(JOBS)
    for w in what:
        JOBS[w]()
