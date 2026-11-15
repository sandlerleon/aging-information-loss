# -*- coding: utf-8 -*-
# Part 0: derived numbers used in the text, front matter and abstract.
import model as mm

_u = np.array(NOM["doses"])
USTAR = NOM["u_star_refined"]
I75_0 = NOM["untreated_I75"]
I75_MIN = min(NOM["I75"])
BENEFIT = NOM["benefit"]
BACKFIRE = NOM["backfire"]
I75_UMAX = NOM["I75"][-1]
AGE_I_OPT = float(np.min(NOM["age_from_I"]))
AGE_SIM_OPT = float(np.min(NOM["age_from_SIM"]))
ARGMIN_SIM = NOM["argmin_SIM"]
ARGMIN_AGESIM = NOM["argmin_ageSIM"]
ARGMIN_SIM_REF = float(mm.refine_opt(np.array([NOM["SIM75"]]), _u)[0])
G0_, G75_ = NOM["untreated"]["G"][0], NOM["untreated"]["G"][-1]

dn = DEC["nominal"]
_x = np.array(dn["doses"])
ROOT = dn["roots"][0]
SLOPE_G = float(np.interp(ROOT, _x, dn["glucose"]))
SLOPE_M = float(np.interp(ROOT, _x, dn["muscle"]))
SLOPE_S = float(np.interp(ROOT, _x, dn["storage"]))
SD = DEC["small_dose"]["0.01"]
SLOPE0 = sum(SD)
SECOND = dn["second_at_root"]

# static criterion
_uu = np.array(STA["doses"])
RHO = {}
for _k in ("I=0.20", "I=0.35", "I=0.50"):
    _rc = np.array(STA[_k]["crit_ratio"], float)
    _fin = np.isfinite(_rc) & (_uu > 0.05)
    _i = int(np.argmin(np.where(_fin, _rc, np.inf)))
    RHO[_k] = dict(min=float(_rc[_i]), at=float(_uu[_i]), umin=STA[_k]["u_min_rate"], at_umax=float(_rc[-1]))
RATIO_NOM = STA["nominal_ratio"]

EB = ENS["base_N40000"]
EO = ENS["original_scheme_N6000"]
PRI = ENS["priors"]
HOR = ENS["horizon"]
ADI = ENS["adiposity"]
ABL = ENS["ablation"]
FRM = ENS["forms"]
_pri_fr = [PRI[k]["frac"] for k in PRI]
_hor_fr = [HOR[k]["frac"] for k in HOR]

ANSWERS = (
    "First, at the nominal parameters the information loss at age 75 is U-shaped in dose with a minimum at $u^* = %.2f$, the slope of the dose–response splits exactly into a glucose, a storage and a muscle pathway, "
    "and a closed-form stationarity condition relates the optimum to the ratio of the muscle and glucose aging couplings (Sections 4.2 and 4.4). "
    "Second, an interior optimum is present in %s of 40,000 sampled parameterizations, with a range of %s–%s across five sampling priors and %s–%s across treatment horizons of 5 to 30 years; "
    "it requires a muscle-deficit cost and a threshold-like muscle-loss response, and it is rare when lean-mass loss is assumed to start proportionally to dose (Sections 4.5–4.9). "
    "Third, an adiposity penalty deepens the benefit but moves the nominal optimum by less than 0.02 (Section 4.7). "
    "Fourth, the index SIM overstates the benefit relative to the information-loss variable and has a different optimum, and the ratio of couplings that controls the optimum is only weakly constrained by glucose, muscle and storage measurements, even in dose-ranging cohorts (Sections 4.3 and 4.11)."
) % (USTAR, pct(EB["frac"], 0), pct(min(_pri_fr), 0), pct(max(_pri_fr), 0), pct(min(_hor_fr), 0), pct(max(_hor_fr), 0))

FATFACTOR = "%.1f" % (ADI["4.0"]["nominal"]["benefit"] / ADI["0.0"]["nominal"]["benefit"])
RHOMIN = ("with $\\rho_{\\min}$ = %.2f, %.2f and %.2f at $I$ = 0.2, 0.35 and 0.5 against a nominal $r_m/r_g$ of %.2f, and the stationary dose of the frozen-$I$ rate at the nominal ratio is $u$ = %.2f–%.2f, against %.2f for the full model"
          % (RHO["I=0.20"]["min"], RHO["I=0.35"]["min"], RHO["I=0.50"]["min"], RATIO_NOM, RHO["I=0.20"]["umin"], RHO["I=0.50"]["umin"], USTAR))

# ---------------------------------------------------------------------------------------------- front matter
add(("title", TITLE))
add(("author", "Leon Sandler"))
add(("small", "Independent Researcher, Northbrook, Illinois, USA"))
add(("small", "ORCID: https://orcid.org/0009-0007-4584-808X"))
add(("small", "Corresponding author: sandler.leon@gmail.com"))
add(("blank",))
add(("h1", "Abstract"))
ABSTRACT = (
    "Glucose-lowering agents that also remove lean mass pose a quantitative question that aging biomarkers cannot answer by themselves: at what dose does the cost to muscle reserve outweigh the benefit of lower glucose? "
    "This study examines the question in a minimal dynamical model in which a latent information-loss variable, a muscle reserve and an energy-storage compartment are coupled through circulating glucose, "
    "and a single dose lowers glucose, removes lean mass and reduces storage. The study is computational: no data were used, no parameter was fitted, and every result is conditional on tabulated assumptions. "
    "At the nominal parameters, information loss at age 75 is U-shaped in dose with a minimum at $u^* = %.2f$. The slope of the dose–response splits exactly into glucose, storage and muscle pathways, and a closed-form stationarity condition relates the optimum to the ratio of the muscle and glucose aging couplings. "
    "An interior optimum is present in %s of 40,000 sampled parameterizations (%s–%s across five sampling priors; %s–%s across treatment horizons of 5 to 30 years). It requires a muscle-deficit cost and a threshold-like muscle-loss response, and it survives an adiposity penalty, which deepens the benefit but moves the nominal optimum by less than 0.02. "
    "The index SIM, glucose times susceptibility per unit muscle, overstates the benefit relative to information loss (metabolic age %.1f against %.1f years) and has a different optimum. "
    "The coupling ratio that controls the optimum is only weakly constrained by glucose, muscle and storage measurements, even in dose-ranging cohorts. A dose-ranging experiment that would refute the model is specified."
) % (USTAR, pct(EB["frac"], 0), pct(min(_pri_fr), 0), pct(max(_pri_fr), 0), pct(min(_hor_fr), 0), pct(max(_hor_fr), 0), AGE_SIM_OPT, AGE_I_OPT)
add(("p0", ABSTRACT))
add(("h1", "Keywords"))
add(("p0", "Theoretical biology; Dynamical systems; Aging; Sarcopenia; Glucose-lowering therapy; Global sensitivity analysis"))
