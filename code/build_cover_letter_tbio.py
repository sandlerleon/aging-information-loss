# -*- coding: utf-8 -*-
"""Cover letter for the submission to Theory in Biosciences. Numbers are read from results/*.json.

    python build_cover_letter_tbio.py   ->  ../manuscript/Cover_Letter_TBIO.docx
"""
import json
import os
import re

import numpy as np
from docx import Document
from docx.shared import Pt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
OUT = os.path.abspath(os.path.join(HERE, "..", "manuscript"))


def L(n):
    return json.load(open(os.path.join(RES, n + ".json"), encoding="utf-8"))


NOM, ENS, DEC, IDE = L("nominal"), L("ensemble"), L("decomposition"), L("identifiability")
pri = [v["frac"] for v in ENS["priors"].values()]
hor = [v["frac"] for v in ENS["horizon"].values()]
sc = [k for k in IDE if isinstance(IDE[k], dict)]
TITLE = "Muscle–glucose trade-offs in metabolic aging: a dynamical-systems model of an interior optimal dose for glucose-lowering interventions"

doc = Document()
st = doc.styles["Normal"]
st.font.name = "Calibri"
st.font.size = Pt(11)
for s in doc.sections:
    s.left_margin = s.right_margin = Pt(72)
    s.top_margin = s.bottom_margin = Pt(64)


def rich(p, text):
    for k, part in enumerate(re.split(r"_\{([^}]*)\}", text)):
        r = p.add_run(part)
        r.font.subscript = k % 2 == 1


def para(text, bold=False, after=6):
    p = doc.add_paragraph()
    if bold:
        p.add_run(text).bold = True
    else:
        rich(p, text)
    p.paragraph_format.space_after = Pt(after)
    return p


def bullet(text):
    p = doc.add_paragraph(style="List Bullet")
    rich(p, text)
    p.paragraph_format.space_after = Pt(3)


for line in ("Leon Sandler", "Independent Researcher, Northbrook, Illinois, USA", "sandler.leon@gmail.com", "ORCID: https://orcid.org/0009-0007-4584-808X"):
    para(line, after=0)
para("")
para("15 November 2026")
para("The Editors\nTheory in Biosciences", after=10)
para("Submission of an Original Article: “%s”" % TITLE, bold=True, after=10)
para("Dear Editors,")
para("I submit the enclosed manuscript, “%s”, for consideration as an Original Article in Theory in Biosciences." % TITLE)
para("Glucose-lowering agents that also suppress appetite lower lean mass as well as fat, and muscle is the main site of glucose disposal. Aging biomarkers describe the state of the "
     "organism but do not say at what dose the cost to muscle reserve outweighs the benefit of lower glucose. The manuscript asks that question in a minimal dynamical model: a latent "
     "information-loss variable, a muscle reserve and a storage compartment are coupled through circulating glucose, and one dose acts on glucose, muscle and storage. "
     "The study is computational. No data were used and no parameter was fitted, and the paper states throughout which quantities are assumed, which are derived and which are only proposed. The dose is a dimensionless model quantity: the paper does not establish a therapeutic dosing threshold for any specific drug.")
para("Main findings", bold=True, after=3)
dn = DEC["nominal"]
bullet("At the nominal parameters, information loss at age 75 is U-shaped in dose with a minimum at u* = %.2f. The slope of the dose–response splits exactly into a glucose, a storage and a muscle "
       "pathway, and a closed-form stationarity condition relates the optimum to the ratio of the muscle and glucose aging couplings." % NOM["u_star_refined"])
bullet("An interior optimum is present in %.0f%% of 40,000 sampled parameterizations, %.0f–%.0f%% across five sampling priors and %.0f–%.0f%% across treatment horizons of 5 to 30 years. "
       "It requires a muscle-deficit cost and a threshold-like muscle-loss response and is rare when lean-mass loss starts in proportion to dose; the paper reports these dependencies as results, not as caveats."
       % (100 * ENS["base_N40000"]["frac"], 100 * min(pri), 100 * max(pri), 100 * min(hor), 100 * max(hor)))
bullet("An adiposity-associated aging cost, which the earlier formulation lacked, deepens the benefit by a factor of %.1f at the largest penalty examined but moves the nominal optimum by less than 0.02."
       % (ENS["adiposity"]["4.0"]["nominal"]["benefit"] / ENS["adiposity"]["0.0"]["nominal"]["benefit"]))
bullet("The index SIM (glucose times susceptibility per unit muscle) overstates the benefit relative to the information-loss variable: the metabolic age at the optimum falls to %.1f years when computed from SIM and to %.1f years when computed "
       "from information loss. SIM is therefore presented as an unvalidated candidate surrogate with a three-level definition, not as an aging clock."
       % (float(np.min(NOM["age_from_SIM"])), float(np.min(NOM["age_from_I"]))))
bullet("Variance-based sensitivity indices and a local identifiability analysis show which parameters control the result and that the ratio of couplings is only weakly constrained by glucose, muscle and storage "
       "measurements even in dose-ranging cohorts (information gain %.2f). A dose-ranging experiment and the outcomes that would refute the model are specified." % IDE[sc[4]]["ratio_gain"])
para("I believe the work suits the journal because it is a mathematical and computational analysis of a biological trade-off, with explicit derivations, tests of robustness and a stated route to falsification.")
para("Previous review and disclosures", bold=True, after=3)
para("An earlier version of this work was submitted to Biological Theory and was not accepted after review by two referees and an associate editor. Their comments concerned the presentation, the "
     "reliability of the references, the introduction of too many ideas at once, the late statement of the central question, the unclear status of SIM and of the information-loss variable, "
     "and the absence of an adiposity cost. The present manuscript is a rewrite and a re-analysis in response. It is organized around one mechanism, the muscle–glucose trade-off; the comparison "
     "between organisms and the synthetic-cell discussion were removed; every biological statement is tied to a reference, and every reference was checked against its Crossref record; the "
     "assumed, derived and proposed quantities are tabulated; the interior optimum is derived and decomposed analytically; and sensitivity to the sampling prior, the treatment horizon, "
     "the functional forms, the removal of each model element and the adiposity penalty is reported, including the cases in which the optimum disappears. The manuscript has not been published "
     "and is not under consideration by any other journal. The use of AI-assisted tools is described in the manuscript. The code, results and scripts that regenerate every number are public "
     "(https://github.com/sandlerleon/aging-information-loss, archived on Zenodo, version 1.0.0, DOI 10.5281/zenodo.24250663). I declare no competing interests and no funding.")
para("Thank you for considering this submission.")
para("Yours sincerely,", after=18)
para("Leon Sandler")
os.makedirs(OUT, exist_ok=True)
out = os.path.join(OUT, "Cover_Letter_TBIO.docx")
from docx.oxml.ns import qn
_z = doc.settings.element.find(qn("w:zoom"))
if _z is not None and _z.get(qn("w:percent")) is None:
    _z.set(qn("w:percent"), "100")
doc.save(out)
print("saved", out)
