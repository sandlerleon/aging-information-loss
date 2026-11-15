# -*- coding: utf-8 -*-
"""Builds the manuscript for Theory in Biosciences from the parts tbio_part1.py (introduction, architecture, model),
tbio_part2.py (results) and tbio_part3.py (discussion, declarations, appendices). Every number in the results is read from results/*.json.

    python build_tbio.py   ->  ../manuscript/Aging_TBIO_Manuscript.docx (+ _no_author_block)
"""
import json
import os
import re
import sys

import numpy as np

import docx_engine as E
from refs_lib import LIB

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))
FIG = os.path.abspath(os.path.join(HERE, "..", "figures"))
OUT = os.path.abspath(os.path.join(HERE, "..", "manuscript"))


def L(n):
    return json.load(open(os.path.join(RES, n + ".json"), encoding="utf-8"))


NOM, DEC, STA, ENS, OAT, PHA, SOB, IDE, VER = (L(n) for n in ("nominal", "decomposition", "static", "ensemble", "oat", "phase", "sobol", "identifiability", "verification"))
if not all(v["passed"] for v in VER.values()):
    sys.exit("verification checks have not all passed; run verify.py")

TITLE = "Muscle–glucose trade-offs in metabolic aging: a dynamical-systems model of an interior optimal dose for glucose-lowering interventions"
B = []
add = B.append


def pct(x, d=1):
    return ("%." + str(d) + "f%%") % (100 * x)


def ci(d, k=1):
    return "%.*f–%.*f%%" % (k, 100 * d["ci95"][0], k, 100 * d["ci95"][1])


# ---------------------------------------------------------------------------------------------- front matter
FRONT = []
for part in ("tbio_part0.py",):
    exec(compile(open(os.path.join(HERE, part), encoding="utf-8").read(), part, "exec"), globals())
for part in ("tbio_part1.py", "tbio_part2.py", "tbio_part3.py"):
    exec(compile(open(os.path.join(HERE, part), encoding="utf-8").read(), part, "exec"), globals())

TOKENS = {"{ANSWERS}": ANSWERS, "{RHOMIN}": RHOMIN, "{FATFACTOR}": FATFACTOR}


def texts(b):
    k = b[0]
    if k in ("p", "p0", "fig"):
        return [b[1]] if k != "fig" else [b[2]]
    if k == "table":
        return [b[1]] + list(b[2]) + [c for r in b[3] for c in r]
    return []


def build(path, author=True):
    cited = []
    for b in B:
        for t in texts(b):
            for mm in re.finditer(r"\[@([^\]]+)\]", t):
                for kk in mm.group(1).split(";"):
                    cited.append(kk.strip().lstrip("@"))
    lib = {k: LIB[k] for k in dict.fromkeys(cited)}
    cit = E.Citations(lib)
    for b in B:
        for t in texts(b):
            cit.scan(t)
    doc = E.Doc(cit)
    for b in B:
        k = b[0]
        if k == "title":
            doc.title(b[1])
        elif k == "author":
            if author:
                doc.centered(b[1], 22)
        elif k == "small":
            if author:
                doc.centered(b[1], 19)
        elif k == "blank":
            doc.blank()
        elif k == "h1":
            doc.h1(b[1])
        elif k == "h2":
            doc.h2(b[1])
        elif k == "p":
            doc.p(b[1])
        elif k == "p0":
            doc.p(b[1], indent=False)
        elif k == "eq":
            doc.eq(b[1], b[2])
        elif k == "fig":
            doc.fig(os.path.join(FIG, b[1]), b[2], b[3])
        elif k == "table":
            doc.table(b[1], b[2], b[3], b[4] if len(b) > 4 else None)
        elif k == "refs":
            doc.h1("References")
            doc.refs(cit.entries())
    doc.save(path)
    return cit


for i, b in enumerate(B):                                    # substitute the numeric tokens
    if b[0] in ("p", "p0"):
        t = b[1]
        for k, v in TOKENS.items():
            t = t.replace(k, v)
        B[i] = (b[0], t)

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    c = build(os.path.join(OUT, "Aging_TBIO_Manuscript.docx"))
    build(os.path.join(OUT, "Aging_TBIO_Manuscript_no_author_block.docx"), author=False)
    print("built; %d references" % len(c.order))
