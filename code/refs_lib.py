# -*- coding: utf-8 -*-
"""Reference library built from the Crossref records fetched and checked by verify_refs.py (results/references_crossref.json).
Entries are formatted as: A. Author, B. Author, Title, Journal volume (year) pages. https://doi.org/..."""
import html
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.abspath(os.path.join(HERE, "..", "results"))


def _initials(given):
    parts = [p for p in re.split(r"[\s]+", given.replace(".", ". ").strip()) if p]
    out = []
    for p in parts:
        if "-" in p:
            out.append("-".join(q[0] + "." for q in p.split("-") if q))
        else:
            out.append(p[0] + ".")
    return "".join(out) if len(out) == 1 else " ".join(out).replace(". ", ".")


def _clean(t):
    t = re.sub(r"<[^>]+>", "", t)
    t = html.unescape(t)
    t = re.sub(r"\s+", " ", t).strip()
    return t.replace("‐", "-").replace("‑", "-")


def _title_case_fix(t):
    return t


def fmt(rec):
    au = rec["authors"]
    names = []
    for fam, giv in au[:6]:
        fam = fam.title() if fam.isupper() else fam
        names.append("%s %s" % (_initials(giv), fam) if giv else fam)
    s = ", ".join(names) + (", et al." if len(au) > 6 else "")
    title = _clean(rec["title"])
    jr = _clean(rec["journal"])
    vol = rec.get("volume")
    pg = rec.get("page") or rec.get("article_number")
    yr = rec.get("pub_year") or rec.get("year")
    ref = "%s, %s, %s" % (s, title, jr)
    if vol:
        ref += " %s" % vol
    ref += " (%s)" % yr
    if pg:
        ref += " %s" % str(pg).replace("-", "–")
    ref += ". https://doi.org/%s" % rec["doi"]
    return ref


def load():
    cache = json.load(open(os.path.join(RES, "references_crossref.json"), encoding="utf-8"))
    return {k: fmt(v) for k, v in cache.items() if v.get("status") == "ok"}


LIB = load()
