# -*- coding: utf-8 -*-
"""Merge the safe part of the Rubriq language edit into the build without letting it change the science, the terms or the format.

Inputs
  base     the version that was sent to Rubriq      (manuscript/rubriq/Aging_TBIO_Manuscript_v1_base.docx)
  rubriq   the edited file that came back           (manuscript/rubriq/Aging_TBIO_Manuscript_v1_rubriq.docx)
  build    the current build from build_tbio.py     (the base for the output; it carries the equations, the current text and the references)

Every word-level difference between base and rubriq is classified. It is applied to the build only if it is
  * a comma or semicolon added after a word (punctuation only), or
  * one of the whitelisted spelling / wording substitutions (SAFE below: US spelling, 'rises' -> 'increases', ...),
and only if the word is found in the build with the same neighbouring words. Everything else is rejected: article and
preposition changes, wording that changes a term ('dynamical-systems', 'stationarity', 'quasi-steady', 'co-intervention'),
reference-list entries, titles, numbers and DOIs. A report of the decisions is written next to the output.

    python merge_rubriq_tbio.py <build.docx> <out.docx> [report.txt]
"""
import difflib
import os
import re
import sys

import docx

HERE = os.path.dirname(os.path.abspath(__file__))
RUB = os.path.abspath(os.path.join(HERE, "..", "manuscript", "rubriq"))
BASE = os.path.join(RUB, "Aging_TBIO_Manuscript_v1_base.docx")
EDIT = os.path.join(RUB, "Aging_TBIO_Manuscript_v1_rubriq.docx")

SAFE = {"rises": "increases", "rise": "increase", "falls": "decreases", "fall": "decrease", "about": "approximately", "largest": "greatest",
        "larger": "greater", "grey": "gray", "Grey": "Gray", "Grey:": "Gray:", "non-negativity": "nonnegativity", "non-negative": "nonnegative",
        "non-decreasing": "nondecreasing", "non-linear": "nonlinear", "non-monotonically": "nonmonotonically", "pseudo-random": "pseudorandom",
        "log-derivative": "log derivative", "grows": "increases", "means": "indicates"}
PUNCT = ",;:.)(]["


def split_tok(t):
    m = re.match(r"^(.*?)([,;:.)\]]*)$", t)
    return m.group(1), m.group(2)


def classify(old, new):
    """Return a replacement for token `old` (or None if the edit is not safe)."""
    co, po = split_tok(old)
    cn, pn = split_tok(new)
    if co == cn and pn in (po + ",", po + ";") and po == "":
        return new                                             # comma / semicolon added after the word
    if co in SAFE and cn == SAFE[co] and po == pn:
        return new
    if old in SAFE and new == SAFE[old]:
        return new
    return None


def body(d):
    return [p for p in d.paragraphs]


def edits(base_p, rub_p):
    A, R = base_p.text.split(), rub_p.text.split()
    sm = difflib.SequenceMatcher(None, A, R, autojunk=False)
    out = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            continue
        old, new = A[i1:i2], R[j1:j2]
        prev = A[i1 - 1] if i1 > 0 else None
        nxt = A[i2] if i2 < len(A) else None
        rep = classify(old[0], new[0]) if (tag == "replace" and len(old) == 1 and len(new) == 1) else None
        out.append(dict(old=old, new=new, prev=prev, next=nxt, safe=rep))
    return out


def apply(par, op):
    """Replace the single word `old` in paragraph `par` (found with the same neighbours) inside one run. Returns True if applied."""
    runs = par.runs
    S = "".join(r.text for r in runs)
    spans = [(m.start(), m.end()) for m in re.finditer(r"\S+", S)]
    words = [S[a:b] for a, b in spans]
    old, rep = op["old"][0], op["safe"]
    for k, w in enumerate(words):
        if w != old:
            continue
        if op["prev"] is not None and (k == 0 or words[k - 1] != op["prev"]):
            continue
        if op["next"] is not None and (k + 1 >= len(words) or words[k + 1] != op["next"]):
            continue
        a, b = spans[k]
        pos = 0
        for r in runs:
            lo, hi = pos, pos + len(r.text)
            if lo <= a and b <= hi:
                r.text = r.text[:a - lo] + rep + r.text[b - lo:]
                return True
            pos = hi
        return False                                           # the word spans several runs (formatting): skip
    return False


def main(build, out, report):
    base, rub, new = docx.Document(BASE), docx.Document(EDIT), docx.Document(build)
    bp, rp, npar = body(base), body(rub), body(new)
    assert len(bp) == len(rp), "the Rubriq file does not align with its base"
    ntext = [p.text for p in npar]
    used = set()
    n_acc = n_rej = n_missing = 0
    log = []
    for i, (a, r) in enumerate(zip(bp, rp)):
        if a.text.strip() == r.text.strip():
            continue
        if re.match(r"^\[\d+\]", a.text.strip()):
            log.append("REJECT whole paragraph (reference-list entry): %s ..." % a.text[:70])
            n_rej += len(edits(a, r))
            continue
        ops = edits(a, r)
        # locate the paragraph in the build
        best, bj = 0.0, None
        for j, t in enumerate(ntext):
            if j in used:
                continue
            if abs(len(t) - len(a.text)) > 0.5 * max(len(a.text), 1) + 40:
                continue
            ratio = difflib.SequenceMatcher(None, a.text, t, autojunk=False).quick_ratio()
            if ratio > best:
                best, bj = ratio, j
        if bj is None or difflib.SequenceMatcher(None, a.text, ntext[bj], autojunk=False).ratio() < 0.75:
            log.append("SKIP paragraph not found in the build: %s ..." % a.text[:70])
            continue
        used.add(bj)
        for op in ops:
            tag = "%s -> %s" % (" ".join(op["old"]), " ".join(op["new"]))
            if op["safe"] is None:
                n_rej += 1
                log.append("  reject  [%s]" % tag)
                continue
            if apply(npar[bj], op):
                n_acc += 1
                log.append("  ACCEPT  [%s]  (%s)" % (tag, a.text[:40]))
            else:
                n_missing += 1
                log.append("  absent  [%s]  not found in the build with the same neighbours" % tag)
    new.save(out)
    summary = "accepted %d, rejected %d, not found in the build %d" % (n_acc, n_rej, n_missing)
    if report:
        open(report, "w", encoding="utf-8").write(summary + "\n\n" + "\n".join(log))
    print(summary)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
