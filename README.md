# Muscle–glucose trade-offs in metabolic aging

Leon Sandler, Independent Researcher — sandler.leon@gmail.com — ORCID [0009-0007-4584-808X](https://orcid.org/0009-0007-4584-808X)

Code, results and manuscript sources for *"Muscle–glucose trade-offs in metabolic aging: a dynamical-systems model of an
interior optimal dose for glucose-lowering interventions"* (prepared for *Theory in Biosciences*). The paper is a theoretical
and computational study: **no human or animal data were used and no parameter was fitted.** Every number is a modeled outcome
conditional on the assumptions in Table 1 of the paper.

## The model

A latent information-loss variable `I`, a muscle reserve `M` and an energy-storage compartment `F` are coupled through
circulating glucose `G` (quasi-steady state, algebraic). A constant dose `u` from age 45 acts on glucose (directly), on
muscle (lean-mass loss) and on storage. The paper asks when an interior optimal dose exists, how robust that is, and what an
index `SIM = G·σ(I)/M` and its "metabolic age" can and cannot say.

```
G     = G0 (1 + aI I + aF F) / (1 + aM M + au u)
σ     = 1 + cσ I
dI/dt = [ r0 + rg G σ + rm (1-M)^pm + rf Φ(F) ] (1 - I)
dM/dt = kM [ (1 - λ I)(1 - φ H(u)) - M ],   H = u^nH / (ku^nH + u^nH)
dF/dt = gF (Fset(I) - F) - cF u F,           Fset = 0.4 + 0.4 I
```

## Layout

| Path | Content |
|---|---|
| `code/model.py`, `code/fastmodel.py` | the model; a vectorised NumPy RK4 reference (checked against LSODA) and a bit-identical compiled (numba) version used for the large ensembles; dose–response and classification |
| `code/analysis.py` | every analysis: nominal case, exact pathway decomposition, static criterion, ensembles (priors, horizon, adiposity, ablation, functional forms), one-at-a-time, phase maps, Sobol indices, identifiability |
| `code/verify.py` | executable checks of the implementation (Appendix B of the paper) → `results/verification.json` |
| `code/verify_refs.py`, `code/refs_lib.py` | every cited DOI is fetched from Crossref and checked (first author, year, title) |
| `code/figures.py` | all figures |
| `code/build_tbio*.py`, `tbio_part*.py`, `docx_engine.py`, `omml.py` | builds the manuscript (Word, native equations) from `results/*.json` |
| `results/` | generated results (JSON) |
| `figures/` | figures (PNG, 300 dpi) |
| `manuscript/` | the manuscript |

## Reproduce

```
pip install -r requirements.txt
cd code
python analysis.py all        # about 30-40 min on a desktop CPU; the ensembles dominate
python verify_refs.py         # needs internet (Crossref)
python verify.py
python figures.py
python build_tbio.py
```

Every ensemble is seeded from `SeedSequence([20261115, job code])`, so each number can be regenerated alone.

## Limits (as stated in the paper)

Mean-field and low-dimensional; parameters illustrative (they reproduce three qualitative targets of the untreated course and
are not fitted); information loss is non-decreasing, so the model can describe slowing but not reversal; SIM is a proposed,
unvalidated surrogate; nothing here is a dosing recommendation.

## License

MIT (see `LICENSE`).
