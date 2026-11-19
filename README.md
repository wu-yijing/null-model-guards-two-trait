# null-model-guards-two-trait

Reproducible, side-by-side recomputation of three cross-study null-model guards for two traits — polycystic ovary syndrome (PCOS, development case) and type 2 diabetes mellitus (T2D, replication case) — from raw GWAS Catalog lead associations. Companion code and results for a two-trait methods study on null-model guards for interpreting GWAS loci.

**Scope:** code and derived numerical results only. No manuscript text, no individual-level data.

## Guards (identical definitions applied to both traits)

- **G1 — power confounding.** Spearman rho between study case count and lead-locus count, computed over all associations and over genome-wide-significant (GWS, p <= 5e-8) associations only. A drop under the fixed threshold flags a threshold artefact.
- **G2 — threshold artefact.** Fraction of catalogued lead associations with p > 5e-8.
- **G3 — cross-study concordance.** Pairwise exact-rsID Jaccard across studies with >= 2 lead loci; reports mean Jaccard and the fraction of study pairs sharing zero rsIDs.

Both traits pass through identical parsing and guard definitions, so the case studies are directly comparable. No number is taken from any manuscript; all values are recomputed from the raw CSVs.

## Inputs (raw, unmodified)

- PCOS: `definitions.csv`, `lead_associations_qc.csv`
- T2D: `t2d_studies_qc.csv`, `t2d_lead_associations_qc.csv`

The unified script reads these from local paths set in `dual_trait_guard_recompute.py` (constants `PCOS` and `T2D`), which point at the `out/` directories produced by the two pipelines below. Adjust those two constants to your own output locations before running.

## Run

The two trait pipelines fetch public GWAS Catalog v2 data and write the QC CSVs that the unified recomputation consumes. Run in order:

```bash
# 1. Development case (PCOS): fetch + QC
python pcos/fetch_p1a.py
python pcos/analyse_p1a.py

# 2. Replication case (T2D): fetch + QC
python t2d/t2a_fetch.py
python t2d/t2a_analyse_local.py

# 3. Unified recomputation (identical guards on both traits)
python dual_trait_guard_recompute.py
```

Standard library only; no third-party dependencies and no credentials (public APIs). Each script uses paths relative to its own directory. The unified step writes `dual_trait_guards.json` and prints it to stdout. Fixed GWS threshold 5e-8. Record a stable commit before running so the analysis can be tied to a specific code version.

## Results at a glance

| Guard | PCOS (case 1) | T2D (case 2) |
|---|---|---|
| Lead associations used | 128 | 9,410 |
| Studies with usable leads | 18 | 235 |
| G1 rho (all -> GWS-only) | 0.766 -> 0.429 | 0.290 -> 0.441 |
| p > 5e-8 | 40.6% | 13.3% |
| G3 study pairs with zero shared rsID | 83.8% (88/105) | 82.9% (14,567/17,578) |

## Files

- `dual_trait_guard_recompute.py` — unified recomputation script (identical guards for both traits)
- `dual_trait_guards.json` — computed guard values for both traits
- `pcos/fetch_p1a.py`, `pcos/analyse_p1a.py` — PCOS (case 1) fetch and QC pipeline
- `t2d/t2a_fetch.py`, `t2d/t2a_analyse_local.py` — T2D (case 2) fetch and QC pipeline
- `README.md` — this file
- `LICENSE` — MIT

## License

MIT
