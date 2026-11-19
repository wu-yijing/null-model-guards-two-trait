"""P1-A minimal analysis: PCOS phenotype-definition sensitivity.

Uses only GWAS Catalog v2 public lead associations (no summary-statistic
download). Pure Python standard library, no numpy/pandas/ggplot required.

Run: python analyse_p1a.py
Outputs: out/*.csv, out/summary.json
"""
import csv
import json
import math
import re
import statistics
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA, OUT = ROOT / "data", ROOT / "out"
OUT.mkdir(exist_ok=True)

STUDIES = json.loads((DATA / "catalog_pcos_studies.json").read_text(encoding="utf-8"))["_embedded"]["studies"]
ASSOC = json.loads((DATA / "catalog_associations.json").read_text(encoding="utf-8"))
SIZES = json.loads((DATA / "sumstats_sizes.json").read_text(encoding="utf-8")) if (DATA / "sumstats_sizes.json").exists() else {}

# Definition classes. Cross-trait pleiotropy studies are excluded from the
# primary definition-sensitivity analysis (their effect sizes are not PCOS main effects).
DEFN = {
    "GCST007089": "clinical_criteria_meta",
    "GCST001634": "clinical_criteria", "GCST002741": "clinical_criteria",
    "GCST000914": "clinical_criteria", "GCST005043": "clinical_criteria",
    "GCST001669": "clinical_criteria", "GCST006659": "clinical_criteria",
    "GCST005273": "clinical_criteria", "GCST003144": "clinical_criteria",
    "GCST010687": "clinical_criteria_multi_ancestry", "GCST90693150": "clinical_criteria",
    "GCST90483500": "clinical_criteria", "GCST90727246": "clinical_criteria_south_asian",
    "GCST90454205": "icd10_E28.2",
    "GCST90435733": "phecode_256.4", "GCST90651355": "phecode_256.4",
    "GCST90477347": "phecode_256.4", "GCST90481599": "phecode_256.4",
    "GCST90481600": "phecode_256.4", "GCST90480881": "phecode_256.4",
    "GCST90077795": "selfreport_UKB", "GCST90081781": "gene_based_UKB",
    "GCST90044902": "age_adjusted", "GCST90044903": "age_BMI_adjusted",
    "GCST90558301": "lean_BMI_le25", "GCST90558302": "obese_BMI_gt25",
    "GCST010570": "subtype_reproductive", "GCST010569": "subtype_metabolic",
    "GCST010568": "subtype_indeterminate",
}
for acc in ["GCST90134404", "GCST90134405", "GCST90134406", "GCST90134407", "GCST90134408"]:
    DEFN[acc] = "cross_trait_pleiotropy"

PAIR_GROUP = {
    "GCST90558301": "BMI strata (PMID 38408933)",
    "GCST90558302": "BMI strata (PMID 38408933)",
    "GCST010568": "Clinical subtypes (PMID 32574161)",
    "GCST010569": "Clinical subtypes (PMID 32574161)",
    "GCST010570": "Clinical subtypes (PMID 32574161)",
    "GCST90044902": "Covariate adjustment (PMID 34791234)",
    "GCST90044903": "Covariate adjustment (PMID 34791234)",
    "GCST007089": "Definition family (clinical vs registry)",
    "GCST90454205": "Definition family (clinical vs registry)",
    "GCST90435733": "Definition family (clinical vs registry)",
}

EXCLUDE = {"cross_trait_pleiotropy", "gene_based_UKB"}


def parse_beta(raw):
    """Return (signed magnitude, direction) or None."""
    if not raw:
        return None
    m = re.search(r"([0-9]*\.?[0-9]+)\s*unit\s*(increase|decrease)", str(raw), re.I)
    if not m:
        return None
    val = float(m.group(1))
    return (-val if m.group(2).lower() == "decrease" else val), m.group(2).lower()


def parse_ci(raw):
    """Return (lo, hi, SE) or None. Handles [a-b], [a-b] nbsp, en dash, bare a-b."""
    if not raw:
        return None
    s = str(raw).replace("\xa0", " ").replace("–", "-").replace("−", "-")
    m = re.search(r"(-?\d+\.?\d*)\s*-\s*(-?\d+\.?\d*)", s)
    if not m:
        return None
    lo, hi = float(m.group(1)), float(m.group(2))
    if hi < lo:
        lo, hi = hi, lo
    return lo, hi, (hi - lo) / (2 * 1.96)


def parse_p(a):
    try:
        return float(a["pValue"]) * 10 ** int(a["pValueExponent"])
    except (KeyError, TypeError, ValueError):
        return None


def parse_cases(desc):
    for d in desc or []:
        m = re.search(r"(\d[\d,]*)[^,]*?\bcases\b", d, re.I)
        if m:
            return int(m.group(1).replace(",", ""))
    return None


def parse_controls(desc):
    for d in desc or []:
        m = re.search(r"(\d[\d,]*)[^,]*?\bcontrols\b", d, re.I)
        if m:
            return int(m.group(1).replace(",", ""))
    return None


def parse_total(ancestry):
    vals = []
    for a in ancestry or []:
        m = re.match(r"\s*(\d[\d,]*)", a)
        if m:
            vals.append(int(m.group(1).replace(",", "")))
    return sum(vals) if vals else None


def parse_ancestry(ancestry):
    return "; ".join(ancestry or []) or "NR"


def parse_locus(loc):
    m = re.match(r"\s*([0-9XY]+):(\d+)", str(loc or ""))
    return (m.group(1), int(m.group(2))) if m else None


def spearman(x, y):
    """Spearman rho with average ranks; returns None if <3 paired points."""
    pairs = [(a, b) for a, b in zip(x, y) if a is not None and b is not None]
    if len(pairs) < 3:
        return None, len(pairs)

    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs])
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return (num / den if den else None), len(pairs)


def write_csv(name, rows, fields):
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {name}: {len(rows)} rows", flush=True)


# ---------- 1. definition catalogue ----------
def_rows = []
for s in STUDIES:
    acc = s["accessionId"]
    cls = DEFN.get(acc, "unclassified")
    ss = s.get("summaryStatistics")
    has_ss = bool(ss) and str(ss).startswith("http")
    main_files = SIZES.get(acc, {})
    nbytes = max([v for v in main_files.values() if v], default=None)
    n_total = parse_total(s.get("discoverySampleAncestry"))
    n_cases = parse_cases(s.get("initialSampleDescription"))
    n_ctrl = parse_controls(s.get("initialSampleDescription"))
    n_derived = (n_total - n_ctrl) if (n_total and n_ctrl and n_total > n_ctrl) else None
    def_rows.append({
        "accession": acc, "definition_class": cls,
        "reported_trait": s.get("reportedTrait"), "title": s.get("title"),
        "pubmed": s.get("pubmedId"), "year": str(s.get("publicationDate", ""))[:4],
        "first_author": s.get("firstAuthor"), "journal": s.get("journal"),
        "n_total_discovery": n_total,
        "n_cases_reported": n_cases,
        "n_controls_reported": n_ctrl,
        "n_cases_derived": n_derived,
        "n_cases_used": n_cases if n_cases is not None else n_derived,
        "n_cases_source": "reported" if n_cases is not None else ("total_minus_controls" if n_derived else "missing"),
        "ancestry": parse_ancestry(s.get("discoverySampleAncestry")),
        "n_catalog_associations": s.get("associationCount") or 0,
        "sumstats_available": has_ss,
        "sumstats_main_file_bytes": nbytes,
        "sumstats_files": "; ".join((json.loads((DATA / "sumstats_directory.json").read_text(encoding="utf-8")).get(acc, []))[:3]),
        "in_primary_analysis": cls not in EXCLUDE and (s.get("associationCount") or 0) > 0,
    })
write_csv("definitions.csv", def_rows, list(def_rows[0].keys()))
META = {r["accession"]: r for r in def_rows}

# ---------- 2. lead association QC ----------
qc_rows = []
for acc, alist in ASSOC.items():
    for a in alist:
        beta = parse_beta(a.get("beta"))
        ci = parse_ci(a.get("ci"))
        p = parse_p(a)
        ra = (a.get("riskAllele") or [{}])[0]
        loc = (a.get("locations") or [None])[0]
        chrom, pos = parse_locus(loc) if loc else (None, None)
        qc_rows.append({
            "accession": acc,
            "definition_class": DEFN.get(acc, "unclassified"),
            "in_primary_analysis": META[acc]["in_primary_analysis"],
            "rsid": ra.get("key"),
            "risk_allele": ra.get("label"),
            "risk_frequency": a.get("riskFrequency"),
            "chrom": chrom, "pos": pos,
            "mapped_genes": "; ".join(a.get("mappedGenes") or []),
            "beta_raw": a.get("beta"), "beta_signed": beta[0] if beta else None,
            "beta_direction": beta[1] if beta else None,
            "ci_raw": a.get("ci"), "ci_lo": ci[0] if ci else None,
            "ci_hi": ci[1] if ci else None, "se_from_ci": ci[2] if ci else None,
            "p_value": p, "neglog10p": (-math.log10(p) if p and p > 0 else None),
            "trait_name": "; ".join(a.get("traitName") or []),
        })
write_csv("lead_associations_qc.csv", qc_rows, list(qc_rows[0].keys()))

PRIMARY = [r for r in qc_rows if r["in_primary_analysis"]]
BY_STUDY = {}
for r in PRIMARY:
    BY_STUDY.setdefault(r["accession"], []).append(r)

# ---------- 3. pairwise definition concordance ----------
def interval_concordance(A, B, window=250000):
    """Fraction of loci in the smaller set with a partner within +/-window."""
    if not A or not B:
        return None, 0
    la = [(r["chrom"], r["pos"]) for r in A if r["chrom"] and r["pos"]]
    lb = [(r["chrom"], r["pos"]) for r in B if r["chrom"] and r["pos"]]
    if not la or not lb:
        return None, 0
    matched_a = sum(1 for c, p in la if any(c == c2 and abs(p - p2) <= window for c2, p2 in lb))
    matched_b = sum(1 for c, p in lb if any(c == c2 and abs(p - p2) <= window for c2, p2 in la))
    return (matched_a + matched_b) / (len(la) + len(lb)), matched_a + matched_b


ov_rows = []
for a, b in combinations(sorted(BY_STUDY), 2):
    A, B = BY_STUDY[a], BY_STUDY[b]
    ra = {r["rsid"] for r in A if r["rsid"]}
    rb = {r["rsid"] for r in B if r["rsid"]}
    ga = {g for r in A for g in r["mapped_genes"].split("; ") if g}
    gb = {g for r in B for g in r["mapped_genes"].split("; ") if g}
    shared = ra & rb
    betas_a = {r["rsid"]: r["beta_signed"] for r in A if r["rsid"] and r["beta_signed"] is not None}
    betas_b = {r["rsid"]: r["beta_signed"] for r in B if r["rsid"] and r["beta_signed"] is not None}
    shared_beta = set(betas_a) & set(betas_b)
    concordant = sum(1 for rs in shared_beta if (betas_a[rs] > 0) == (betas_b[rs] > 0))
    ic, _ = interval_concordance(A, B)
    ov_rows.append({
        "study_a": a, "study_b": b,
        "class_a": DEFN.get(a), "class_b": DEFN.get(b),
        "n_loci_a": len(A), "n_loci_b": len(B),
        "shared_rsid": len(shared),
        "jaccard_rsid": round(len(shared) / len(ra | rb), 3) if (ra | rb) else None,
        "shared_genes": len(ga & gb),
        "jaccard_genes": round(len(ga & gb) / len(ga | gb), 3) if (ga | gb) else None,
        "interval_concordance_250kb": round(ic, 3) if ic is not None else None,
        "n_shared_with_beta": len(shared_beta),
        "direction_concordance": round(concordant / len(shared_beta), 3) if shared_beta else None,
        "same_pmid": META[a]["pubmed"] == META[b]["pubmed"],
    })
write_csv("definition_overlap.csv", ov_rows, list(ov_rows[0].keys()))

# ---------- 4. paired within-study contrasts ----------
pair_rows = []
for acc, rows in sorted(BY_STUDY.items()):
    grp = PAIR_GROUP.get(acc)
    if not grp:
        continue
    betas = [abs(r["beta_signed"]) for r in rows if r["beta_signed"] is not None]
    ps = [r["p_value"] for r in rows if r["p_value"]]
    pair_rows.append({
        "paired_group": grp, "accession": acc, "definition_class": DEFN.get(acc),
        "reported_trait": META[acc]["reported_trait"], "n_cases": META[acc]["n_cases_used"],
        "n_cases_source": META[acc]["n_cases_source"],
        "n_total": META[acc]["n_total_discovery"],
        "n_lead_loci": len(rows), "n_with_beta": len(betas),
        "median_abs_beta": round(statistics.median(betas), 4) if betas else None,
        "min_p": min(ps) if ps else None,
        "max_neglog10p": round(max(-math.log10(p) for p in ps), 2) if ps else None,
        "lead_rsids": "; ".join(str(r["rsid"]) for r in rows[:8]),
    })
write_csv("paired_contrasts.csv", pair_rows, list(pair_rows[0].keys()))

# ---------- 5. specification curve data ----------
spec_rows = []
for acc, rows in sorted(BY_STUDY.items()):
    m = META[acc]
    betas = [abs(r["beta_signed"]) for r in rows if r["beta_signed"] is not None]
    ps = [r["p_value"] for r in rows if r["p_value"]]
    spec_rows.append({
        "accession": acc, "definition_class": DEFN.get(acc),
        "label": f"{acc} {DEFN.get(acc)}",
        "pubmed": m["pubmed"], "year": m["year"],
        "n_cases": m["n_cases_used"], "n_cases_source": m["n_cases_source"],
        "n_total": m["n_total_discovery"],
        "ancestry": m["ancestry"], "n_lead_loci": len(rows),
        "n_with_beta": len(betas),
        "median_abs_beta": round(statistics.median(betas), 4) if betas else None,
        "max_neglog10p": round(max(-math.log10(p) for p in ps), 2) if ps else None,
        "sumstats_available": m["sumstats_available"],
        "paired_group": PAIR_GROUP.get(acc, "none"),
    })
spec_rows.sort(key=lambda r: (r["definition_class"], -(r["n_cases"] or 0)))
write_csv("specification_curve.csv", spec_rows, list(spec_rows[0].keys()))

# Power confound guard: are locus counts driven by sample size rather than definition?
rho, npairs = spearman([r["n_cases"] for r in spec_rows], [r["n_lead_loci"] for r in spec_rows])

# ---------- 6. literature clinical anchors (NOT computed in this run) ----------
anchors = [
    {"outcome": "Incident CVD", "source": "Zhao 2024/25 UK Biobank", "pmid": "39563606",
     "definition": "UKB PCOS (645 cases, matched cohort)", "effect": 1.77, "lo": 1.19, "hi": 2.65,
     "n_cases": 645, "note": "Literature-extracted; Cox HR"},
    {"outcome": "IHD/cerebrovascular events", "source": "BJOG 2024 abstract UK Biobank", "pmid": "",
     "definition": "UKB PCOS (15,747 cases, 75,142 women)", "effect": 0.90, "lo": 0.79, "hi": 1.02,
     "n_cases": 15747, "note": "Literature-extracted; Cox HR"},
    {"outcome": "All-cause dementia", "source": "Henney 2025 UK Biobank", "pmid": "",
     "definition": "UKB PCOS (1,008 cases, BMI-matched)", "effect": 0.96, "lo": 0.32, "hi": 2.89,
     "n_cases": 1008, "note": "Literature-extracted; HR, wide CI"},
    {"outcome": "Dementia prevalence", "source": "Vine 2024 Canada cohort", "pmid": "",
     "definition": "Administrative PCOS (16,531 cases)", "effect": 2.19, "lo": None, "hi": None,
     "n_cases": 16531, "note": "Literature-extracted; crude prevalence ratio 0.48%/0.22%, NOT an adjusted HR"},
]
write_csv("clinical_anchors.csv", anchors, list(anchors[0].keys()))

# ---------- 7. run summary ----------
summary = {
    "run_date": "2026-11-11",
    "data_source": "GWAS Catalog v2 public API (studies + associations + FTP directory listings)",
    "n_studies": len(STUDIES),
    "n_studies_with_associations": len(BY_STUDY),
    "n_lead_associations_total": len(qc_rows),
    "n_lead_associations_primary": len(PRIMARY),
    "n_with_parsable_beta": sum(1 for r in PRIMARY if r["beta_signed"] is not None),
    "n_with_parsable_ci_se": sum(1 for r in PRIMARY if r["se_from_ci"] is not None),
    "n_sumstats_directories": sum(1 for s in STUDIES if s.get("summaryStatistics") and str(s["summaryStatistics"]).startswith("http")),
    "sumstats_main_bytes_min": min([v for d in SIZES.values() for v in d.values() if v], default=None),
    "sumstats_main_bytes_max": max([v for d in SIZES.values() for v in d.values() if v], default=None),
    "spearman_ncases_vs_nloci": round(rho, 3) if rho is not None else None,
    "spearman_n_pairs": npairs,
    "n_specs_with_case_count": sum(1 for r in spec_rows if r["n_cases"] is not None),
    "n_case_counts_derived": sum(1 for r in spec_rows if r["n_cases_source"] == "total_minus_controls"),
    "strict_ldsc_feasible_this_run": False,
    "strict_ldsc_blockers": [
        "GCST007089 ships only PCOS_top_10000_SNPS.txt (top-10k subset, not genome-wide); LDSC on selected SNPs is biased",
        "Full genome-wide sumstats are 4.6 MB - 2.1 GB; EBI FTP throughput on this machine is known to be very slow",
        "No 1000G LD-score reference (.l2.ldscore.gz) present locally",
        "Python 3.13 has no numpy/scipy/pandas; ldsc cannot run without them",
        "IEU OpenGWAS token expired (HTTP 401 verified this run), so alternate API access is unavailable",
    ],
    "novelty_caveat": "GCST007089 (Day 2018) is titled '...suggests shared genetic architecture for different diagnosis criteria': the qualitative claim that criteria share genetic architecture is NOT new. P1-A's contribution is quantitative specification-curve drift across biobank definitions, subtypes, BMI strata and clinical outcomes.",
    "outputs": ["definitions.csv", "lead_associations_qc.csv", "definition_overlap.csv",
                "paired_contrasts.csv", "specification_curve.csv", "clinical_anchors.csv"],
}
# ---------- 8. Stage 2 LDSC input inventory (real URLs + byte sizes) ----------
DIRLIST = json.loads((DATA / "sumstats_directory.json").read_text(encoding="utf-8")) if (DATA / "sumstats_directory.json").exists() else {}
ss_rows = []
for s in STUDIES:
    link = s.get("summaryStatistics")
    if not link or not str(link).startswith("http"):
        continue
    acc = s["accessionId"]
    base = str(link).rstrip("/") + "/"
    if base.startswith("http://"):
        base = "https://" + base[7:]
    files = DIRLIST.get(acc, [])
    has_harm = "harmonised/" in files
    main_files = SIZES.get(acc, {})
    if main_files:
        for fname, nbytes in main_files.items():
            ss_rows.append({
                "accession": acc, "definition_class": DEFN.get(acc, "unclassified"),
                "reported_trait": s.get("reportedTrait"), "pubmed": s.get("pubmedId"),
                "ancestry": parse_ancestry(s.get("discoverySampleAncestry")),
                "n_cases_used": META[acc]["n_cases_used"],
                "main_file": fname, "bytes": nbytes,
                "download_url": base + fname,
                "has_harmonised_dir": has_harm,
                "ldsc_eligible_candidate": "yes (check build, ancestry, SNP coverage)",
            })
    else:
        ss_rows.append({
            "accession": acc, "definition_class": DEFN.get(acc, "unclassified"),
            "reported_trait": s.get("reportedTrait"), "pubmed": s.get("pubmedId"),
            "ancestry": parse_ancestry(s.get("discoverySampleAncestry")),
            "n_cases_used": META[acc]["n_cases_used"],
            "main_file": "; ".join(files[:3]), "bytes": None,
            "download_url": base, "has_harmonised_dir": has_harm,
            "ldsc_eligible_candidate": "NO - top-10k SNP subset only, not genome-wide",
        })
write_csv("ldsc_stage2_inputs.csv", ss_rows, list(ss_rows[0].keys()))

(OUT / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False), encoding="utf-8")
print(json.dumps(summary, indent=1, ensure_ascii=False), flush=True)
