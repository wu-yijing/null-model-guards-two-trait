"""Unified, side-by-side recomputation of the three cross-study guards for
BOTH case-study traits (PCOS and type 2 diabetes) from raw deposited CSVs.

The point of this script is methodological: both traits are passed through
IDENTICAL guard definitions and IDENTICAL parsing, so the two case studies are
comparable. It does not reuse any number printed in either manuscript.

Guards (applied identically to both traits):
  G1 power confounding   : Spearman(n_cases, n_lead_loci) over studies; all sets
                           vs GWS-only (p<=5e-8) loci.
  G2 threshold artefact  : fraction of catalogued lead associations with p>5e-8.
  G3 cross-study concord : pairwise exact-rsID Jaccard over studies with >=2 loci;
                           mean Jaccard and fraction of pairs with zero shared rsID.

Inputs (raw, unmodified):
  PCOS : .../2026-11-19-16-25-03/p1a/out/{definitions.csv,lead_associations_qc.csv}
  T2D  : .../2026-11-19-18-12-17/t2a/out/{t2d_studies_qc.csv,t2d_lead_associations_qc.csv}

Output: dual_trait_guards.json (next to this script).
"""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

GW = 5e-8
PCOS = Path("E:/workbuddy/2026-11-19-16-25-03/p1a/out")
T2D = Path("E:/workbuddy/2026-11-19-18-12-17/t2a/out")
HERE = Path(__file__).resolve().parent


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def parse_p(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def truthy(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def spearman(x, y):
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


def run_trait(name, assoc_rows, ncase_by_acc, primary_only):
    """assoc_rows: list of dict with accession, rsid, p_value (and optionally
    in_primary_analysis). ncase_by_acc: dict accession -> n_cases or None."""
    if primary_only:
        assoc_rows = [r for r in assoc_rows if truthy(r.get("in_primary_analysis"))]

    # G2 threshold artefact
    ps = [parse_p(r.get("p_value")) for r in assoc_rows]
    ps = [p for p in ps if p is not None]
    n_nongws = sum(1 for p in ps if p > GW)
    g2 = {
        "n_assoc_with_p": len(ps),
        "n_p_gt_5e8": n_nongws,
        "pct_p_gt_5e8": round(100 * n_nongws / len(ps), 1) if ps else None,
    }

    # per-study lead loci (all and GWS-only)
    loci_all = defaultdict(set)
    loci_gws = defaultdict(set)
    for r in assoc_rows:
        rsid = (r.get("rsid") or "").strip()
        p = parse_p(r.get("p_value"))
        if not rsid or p is None:
            continue
        loci_all[r["accession"]].add(rsid)
        if p <= GW:
            loci_gws[r["accession"]].add(rsid)

    # G1 power confounding
    accs = [a for a in loci_all if ncase_by_acc.get(a) is not None]
    rho_all, n_all = spearman([ncase_by_acc[a] for a in accs],
                              [len(loci_all[a]) for a in accs])
    accs_g = [a for a in loci_gws if loci_gws[a] and ncase_by_acc.get(a) is not None]
    rho_gws, n_gws = spearman([ncase_by_acc[a] for a in accs_g],
                              [len(loci_gws[a]) for a in accs_g])
    g1 = {
        "rho_all": round(rho_all, 3) if rho_all is not None else None,
        "n_studies_all": n_all,
        "rho_gws_only": round(rho_gws, 3) if rho_gws is not None else None,
        "n_studies_gws": n_gws,
        "n_studies_with_case_count": sum(
            1 for a in set(loci_all) | set(loci_gws)
            if ncase_by_acc.get(a) is not None
        ),
    }

    # G3 cross-study concordance (studies with >=2 lead loci)
    big = [a for a in loci_all if len(loci_all[a]) >= 2]
    js, zero = [], 0
    for i in range(len(big)):
        for j in range(i + 1, len(big)):
            sa, sb = loci_all[big[i]], loci_all[big[j]]
            inter = sa & sb
            if not inter:
                zero += 1
            js.append(len(inter) / len(sa | sb) if (sa | sb) else 0.0)
    npairs = len(js)
    g3 = {
        "n_studies_ge2_loci": len(big),
        "n_pairwise": npairs,
        "mean_jaccard": round(statistics.mean(js), 4) if js else None,
        "n_pairs_zero_shared": zero,
        "pct_pairs_zero_shared": round(100 * zero / npairs, 1) if npairs else None,
    }

    return {
        "trait": name,
        "n_lead_assoc_used": len(assoc_rows),
        "n_studies_with_leads": len(set(loci_all) | set(loci_gws)),
        "G1_power_confounding": g1,
        "G2_threshold_artefact": g2,
        "G3_cross_study_concordance": g3,
    }


def main():
    # ---- PCOS ----
    pdef = read_csv(PCOS / "definitions.csv")
    ncase_p = {}
    primary_accs = set()
    for r in pdef:
        ncase_p[r["accession"]] = (int(r["n_cases_used"])
                                   if r.get("n_cases_used", "").strip() else None)
        if truthy(r.get("in_primary_analysis")):
            primary_accs.add(r["accession"])
    passoc = read_csv(PCOS / "lead_associations_qc.csv")
    # restrict to primary definitions for the case-count map, keep association
    # primary flag as deposited
    pcos = run_trait("polycystic ovary syndrome", passoc, ncase_p, primary_only=True)

    # ---- T2D ----
    tstud = read_csv(T2D / "t2d_studies_qc.csv")
    ncase_t = {}
    for r in tstud:
        ncase_t[r["accession"]] = (int(r["n_cases_used"])
                                   if r.get("n_cases_used", "").strip() else None)
    tassoc = read_csv(T2D / "t2d_lead_associations_qc.csv")
    t2d = run_trait("type 2 diabetes mellitus", tassoc, ncase_t, primary_only=False)

    out = {
        "note": "Both traits computed with identical guard definitions from raw deposited CSVs; no manuscript numbers reused.",
        "gws_threshold": GW,
        "PCOS_case1": pcos,
        "T2D_case2": t2d,
    }
    (HERE / "dual_trait_guards.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
