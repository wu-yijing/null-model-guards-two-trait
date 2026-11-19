"""T2D case study: QC + local guards (power confounding, cross-study concordance).
Mirrors p1a/analyse_p1a.py parsing exactly. Public data only. No credentials.

Run: python t2a_analyse_local.py
Outputs: t2a/out/t2d_lead_associations_qc.csv, t2d_studies_qc.csv,
         t2d_pairwise_concordance.csv, t2d_local_guards.json
"""
import csv
import json
import math
import re
import statistics
import urllib.parse
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA, OUT = ROOT / "data", ROOT / "out"
OUT.mkdir(exist_ok=True)

STUDIES = json.load(open(DATA / "catalog_t2d_studies.json", encoding="utf-8"))
ASSOC = json.load(open(DATA / "catalog_t2d_associations.json", encoding="utf-8"))
GW = 5e-8


# ---- parsers (verbatim from p1a/analyse_p1a.py) ----
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


def rsid_of(a):
    for r in a.get("riskAllele", []) or []:
        k = str(r.get("key", ""))
        if re.match(r"^rs\d+$", k):
            return k
    return None


def loc_of(a):
    for loc in a.get("locations", []) or []:
        m = re.match(r"\s*([0-9XY]+):(\d+)", str(loc))
        if m:
            return m.group(1), int(m.group(2))
    return None


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


def write_csv(name, rows, fields=None):
    if not rows:
        return
    fields = fields or list(rows[0].keys())
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    # meta per study
    meta = {}
    for s in STUDIES:
        acc = s["accessionId"]
        cases = parse_cases(s.get("initialSampleDescription"))
        ctrl = parse_controls(s.get("initialSampleDescription"))
        total = parse_total(s.get("discoverySampleAncestry"))
        derived = None
        if cases is None and total is not None and ctrl is not None:
            derived = total - ctrl
        meta[acc] = {
            "accession": acc,
            "reported_trait": s.get("reportedTrait") or "",
            "pubmed": s.get("pubmedId") or "",
            "pub_date": s.get("publicationDate") or "",
            "first_author": s.get("firstAuthor") or "",
            "journal": s.get("journal") or "",
            "association_count_reported": s.get("associationCount") or 0,
            "n_cases_reported": cases,
            "n_cases_derived": derived,
            "n_cases_used": cases if cases is not None else derived,
            "n_cases_source": "reported" if cases is not None else ("total_minus_controls" if derived is not None else "missing"),
            "n_controls": ctrl,
            "ancestry": parse_ancestry(s.get("discoverySampleAncestry")),
        }

    # QC lead associations
    lead_rows = []
    for acc, arr in ASSOC.items():
        for a in arr or []:
            rsid = rsid_of(a)
            p = parse_p(a)
            loc = loc_of(a)
            lead_rows.append({
                "accession": acc,
                "rsid": rsid or "",
                "p_value": f"{p:.6g}" if p is not None else "",
                "p_is_gws": (p is not None and p <= GW),
                "chrom": loc[0] if loc else "",
                "pos_grch38": loc[1] if loc else "",
                "reported_trait": meta[acc]["reported_trait"],
                "n_cases_used": meta[acc]["n_cases_used"],
                "n_cases_source": meta[acc]["n_cases_source"],
                "mapped_genes": ";".join(a.get("mappedGenes", []) or []),
            })
    write_csv("t2d_lead_associations_qc.csv", lead_rows)
    print(f"lead associations QC: {len(lead_rows)} rows", flush=True)

    # association sets = studies with >=1 usable lead association
    sets = defaultdict(list)
    for r in lead_rows:
        if r["rsid"] and r["p_value"]:
            sets[r["accession"]].append(r)

    set_rows = []
    for acc, rows in sets.items():
        rsids = set(r["rsid"] for r in rows)
        gws_rsids = set(r["rsid"] for r in rows if r["p_is_gws"] == "True" or r["p_is_gws"] is True)
        m = meta[acc]
        set_rows.append({
            "accession": acc,
            "reported_trait": m["reported_trait"],
            "n_cases_used": m["n_cases_used"],
            "n_cases_source": m["n_cases_source"],
            "ancestry": m["ancestry"],
            "pub_date": m["pub_date"],
            "n_lead_loci": len(rsids),
            "n_lead_loci_gws": len(gws_rsids),
            "n_assoc": len(rows),
            "n_gws_assoc": sum(1 for r in rows if r["p_is_gws"] == "True" or r["p_is_gws"] is True),
            "lead_rsids": ";".join(sorted(rsids)),
        })
    write_csv("t2d_studies_qc.csv", set_rows)
    print(f"association sets with usable leads: {len(set_rows)}", flush=True)

    # ---- Guard 1: power confounding ----
    # n_cases vs n_lead_loci, all sets vs GWS-only sets
    rho_all, n_all = spearman([r["n_cases_used"] for r in set_rows],
                              [r["n_lead_loci"] for r in set_rows])
    rho_gws, n_gws_sets = spearman([r["n_cases_used"] for r in set_rows if r["n_lead_loci_gws"] > 0],
                                   [r["n_lead_loci_gws"] for r in set_rows if r["n_lead_loci_gws"] > 0])
    # also: threshold artefact — fraction of associations not GWS
    n_withp = sum(1 for r in lead_rows if r["p_value"])
    n_nongws = sum(1 for r in lead_rows if r["p_value"] and not (r["p_is_gws"] == "True" or r["p_is_gws"] is True))

    # ---- Guard 2: cross-study concordance ----
    # pairwise exact rsID overlap among sets with >=2 loci
    big = [r for r in set_rows if r["n_lead_loci"] >= 2]
    pairwise = []
    for i in range(len(big)):
        for j in range(i + 1, len(big)):
            a, b = big[i], big[j]
            sa = set(a["lead_rsids"].split(";"))
            sb = set(b["lead_rsids"].split(";"))
            inter = sa & sb
            pairwise.append({
                "set_a": a["accession"], "set_b": b["accession"],
                "n_a": len(sa), "n_b": len(sb),
                "n_exact_shared": len(inter),
                "jaccard": round(len(inter) / len(sa | sb), 4) if (sa | sb) else 0,
            })
    write_csv("t2d_pairwise_concordance.csv", pairwise)
    mean_j = statistics.mean(p["jaccard"] for p in pairwise) if pairwise else None
    n_zero_shared = sum(1 for p in pairwise if p["n_exact_shared"] == 0)

    # select a representative locus set for downstream guards:
    # largest study by GWS loci, take up to 30 lead SNPs with GRCh38 coords
    cand = [r for r in set_rows if r["n_lead_loci_gws"] > 0]
    cand.sort(key=lambda r: (-r["n_lead_loci_gws"], -(r["n_cases_used"] or 0)))
    chosen = []
    seen = set()
    for r in cand:
        for row in sets[r["accession"]]:
            if row["rsid"] and row["chrom"] and row["pos_grch38"]:
                if row["rsid"] not in seen:
                    seen.add(row["rsid"])
                    rr = dict(row)
                    rr["source_accession"] = r["accession"]
                    chosen.append(rr)
        if len(chosen) >= 30:
            break
    chosen = chosen[:30]
    write_csv("t2d_representative_loci.csv", chosen)
    print(f"representative loci for orthogonal guards: {len(chosen)}", flush=True)

    guards = {
        "trait": "type 2 diabetes mellitus",
        "n_studies_total": len(STUDIES),
        "n_studies_with_usable_leads": len(set_rows),
        "n_lead_associations_qc": len(lead_rows),
        "n_assoc_with_p": n_withp,
        "n_assoc_nongws": n_nongws,
        "pct_assoc_nongws": round(100 * n_nongws / n_withp, 1) if n_withp else None,
        "guard1_power_confounding": {
            "rho_ncases_vs_nloci_all": round(rho_all, 3) if rho_all is not None else None,
            "n_pairs_all": n_all,
            "rho_ncases_vs_nloci_gws_only": round(rho_gws, 3) if rho_gws is not None else None,
            "n_pairs_gws": n_gws_sets,
            "n_sets_with_case_count": sum(1 for r in set_rows if r["n_cases_used"] is not None),
        },
        "guard2_cross_study_concordance": {
            "n_sets_with_ge2_loci": len(big),
            "n_pairwise_comparisons": len(pairwise),
            "mean_jaccard": round(mean_j, 4) if mean_j is not None else None,
            "n_pairs_zero_exact_shared": n_zero_shared,
            "pct_pairs_zero_shared": round(100 * n_zero_shared / len(pairwise), 1) if pairwise else None,
        },
        "representative_loci_n": len(chosen),
        "representative_source_accessions": sorted(set(r["source_accession"] for r in chosen)),
    }
    (OUT / "t2d_local_guards.json").write_text(json.dumps(guards, indent=1), encoding="utf-8")
    print("LOCAL_GUARDS", json.dumps(guards, indent=1), flush=True)


if __name__ == "__main__":
    main()
