"""P1-A minimal: fetch PCOS GWAS Catalog v2 metadata, lead associations and
public summary-statistics directory listings. No credentials required.

Run: python fetch_p1a.py
Outputs: data/*.json
"""
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
API = "https://www.ebi.ac.uk/gwas/api/v2"
TRAIT = "polycystic ovary syndrome"
UA = {"Accept": "application/json", "User-Agent": "P1A-repro/1.0"}


def get_json(url, timeout=25):
    req = urllib.request.Request(url, headers=UA)
    return json.load(urllib.request.urlopen(req, timeout=timeout))


def main():
    # 1. studies
    url = f"{API}/studies?efoTrait={urllib.parse.quote(TRAIT)}&size=100"
    studies = get_json(url)
    (DATA / "catalog_pcos_studies.json").write_text(
        json.dumps(studies, indent=1, ensure_ascii=False), encoding="utf-8")
    st = studies["_embedded"]["studies"]
    print(f"studies: {len(st)}", flush=True)

    # 2. lead associations for every study with >=1 catalogue association
    assoc = {}
    for s in st:
        acc = s["accessionId"]
        if not s.get("associationCount"):
            assoc[acc] = []
            continue
        try:
            d = get_json(f"{API}/studies/{acc}/associations?size=1000")
            assoc[acc] = d["_embedded"]["associations"]
            print(f"  {acc}: {len(assoc[acc])} associations", flush=True)
        except Exception as e:
            assoc[acc] = []
            print(f"  {acc}: FAILED {e}", flush=True)
        time.sleep(0.25)
    (DATA / "catalog_associations.json").write_text(
        json.dumps(assoc, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"total associations: {sum(len(v) for v in assoc.values())}", flush=True)

    # 3. summary-statistics directory listing + file sizes (HEAD only, no download)
    dirs, sizes = {}, {}
    for s in st:
        link = s.get("summaryStatistics")
        if not link or not str(link).startswith("http"):
            continue
        base = str(link).rstrip("/") + "/"
        if base.startswith("http://"):
            base = "https://" + base[7:]
        acc = s["accessionId"]
        try:
            html = urllib.request.urlopen(base, timeout=10).read().decode("utf-8", "ignore")
            files = [f for f in re.findall(r'href="([^"?][^"]*)"', html)
                     if not f.startswith("/") and f != "../"]
            dirs[acc] = files
            for f in files:
                if f.startswith("GCST") and f.endswith((".tsv", ".tsv.gz")):
                    try:
                        r = urllib.request.urlopen(
                            urllib.request.Request(base + f, method="HEAD"), timeout=10)
                        sizes.setdefault(acc, {})[f] = int(r.headers.get("Content-Length", 0))
                    except Exception:
                        sizes.setdefault(acc, {})[f] = None
            print(f"  {acc}: {len(files)} files", flush=True)
        except Exception as e:
            dirs[acc] = []
            print(f"  {acc}: DIRECTORY FAILED {e}", flush=True)
        time.sleep(0.2)
    (DATA / "sumstats_directory.json").write_text(json.dumps(dirs, indent=1), encoding="utf-8")
    (DATA / "sumstats_sizes.json").write_text(json.dumps(sizes, indent=1), encoding="utf-8")

    # 4. Nature Genetics 2025 (PMID 41188533) public Supplementary Tables XLSX.
    # Contains Supplementary Table 2 (94 lead loci) and Table 3 (prior-report comparison).
    supp = ("https://static-content.springer.com/esm/"
            "art%3A10.1038%2Fs41588-025-02393-x/MediaObjects/"
            "41588_2025_2393_MOESM3_ESM.xlsx")
    target = DATA / "ng2025_supp_tables.xlsx"
    if not target.exists():
        req = urllib.request.Request(supp, headers={"User-Agent": "Mozilla/5.0"})
        payload = urllib.request.urlopen(req, timeout=90).read()
        target.write_bytes(payload)
        print(f"  downloaded NG2025 supplementary tables: {len(payload)} bytes", flush=True)
    else:
        print("  NG2025 supplementary tables already present", flush=True)

    # 5. Source Data Fig. 3 (public GEO GSE118127 single-cell analysis) and Fig. 4.
    for n, nm in [(7, "ng2025_fig3_sourcedata.xlsx"), (8, "ng2025_fig4_sourcedata.xlsx")]:
        target = DATA / nm
        if target.exists():
            continue
        u = ("https://static-content.springer.com/esm/"
             "art%3A10.1038%2Fs41588-025-02393-x/MediaObjects/"
             f"41588_2025_2393_MOESM{n}_ESM.xlsx")
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            payload = urllib.request.urlopen(req, timeout=90).read()
            target.write_bytes(payload)
            print(f"  downloaded {nm}: {len(payload)} bytes", flush=True)
        except Exception as e:
            print(f"  WARNING could not download {nm}: {e}", flush=True)


if __name__ == "__main__":
    main()
