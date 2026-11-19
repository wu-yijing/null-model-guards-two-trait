"""T2D case study: fetch GWAS Catalog v2 metadata + lead associations for
type 2 diabetes mellitus. Public data, no credentials. Mirrors p1a/fetch_p1a.py.

Run: python t2a_fetch.py
Outputs: t2a/data/catalog_t2d_studies.json, catalog_t2d_associations.json
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
API = "https://www.ebi.ac.uk/gwas/api/v2"
TRAIT = "type 2 diabetes mellitus"
UA = {"Accept": "application/json", "User-Agent": "T2D-guard-case/1.0"}


def get_json(url, timeout=45, retries=3):
    last = None
    for k in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.load(urllib.request.urlopen(req, timeout=timeout))
        except Exception as e:
            last = e
            time.sleep(1.5 * (k + 1))
    raise last


def main():
    # 1. all studies, paginated (size=100)
    studies = []
    page = 0
    while True:
        url = f"{API}/studies?efoTrait={urllib.parse.quote(TRAIT)}&size=100&page={page}"
        d = get_json(url)
        chunk = d["_embedded"]["studies"]
        studies += chunk
        tp = d["page"]["totalPages"]
        print(f"studies page {page+1}/{tp}: +{len(chunk)} (total {len(studies)})", flush=True)
        if page >= tp - 1:
            break
        page += 1
        time.sleep(0.25)
    (DATA / "catalog_t2d_studies.json").write_text(
        json.dumps(studies, indent=1, ensure_ascii=False), encoding="utf-8")

    # 2. lead associations for every study with >=1 catalogue association
    assoc = {}
    for i, s in enumerate(studies):
        acc = s["accessionId"]
        if not s.get("associationCount"):
            assoc[acc] = []
            continue
        try:
            d = get_json(f"{API}/studies/{acc}/associations?size=1000")
            assoc[acc] = d["_embedded"]["associations"]
        except Exception as e:
            assoc[acc] = []
            print(f"  {acc}: FAILED {e}", flush=True)
        if (i + 1) % 25 == 0:
            print(f"  associations fetched for {i+1}/{len(studies)} studies", flush=True)
        time.sleep(0.15)
    (DATA / "catalog_t2d_associations.json").write_text(
        json.dumps(assoc, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"total associations: {sum(len(v) for v in assoc.values())}", flush=True)


if __name__ == "__main__":
    main()
