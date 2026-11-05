"""Fixed seeds 100,...,1000: region union coverage of frozen gene panels."""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(repro_path("scz_validation_workspace"))
OUT = HERE / "first10_seed_union_gene_overlap"
METHODS = ["GhostKnockoff", "AnnoGK", "Transfer", "AdaKn"]
SEEDS = set(range(100, 1001, 100))
REGIONS = HERE / "first10_seed_region_stability_q0p1/region_seed_recurrence.tsv"
DETECTIONS = HERE / "hundred_seed_region_stability_q0p1/region_by_seed_detections.tsv"
GC = HERE / "genecards_relevance_gt10_seed_ge2/reference_gene_coordinates.tsv"
SCHEMA = HERE / "schema_silver_standard_seed_ge2/reference_gene_coordinates.tsv"


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(name, rows):
    with (OUT / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


regions = read(REGIONS)
assert len(regions) == 338
counts = defaultdict(set)
seen = defaultdict(set)
for r in read(DETECTIONS):
    seed = int(r["seed"])
    if seed not in SEEDS:
        continue
    seen[r["ancestry"], r["region_id"]].add(seed)
    for m in METHODS:
        if r[m] == "TRUE":
            counts[r["ancestry"], r["region_id"], m].add(seed)
for r in regions:
    assert seen[r["ancestry"], r["region_id"]] == SEEDS
    for m in METHODS:
        assert int(r[m]) == len(counts[r["ancestry"], r["region_id"], m])
gc = read(GC)
schema = read(SCHEMA)
for g in schema:
    g["start"], g["end"] = g["gene_start"], g["gene_end"]
panels = {
    "GeneCards_gt10_GENCODE_only": [g for g in gc if g["source"].startswith("GENCODE")],
    "GeneCards_gt10_with_NCBI_RNA": gc,
    "SCHEMA_strict10": [g for g in schema if g["strict10"] == "True"],
    "SCHEMA_expanded32": schema,
}
summary, pairs = [], []
for panel, genes in panels.items():
    for a in ["EUR", "EAS"]:
        for m in METHODS:
            selected = [r for r in regions if r["ancestry"] == a and int(r[m]) >= 1]
            hits = []
            for r in selected:
                for g in genes:
                    if r["chr"] == g["chr"] and int(r["start"]) <= int(g["start"]) and int(g["end"]) <= int(r["end"]):
                        hits.append(g)
                        pairs.append(dict(panel=panel, ancestry=a, method=m,
                            region_id=r["region_id"], gene_id=g["gene_id"], symbol=g["symbol"],
                            gene_start=g["start"], gene_end=g["end"],
                            **{method: int(r[method]) for method in METHODS}))
            ids = {g["gene_id"] for g in hits}
            summary.append(dict(panel=panel, ancestry=a, method=m,
                union_regions=len(selected), reference_genes=len(genes), covered_genes=len(ids),
                genes=";".join(sorted({g["symbol"] for g in hits}))))
OUT.mkdir(exist_ok=True)
write("method_gene_summary.tsv", summary)
write("selected_strict_region_gene_pairs.tsv", pairs)
metadata = dict(seeds=sorted(SEEDS), region_rule="Union within the ten fixed seeds; >=1/10",
    M=1, q_threshold=0.1, chromosomes=[i for i in range(1, 23) if i != 21],
    gene_rule="Full gene body contained within a single selected region on same chromosome; GRCh37/hg19",
    deduplication="Gene ID per method/ancestry/panel; shared regions count for all selecting methods",
    genecards_reference="38 entries with relevance >10 in supplied top100 CSV; 5 functional elements excluded; 33 genes (31 GENCODE v19 plus 2 NCBI RNA genes)",
    schema_reference="Singh 2022 Table S5: 10 exome-wide significant genes, 32 Q meta <0.05 genes including the 10",
    coordinate_caveat="Original PLINK reference-panel build remains unaudited; follows manuscript GRCh37 statement",
    validation="All 338 region/ancestry rows and four method counts match distinct-seed reconstruction from detection table",
    sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [REGIONS, DETECTIONS, GC, SCHEMA, Path(__file__)]})
(OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
for row in summary:
    print(row)
print("Distinct hit-region recurrence:")
printed = set()
for row in pairs:
    key = row["ancestry"], row["region_id"], row["symbol"]
    if key not in printed:
        printed.add(key)
        print({k:row[k] for k in ["ancestry", "region_id", "symbol"] + METHODS})
