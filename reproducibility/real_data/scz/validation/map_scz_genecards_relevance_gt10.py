"""Map the user-supplied SCZ reference panel to >=2/100-seed regions.

Primary coordinate source matches Height: GENCODE v19 GRCh37 gene features.
Report a separate analysis supplementing two missing RNA genes with NCBI
GRCh37 coordinates. Never substitute MIR137HG for MIR137.
"""
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

HERE = Path(repro_path("scz_validation_workspace"))
REF = Path(repro_path("height_validation_workspace", "reference"))
INPUT = HERE / "genecards_SCZ_top100.csv"
REGIONS = HERE / "hundred_seed_region_stability_q0p1/region_seed_recurrence.tsv"
GTF = REF / "gencode.v19.annotation.gtf.gz"
OUT = HERE / "genecards_relevance_gt10_seed_ge2"
METHODS = ["GhostKnockoff", "AnnoGK", "Transfer", "AdaKn"]
# Frozen official NCBI Gene genomic-context records, 1-based closed coordinates.
SUPPLEMENT = {
    "MIR137": dict(chr="1", start=98511626, end=98511727, gene_id="NCBI:406928",
                   source="NCBI GRCh37.p13 RS_2024_09", source_url=
                   "https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&list_uids=406928"),
    "TRA-TGC7-1": dict(chr="6", start=28770577, end=28770647, gene_id="NCBI:7154",
                       source="NCBI GRCh37.p13 RS_2024_09", source_url=
                       "https://www.ncbi.nlm.nih.gov/gene/7154"),
}


def read(path, delimiter=","):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write(name, rows):
    assert rows
    with (OUT / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


assert digest(GTF, "md5") == "bd83e28270e595d3bde6bfcb21c9748f"
raw = read(INPUT)
assert len(raw) == 100 and [int(r["Rank"]) for r in raw] == list(range(1, 101))
selected = [r for r in raw if float(r["Relevance Score"]) > 10]
assert len(selected) == 38 and len({r["Symbol"] for r in selected}) == 38
models = {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        z = line.rstrip().split("\t")
        if z[2] != "gene":
            continue
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', z[8]))
        models.setdefault(attrs["gene_name"], []).append(dict(
            chr=z[0].removeprefix("chr"), start=int(z[3]), end=int(z[4]),
            gene_id=attrs["gene_id"].split(".")[0], source="GENCODE v19 GRCh37",
            source_url="https://www.gencodegenes.org/human/release_19.html"))

panel, coordinates = [], []
for r in selected:
    symbol = r["Symbol"]
    p = dict(rank=int(r["Rank"]), symbol=symbol, name=r["Name"],
             type=r["Type"], relevance_score=float(r["Relevance Score"]))
    if r["Type"] == "Functional Element":
        panel.append(dict(**p, status="excluded_functional_element"))
        continue
    assert r["Type"] in {"Protein Coding", "RNA Gene"}
    hits = models.get(symbol, [])
    if hits:
        assert len(hits) == 1, "Review multi-locus mapping explicitly"
        status = "GENCODE_v19_exact_gene_name"
    else:
        assert symbol in SUPPLEMENT, f"Unmapped gene: {symbol}"
        hits = [SUPPLEMENT[symbol]]
        status = "supplemental_NCBI_GRCh37"
    panel.append(dict(**p, status=status))
    for hit in hits:
        coordinates.append(dict(**p, **hit, status=status,
                                chromosome_analyzed=hit["chr"] in {str(i) for i in range(1, 23) if i != 21}))

assert len(coordinates) == 33
regions = read(REGIONS, "\t")
assert all(int(r["n_seeds_observed"]) == 100 for r in regions)
pairs = []
for r in regions:
    for g in coordinates:
        if r["chr"] != g["chr"]:
            continue
        if int(r["start"]) <= g["end"] and g["start"] <= int(r["end"]):
            strict = int(r["start"]) <= g["start"] and g["end"] <= int(r["end"])
            pairs.append(dict(ancestry=r["ancestry"], region_id=r["region_id"],
                region_start=int(r["start"]), region_end=int(r["end"]), **g,
                fully_contained=strict, **{m: int(r[m]) for m in METHODS}))

summary, method_genes = [], []
for scope in ["GENCODE_v19_only", "GENCODE_v19_plus_NCBI_RNA"]:
    for ancestry in ["EUR", "EAS"]:
        for method in METHODS:
            hits = [p for p in pairs if p["ancestry"] == ancestry and p["fully_contained"]
                    and p[method] >= 2 and (scope.endswith("RNA") or p["source"].startswith("GENCODE"))]
            ids = {p["gene_id"] for p in hits}
            symbols = sorted({p["symbol"] for p in hits})
            assert len(ids) == len(symbols)
            summary.append(dict(coordinate_scope=scope, ancestry=ancestry, method=method,
                selected_regions=sum(r["ancestry"] == ancestry and int(r[method]) >= 2 for r in regions),
                regions_with_reference_gene=len({p["region_id"] for p in hits}),
                distinct_genes=len(ids), genes=";".join(symbols)))
            for p in hits:
                method_genes.append(dict(coordinate_scope=scope, method=method, **p))

OUT.mkdir(exist_ok=True)
write("reference_entries_gt10.tsv", panel)
write("reference_gene_coordinates.tsv", coordinates)
write("all_candidate_region_gene_overlaps.tsv", pairs)
write("method_gene_summary.tsv", summary)
write("selected_strict_region_gene_pairs.tsv", method_genes)
metadata = dict(input_rows=100, score_rule="strictly >10", qualifying_entries=38,
    excluded_functional_elements=[p["symbol"] for p in panel if p["status"].startswith("excluded")],
    gene_entries=33, gencode_mapped=31, supplemental_NCBI_RNA=list(SUPPLEMENT),
    genes_on_excluded_chromosomes=[g["symbol"] for g in coordinates if not g["chromosome_analyzed"]],
    region_rule="At least 2 distinct seeds out of 100; M=1; q<=0.1; chr1-22 excluding chr21",
    AdaKn_rule="Use stored rejection indices as in original analysis",
    mapping_rule="Same chromosome and full gene body contained in a single region; 1-based closed coordinates",
    deduplication="Stable gene identifier within each method and ancestry; shared regions count for each method",
    build_evidence="Manuscript JASA_AnnoKn_v1_supple.tex states SCZ summary statistics use GRCh37/hg19; SCZ_new.R directly joins positions. Original PLINK reference-panel build not independently resolved.",
    supplementary_coordinates=SUPPLEMENT,
    limitations="GeneCards relevance is a ranking, not causal evidence or statistical significance. Region containment does not assign selected SNPs to causal genes. The 38 entries are a thresholded subset of the user-supplied top100 export.",
    sha256={str(p): digest(p) for p in [INPUT, REGIONS, GTF, Path(__file__)]})
(OUT / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
for row in summary:
    print(row)
print("Outputs:", OUT)
