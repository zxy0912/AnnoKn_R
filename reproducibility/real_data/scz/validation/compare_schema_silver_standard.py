"""Compare SCHEMA Table S5 panels with SCZ >=2/100-seed region discoveries."""
import csv
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

HERE = Path(repro_path("scz_validation_workspace"))
BOOK = HERE / "reference/41586_2022_4556_MOESM3_ESM.xlsx"
PANEL = HERE / "reference/Singh_2022_TableS5_FDR_lt0p05_32genes.tsv"
GTF = Path(repro_path("height_validation_workspace", "reference/gencode.v19.annotation.gtf.gz"))
REGIONS = HERE / "hundred_seed_region_stability_q0p1/region_seed_recurrence.tsv"
OUT = HERE / "schema_silver_standard_seed_ge2"
METHODS = ["GhostKnockoff", "AnnoGK", "Transfer", "AdaKn"]
CHROMS = {str(i) for i in range(1, 23) if i != 21}


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(name, rows, fields=None):
    with (OUT / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def coverage_bp(regions):
    """Length of interval union on analyzed chromosomes (1-based closed)."""
    total = 0
    for chrom in CHROMS:
        intervals = sorted((max(1, int(r["start"])), int(r["end"]))
                           for r in regions if r["chr"] == chrom)
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        total += sum(end - start + 1 for start, end in merged)
    return total


# Independently check the frozen extracted panel against the official workbook.
ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
with zipfile.ZipFile(BOOK) as z:
    strings = ["".join(e.itertext()) for e in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    workbook = ET.fromstring(z.read("xl/workbook.xml"))
    assert workbook.find("x:sheets", ns)[2].get("name") == "Table S5 - Gene Results"
    source = {}
    for row in ET.fromstring(z.read("xl/worksheets/sheet3.xml")).findall("x:sheetData/x:row", ns)[1:]:
        cells = {}
        for cell in row.findall("x:c", ns):
            v = cell.find("x:v", ns)
            if v is not None:
                cells[re.sub(r"\d", "", cell.get("r"))] = strings[int(v.text)] if cell.get("t") == "s" else v.text
        if "Q" in cells and float(cells["Q"]) < 0.05:
            source[cells["A"]] = (float(cells["P"]), float(cells["Q"]))
panel = read(PANEL)
assert len(panel) == len(source) == 32
assert {r["symbol"] for r in panel} == set(source)
for r in panel:
    assert (float(r["p_meta"]), float(r["q_meta"])) == source[r["symbol"]]
    r["strict10"] = float(r["p_meta"]) < 2.14e-6
assert sum(r["strict10"] for r in panel) == 10
models = {}
with gzip.open(GTF, "rt") as f:
    headers = []
    for line in f:
        if line.startswith("#"):
            headers.append(line.strip())
            continue
        z = line.rstrip().split("\t")
        if z[2] != "gene":
            continue
        a = dict(re.findall(r'(\w+) "([^"]*)"', z[8]))
        if a["gene_name"] in source:
            models.setdefault(a["gene_name"], []).append(dict(chr=z[0].removeprefix("chr"),
                gene_start=int(z[3]), gene_end=int(z[4]), gene_id=a["gene_id"].split(".")[0]))
assert "GRCh37" in " ".join(headers)
genes = []
for r in panel:
    hits = models.get(r["symbol"], [])
    assert len(hits) == 1, f"Review ambiguous/unmapped gene {r['symbol']}"
    genes.append(dict(**r, **hits[0], analyzed_chromosome=hits[0]["chr"] in CHROMS))
assert len({g["gene_id"] for g in genes}) == 32
regions = read(REGIONS)
assert all(int(r["n_seeds_observed"]) == 100 for r in regions)
pairs, statuses = [], []
for ancestry in ["EUR", "EAS"]:
    rr = [r for r in regions if r["ancestry"] == ancestry]
    for gene in genes:
        overlapping = [r for r in rr if r["chr"] == gene["chr"] and
                       int(r["start"]) <= gene["gene_end"] and gene["gene_start"] <= int(r["end"])]
        contained = [r for r in overlapping if int(r["start"]) <= gene["gene_start"] and
                     gene["gene_end"] <= int(r["end"])]
        status = "outside_analyzed_chromosomes" if not gene["analyzed_chromosome"] else (
            "fully_contained_in_candidate_region" if contained else
            "partial_overlap_only" if overlapping else "outside_candidate_regions")
        statuses.append(dict(ancestry=ancestry, **gene, candidate_region_status=status,
                            contained_region_ids=";".join(r["region_id"] for r in contained)))
        for r in overlapping:
            pairs.append(dict(ancestry=ancestry, region_id=r["region_id"],
                region_start=int(r["start"]), region_end=int(r["end"]), **gene,
                fully_contained=r in contained, **{m: int(r[m]) for m in METHODS}))
summary, selected_pairs = [], []
for ancestry in ["EUR", "EAS"]:
    for method in METHODS:
        selected = [r for r in regions if r["ancestry"] == ancestry and int(r[method]) >= 2]
        for panel_name in ["strict10", "expanded32"]:
            gg = [g for g in genes if panel_name == "expanded32" or g["strict10"]]
            hits = [p for p in pairs if p["ancestry"] == ancestry and p["fully_contained"]
                    and p[method] >= 2 and (panel_name == "expanded32" or p["strict10"])]
            ids = {p["gene_id"] for p in hits}
            summary.append(dict(ancestry=ancestry, method=method, reference_panel=panel_name,
                total_reference_genes=len(gg), reference_genes_on_analyzed_chromosomes=sum(g["analyzed_chromosome"] for g in gg),
                reference_genes_contained_in_candidate_regions=sum(s["ancestry"] == ancestry and
                    s["candidate_region_status"] == "fully_contained_in_candidate_region" and
                    (panel_name == "expanded32" or s["strict10"]) for s in statuses),
                selected_regions=len(selected), selected_region_union_bp=coverage_bp(selected),
                regions_with_reference_gene=len({p["region_id"] for p in hits}),
                covered_reference_genes=len(ids), genes=";".join(sorted({p["symbol"] for p in hits}))))
            selected_pairs.extend(dict(method=method, reference_panel=panel_name, **p) for p in hits)
OUT.mkdir(exist_ok=True)
write("reference_gene_coordinates.tsv", genes)
write("reference_gene_candidate_region_status.tsv", statuses)
write("all_candidate_region_gene_overlaps.tsv", pairs)
write("method_gene_summary.tsv", summary)
write("selected_region_gene_pairs.tsv", selected_pairs,
      ["method", "reference_panel"] + list(pairs[0]))
metadata = dict(reference_doi="10.1038/s41586-022-04556-w", reference_table="Table S5 - Gene Results",
    strict10_rule="P meta < 2.14e-6", expanded32_rule="Q meta < 0.05; includes strict10",
    coordinate_source="GENCODE v19 gene features; exact original gene symbols; GRCh37/hg19; 1-based closed",
    discovery_rule="Region detected in >=2 of 100 seeds; M=1; q<=0.1; stored AdaKn rejection indices",
    seed_values=list(range(100, 10001, 100)), chromosomes=sorted(CHROMS, key=int),
    containment_rule="Same chromosome, region_start <= gene_start and gene_end <= region_end",
    deduplication="Stable Ensembl gene ID within method/ancestry/panel; shared regions included for each method",
    build_evidence="SCZ genome build follows manuscript GRCh37 statement and direct-position matching in SCZ_new.R; original PLINK reference build remains unaudited",
    excluded_chromosome_genes=[g["symbol"] for g in genes if not g["analyzed_chromosome"]],
    interpretation="Positional reference-gene coverage, not causal gene identification, formal enrichment, FDR estimation, or necessarily independent replication",
    sha256={str(p):sha256(p) for p in [BOOK, PANEL, GTF, REGIONS, Path(__file__)]})
(OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
for row in summary:
    print(row)
print("Selected reference-gene region recurrence:")
for p in pairs:
    if p["fully_contained"] and any(p[m] >= 2 for m in METHODS):
        print({k:p[k] for k in ["ancestry", "symbol", "strict10", "region_id"] + METHODS})
print("Outside analyzed chromosomes:", metadata["excluded_chromosome_genes"])
print("Outputs:", OUT)
