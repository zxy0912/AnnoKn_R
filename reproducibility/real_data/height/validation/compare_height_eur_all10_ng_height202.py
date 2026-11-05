#!/usr/bin/env python3
"""Count the user-approved 202 reference genes fully inside EUR regions.

Gene coordinates come from GENCODE v19 gene records matched by stable
Ensembl ID. Table 16 DEPICT locus bounds are not used as gene coordinates.
"""

import csv
import gzip
import hashlib
import io
import json
import re
from pathlib import Path


BASE = Path(repro_path("height_validation_workspace"))
REFERENCE = BASE / "reference"
STABILITY = BASE / "height_ten_seed_region_stability_q0p1"
PANEL_FILE = REFERENCE / "ng_height_gene_evidence_audit_2026-09-07" / "ng_height_202_candidate_from_table9_stars_NOT_CONFIRMED.csv"
REGION_FILE = STABILITY / "EUR_all10_union_regions.tsv"
GTF_FILE = REFERENCE / "gencode.v19.annotation.gtf.gz"
PREFIX = "EUR_all10_ng_height202_strict"
CLASSES = ("Both", "AnnoGK_only", "GhostKnockoff_only")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_table(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_text(path, content):
    if path.exists():
        require(path.read_text() == content, f"Existing output differs: {path}")
    else:
        path.write_text(content, encoding="utf-8")


def write_table(path, rows, fields):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    write_text(path, buffer.getvalue())


def main():
    panel = read_table(PANEL_FILE, ",")
    require(len(panel) == 202, "Expected the authorized 202-gene panel")
    targets = {r["ensembl_gene_id"]: r for r in panel}
    require(len(targets) == 202, "Reference panel has duplicate Ensembl IDs")
    require(len({r["gene_symbol"] for r in panel}) == 202, "Reference symbols are not unique")
    regions = read_table(REGION_FILE)
    require(len({r["region_id"] for r in regions}) == len(regions), "Duplicate EUR region")
    for r in regions:
        for key in ("chromosome", "start", "end", "ghost_n_seeds", "annogk_n_seeds"):
            r[key] = int(r[key])
        require(r["ancestry"] == "EUR" and 1 <= r["chromosome"] <= 22, "Unexpected ancestry/chromosome")
        require(r["start"] <= r["end"], "Invalid region bounds")
        require(r["region_id"] == f"chr{r['chromosome']}:{r['start']}-{r['end']}", "Region ID/bounds mismatch")
        require(0 <= r["ghost_n_seeds"] <= 10 and 0 <= r["annogk_n_seeds"] <= 10, "Invalid seed count")
        ghost = r["ghost_n_seeds"] == 10
        anno = r["annogk_n_seeds"] == 10
        require(ghost or anno, "Region fails both all-ten criteria")
        expected_class = "Both" if ghost and anno else "AnnoGK_only" if anno else "GhostKnockoff_only"
        require(r["all_10_class"] == expected_class, "Region class/recurrence mismatch")
    require(sum(r["annogk_n_seeds"] == 10 for r in regions) == 192, "EUR AnnoGK handoff mismatch")
    require(sum(r["ghost_n_seeds"] == 10 for r in regions) == 99, "EUR GhostKnockoff handoff mismatch")

    # Check the analysis regions against the stored EUR risk-region files.
    for chromosome in {r["chromosome"] for r in regions}:
        risk_file = BASE / "annot_pvalus" / "risk_region" / f"EUR_chr_{chromosome}.txt"
        risk_intervals = {tuple(int(float(x)) for x in line.split())
                          for line in risk_file.read_text().splitlines() if line.strip()}
        require(all((r["start"], r["end"]) in risk_intervals
                    for r in regions if r["chromosome"] == chromosome),
                f"Stored EUR risk-region bounds differ on chromosome {chromosome}")

    gene_models = {}
    with gzip.open(GTF_FILE, "rt") as handle:
        for line in handle:
            if "\tgene\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            stable_id = attrs["gene_id"].split(".")[0]
            if stable_id not in targets:
                continue
            require(stable_id not in gene_models, f"Multiple GENCODE models for {stable_id}")
            target = targets[stable_id]
            require(fields[0] == "chr" + target["chromosome"], f"Gene chromosome mismatch: {stable_id}")
            require(attrs["gene_name"] == target["gene_symbol"], f"Gene symbol mismatch: {stable_id}")
            gene_models[stable_id] = {
                "ensembl_gene_id": stable_id,
                "gene_symbol": target["gene_symbol"],
                "gencode_gene_id_versioned": attrs["gene_id"],
                "gencode_gene_type": attrs["gene_type"],
                "gene_chromosome": int(fields[0][3:]),
                "gene_start": int(fields[3]),
                "gene_end": int(fields[4]),
                "gene_strand": fields[6],
                "mapping_rule": "exact_stable_Ensembl_ID",
                "coordinate_system": "GRCh37/hg19",
                "interval_convention": "1-based closed",
            }
    require(set(gene_models) == set(targets), "Unmapped reference genes: " + ",".join(sorted(set(targets) - set(gene_models))))
    mapped = [gene_models[r["ensembl_gene_id"]] for r in panel]
    write_table(STABILITY / (PREFIX + "_reference_gene_coordinates.tsv"), mapped, list(mapped[0]))

    gene_fields = list(mapped[0])
    region_fields = ["region_id", "chromosome", "start", "end", "all_10_class", "ghost_n_seeds", "annogk_n_seeds"]
    pairs, per_region, boundaries = [], [], []
    for r in regions:
        contained = []
        for g in mapped:
            same_chr = g["gene_chromosome"] == r["chromosome"]
            full = same_chr and r["start"] <= g["gene_start"] <= g["gene_end"] <= r["end"]
            overlap = same_chr and g["gene_start"] <= r["end"] and g["gene_end"] >= r["start"]
            combined = {**{k: r[k] for k in region_fields}, **g}
            if full:
                contained.append(g)
                pairs.append(combined)
            elif overlap:
                boundaries.append(combined)
        per_region.append({
            **{k: r[k] for k in region_fields},
            "n_reference_genes_fully_contained": len({g["ensembl_gene_id"] for g in contained}),
            "gene_symbols": ";".join(sorted(g["gene_symbol"] for g in contained)),
            "ensembl_gene_ids": ";".join(sorted(g["ensembl_gene_id"] for g in contained)),
        })
    require(len({(p["region_id"], p["ensembl_gene_id"]) for p in pairs}) == len(pairs), "Duplicate region-gene pair")
    summaries = []
    for label in CLASSES:
        these_regions = [r for r in per_region if r["all_10_class"] == label]
        these_pairs = [r for r in pairs if r["all_10_class"] == label]
        ids = {r["ensembl_gene_id"] for r in these_pairs}
        summaries.append({
            "region_class": label,
            "total_regions": len(these_regions),
            "regions_with_reference_gene": sum(r["n_reference_genes_fully_contained"] > 0 for r in these_regions),
            "regions_without_reference_gene": sum(r["n_reference_genes_fully_contained"] == 0 for r in these_regions),
            "region_gene_pairs": len(these_pairs),
            "unique_reference_genes": len(ids),
            "gene_symbols": ";".join(sorted(gene_models[x]["gene_symbol"] for x in ids)),
        })
        if label != "Both":
            genes = [{**gene_models[x], "region_ids": ";".join(sorted({p["region_id"] for p in these_pairs if p["ensembl_gene_id"] == x}))}
                     for x in sorted(ids, key=lambda ident: gene_models[ident]["gene_symbol"])]
            write_table(STABILITY / (PREFIX + f"_{label}_genes.tsv"), genes, gene_fields + ["region_ids"])
    write_table(STABILITY / (PREFIX + "_region_gene_pairs.tsv"), pairs, region_fields + gene_fields)
    write_table(STABILITY / (PREFIX + "_excluded_boundary_overlaps.tsv"), boundaries, region_fields + gene_fields)
    write_table(STABILITY / (PREFIX + "_per_region.tsv"), per_region, list(per_region[0]))
    write_table(STABILITY / (PREFIX + "_class_summary.tsv"), summaries, list(summaries[0]))

    membership = []
    for g in mapped:
        hits = [r for r in pairs if r["ensembl_gene_id"] == g["ensembl_gene_id"]]
        membership.append({**g, **{label + "_regions": ";".join(sorted(r["region_id"] for r in hits if r["all_10_class"] == label)) for label in CLASSES}})
    write_table(STABILITY / (PREFIX + "_all202_gene_membership.tsv"), membership, list(membership[0]))

    metadata = {
        "date": "2026-09-07",
        "reference_set": "User-approved 202 candidates reconstructed from Table 16 GRAIL and Table 9 stars",
        "reference_set_source": str(PANEL_FILE), "reference_set_sha256": sha256(PANEL_FILE),
        "reference_genes": len(panel), "mapped_genes": len(mapped),
        "gene_annotation": "GENCODE v19 GRCh37.p13, Ensembl 74",
        "gene_annotation_sha256": sha256(GTF_FILE),
        "regions_source": str(REGION_FILE), "regions_sha256": sha256(REGION_FILE),
        "region_class_definition": "Both methods separately require detection in all ten M=1 seeds 100,...,1000 at q<=0.1; only means the other method fails the all-ten criterion",
        "gene_membership_rule": "same chromosome AND region_start <= gene_start AND gene_end <= region_end",
        "gene_count_rule": "unique stable Ensembl IDs within each region class",
        "source_EUR_risk_region_bounds_match": True,
        "class_summary": summaries,
        "boundary_overlap_pairs_excluded": len(boundaries),
        "strict_genes_shared_across_classes": [g["gene_symbol"] for g in membership
                                                if sum(bool(g[label + "_regions"]) for label in CLASSES) > 1],
        "script_sha256": sha256(Path(__file__)),
    }
    write_text(STABILITY / (PREFIX + "_metadata.json"), json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")

    names = {"Both": "Both", "AnnoGK_only": "AnnoGK-only", "GhostKnockoff_only": "Ghost-only"}
    table = "\n".join(f"| {names[r['region_class']]} | {r['total_regions']} | {r['regions_with_reference_gene']} | {r['unique_reference_genes']} |"
                      for r in summaries)
    lists = "\n\n".join(f"**{names[r['region_class']]}:** {r['gene_symbols'].replace(';', ', ') or 'None'}." for r in summaries if r["region_class"] != "Both")
    report = f"""# EUR all-ten-seed regions versus the authorized 202 height genes

Date: 2026-09-07

The user chose the 202-gene candidate set reconstructed using Table 16
GRAIL support or Table 9 star markers. The provenance and original
202/204 discrepancy remain documented in
`Height_EUR_ng_height_supported_gene_audit_2026-09-07.md`.

All 220 EUR M=1 result files for seeds 100,200,...,1000 and chromosomes
1--22 were read. At every seed, a region is detected when at least one SNP
has a finite q-value <= 0.1. The complete recurrence curve reproduced the
August 20 handoff. At the all-ten criterion there are 192 AnnoGK regions
and 99 GhostKnockoff regions: 87 shared, 105 AnnoGK-only, and 12 Ghost-only.
Only refers to meeting this recurrence criterion, not absence of detection
by the other method under every individual seed.

All 202 reference genes matched a single GENCODE v19 gene record by stable
Ensembl ID, with identical gene symbols and chromosomes. Gene bodies and
regions use 1-based closed GRCh37/hg19 coordinates. A gene is counted only
if its entire body lies within the region. Table 16 locus bounds were not
used as gene coordinates. Each gene is counted once per region class.

| Region class | Regions | Regions containing a reference gene | Distinct reference genes fully contained |
|---|---:|---:|---:|
{table}

{lists}

The comparison requested by the user is between the two method-specific
rows; the shared regions are shown separately. This is positional gene
coverage within stable regions, not a causal-gene assignment.

## Reproducibility

```bash
Rscript reproducibility/run.R real_data/height/validation/extract_height_eur_ten_seed_stable_regions.R
Rscript reproducibility/run.R real_data/height/validation/compare_height_eur_all10_ng_height202.py
```

Outputs are under `height_ten_seed_region_stability_q0p1/`:

- `EUR_all10_union_regions.tsv`: complete stable-region union and classes.
- `EUR_ten_seed_input_manifest.tsv`: all 220 source files with MD5 checksums.
- `{PREFIX}_reference_gene_coordinates.tsv`: exact gene coordinates.
- `{PREFIX}_class_summary.tsv`: counts and gene lists by region class.
- `{PREFIX}_per_region.tsv`: one row per region, including zero-hit regions.
- `{PREFIX}_region_gene_pairs.tsv`: every strict region-gene match.
- `{PREFIX}_AnnoGK_only_genes.tsv` and `{PREFIX}_GhostKnockoff_only_genes.tsv`: distinct gene lists with their matched regions.
- `{PREFIX}_all202_gene_membership.tsv`: all reference genes, including those not contained in any stable region.
- `{PREFIX}_excluded_boundary_overlaps.tsv`: partial gene-body overlaps excluded by the strict rule.
- `{PREFIX}_metadata.json`: rules, input hashes, mapping checks and counts.
"""
    report_path = BASE / "Height_EUR_all10_ng_height202_strict_gene_comparison_2026-09-07.md"
    write_text(report_path, report)
    print(json.dumps({"mapped_reference_genes": len(mapped), "class_summary": summaries,
                      "strict_genes_shared_across_classes": metadata["strict_genes_shared_across_classes"],
                      "report": str(report_path)}, indent=2))


if __name__ == "__main__":
    main()
