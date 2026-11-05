#!/usr/bin/env python3
"""Apply the frozen AFR GeneCards full-containment analysis to EUR."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

BASE = Path(repro_path("height_validation_workspace"))
OUT = BASE / "height_ten_seed_region_stability_q0p1"
GENES = BASE / "reference/genecards_height_top1000_hg19_gene_positions.csv"


def read(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    inputs = [GENES]
    source = read(GENES, ",")
    assert {int(g["genecards_rank"]) for g in source} == set(range(1, 1001))
    genes = [g for g in source if g["mapping_status"].startswith("mapped_")
             and g["autosomal_1_22"] == "TRUE"]
    assert len(genes) == 936
    for g in genes:
        assert g["coordinate_system"] == "GRCh37/hg19"
        assert g["interval_convention"] == "1-based closed"
        assert g["genecards_type"] != "Functional Element"
        assert 1 <= int(g["start_1based"]) <= int(g["end_1based"])
        assert g["ensembl_gene_id"]

    results = {}
    pairs_by_ancestry = {}
    for ancestry, expected in [("AFR", (9, 12, 6)), ("EUR", (87, 105, 12))]:
        path = OUT / f"{ancestry}_all10_union_regions.tsv"
        inputs.append(path)
        regions = read(path)
        classes = ["Both", "AnnoGK_only", "GhostKnockoff_only"]
        assert Counter(r["all_10_class"] for r in regions) == dict(zip(classes, expected))
        assert len({r["region_id"] for r in regions}) == len(regions)
        pairs = []
        for r in regions:
            assert r["ancestry"] == ancestry
            assert int(r["seeds_with_region_object"]) == 10
            flags = [int(r[k + "_n_seeds"]) == 10 for k in ["annogk", "ghost"]]
            assert any(flags)
            assert r["all_10_class"] == ("Both" if all(flags) else
                                         "AnnoGK_only" if flags[0] else "GhostKnockoff_only")
            for method, selected in zip(["annogk", "ghost"], flags):
                assert (r[method + "_all_10"] == "TRUE") == selected
                if selected:
                    assert float(r[method + "_worst_detected_min_q"]) <= .1
            # Preserve original padded region identity, including a EUR chr11
            # start below 1. All gene starts are >=1, so clipping that start
            # to 1 would produce exactly the same containment memberships.
            assert int(r["start"]) <= int(r["end"]) and int(r["end"]) >= 1
            for g in genes:
                if (g["chromosome"] == r["chromosome"]
                        and int(r["start"]) <= int(g["start_1based"])
                        and int(g["end_1based"]) <= int(r["end"])):
                    pairs.append({
                        "region_id": r["region_id"], "all_10_class": r["all_10_class"],
                        "chromosome": r["chromosome"], "region_start": r["start"],
                        "region_end": r["end"], "annogk_all_10": flags[0],
                        "ghost_all_10": flags[1], "genecards_rank": int(g["genecards_rank"]),
                        "genecards_symbol": g["genecards_symbol"],
                        "ensembl_gene_id": g["ensembl_gene_id"],
                        "gene_start": g["start_1based"], "gene_end": g["end_1based"],
                        "genecards_type": g["genecards_type"],
                    })
        counts = []
        for n in range(100, 1001, 100):
            hits = [p for p in pairs if p["genecards_rank"] <= n]
            sets = {m: {p["ensembl_gene_id"] for p in hits if p[m + "_all_10"]}
                    for m in ["annogk", "ghost"]}
            # Verify rank-, symbol-, and Ensembl-based deduplication agree.
            for m in sets:
                for key in ["genecards_rank", "genecards_symbol"]:
                    assert len(sets[m]) == len({p[key] for p in hits if p[m + "_all_10"]})
            a, g = sets["annogk"], sets["ghost"]
            counts.append(dict(top_n=n, overlap_rule="fully_contained_gene",
                               shared_genes=len(a & g), annogk_only_genes=len(a - g),
                               ghost_only_genes=len(g - a), annogk_total_genes=len(a),
                               ghost_total_genes=len(g), annogk_minus_ghost=len(a)-len(g)))
        results[ancestry] = counts
        pairs_by_ancestry[ancestry] = pairs

    old_path = OUT / "AFR_all10_genecards_height_topN_strict_gene_counts.tsv"
    inputs.append(old_path)
    assert [{k: str(v) for k, v in row.items()} for row in results["AFR"]] == read(old_path)
    # Also reproduce all stored AFR strict region-gene memberships.
    old_pairs_path = OUT / "AFR_all10_union_regions_genecards_height_top1000_overlap_long.tsv"
    inputs.append(old_pairs_path)
    expected_pairs = {(p["region_id"], int(p["genecards_rank"])) for p in read(old_pairs_path)
                      if p["fully_contained_in_region"] == "TRUE"}
    assert expected_pairs == {(p["region_id"], p["genecards_rank"])
                              for p in pairs_by_ancestry["AFR"]}

    prefix = "EUR_all10_genecards_height_topN_strict"
    write(OUT / f"{prefix}_gene_counts.tsv", results["EUR"])
    write(OUT / f"{prefix}_region_gene_pairs.tsv", pairs_by_ancestry["EUR"])
    membership = {}
    for p in pairs_by_ancestry["EUR"]:
        for method in ["annogk", "ghost"]:
            if p[method + "_all_10"]:
                key = (method, p["ensembl_gene_id"])
                row = membership.setdefault(key, dict(method=method,
                    ensembl_gene_id=p["ensembl_gene_id"], genecards_symbol=p["genecards_symbol"],
                    genecards_rank=p["genecards_rank"], region_ids=set()))
                row["region_ids"].add(p["region_id"])
    rows = sorted(membership.values(), key=lambda r: (r["method"], r["genecards_rank"]))
    for row in rows:
        row["region_ids"] = ";".join(sorted(row["region_ids"]))
    write(OUT / f"{prefix}_method_genes.tsv", rows)
    metadata = dict(
        analysis="EUR all-ten-seed GeneCards Height Top-N full-containment coverage",
        seeds=list(range(100, 1001, 100)), M=1, q=.1,
        regions=dict(AnnoGK=192, GhostKnockoff=99, shared=87),
        gene_reference="Frozen 2026-08-30 GeneCards input and GENCODE v19 GRCh37/hg19",
        rule="Same chromosome; region_start <= gene_start and gene_end <= region_end",
        deduplication="Stable Ensembl gene ID within method; shared regions included in both",
        cutoffs="Original GeneCards ranks 1..N, without reranking autosomal mapped genes",
        validation="All AFR strict counts and region-gene memberships reproduced exactly; rank/symbol/ID counts agree",
        boundary_note="Original EUR chr11:-46212-3368517 retained; clipping its start to 1 leaves gene containment unchanged",
        input_sha256={str(p.relative_to(BASE)) if p.is_relative_to(BASE) else p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (OUT / f"{prefix}_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(results["EUR"], indent=2))


if __name__ == "__main__":
    main()
