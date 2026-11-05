#!/usr/bin/env python3
"""Compare GTEx v8 support for old-seed AnnoKn-only and Knockoff-only eGenes.

The two input tables must have been extracted under the same analysis setting.
For the current revision analysis that setting is old seed, promoter + enhancer
+ TFBS annotations, knockoff+ with offset 1, and q <= 0.10.
"""

from __future__ import annotations

import os
os.chdir(repro_path("gtex_validation_workspace"))

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


API_URL = "https://gtexportal.org/api/v2/association/egene"
DEFAULT_CANDIDATE_DIR = Path(
    "gtex_annokn_only_eGenes_q01_initial_seed_PET"
)
DEFAULT_OUTPUT_DIR = Path(
    "gtex_method_specific_eGene_gtex_v8_comparison_q01_old_seed_PET"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--annokn-only",
        type=Path,
        default=DEFAULT_CANDIDATE_DIR / "annokn_only_egenes.tsv",
    )
    parser.add_argument(
        "--knockoff-only",
        type=Path,
        default=DEFAULT_CANDIDATE_DIR / "knockoff_only_egenes.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def stable_gene_id(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def request_egenes(tissue_id: str) -> list[dict]:
    query = urllib.parse.urlencode(
        {
            "tissueSiteDetailId": tissue_id,
            "datasetId": "gtex_v8",
            "itemsPerPage": 100000,
            "page": 0,
        }
    )
    request = urllib.request.Request(
        f"{API_URL}?{query}",
        headers={"User-Agent": "AnnoKn-method-specific-eGene-audit/1.0"},
    )
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.load(response)
            records = payload.get("data", [])
            paging = payload.get("paging_info", {})
            expected = paging.get("totalNumberOfItems", len(records))
            if len(records) != expected:
                raise RuntimeError(
                    f"Incomplete API response for {tissue_id}: "
                    f"received {len(records)} of {expected} records"
                )
            return records
        except Exception as error:  # retry transient API/network failures
            last_error = error
            time.sleep(2**attempt)
    raise RuntimeError(f"GTEx API failed for {tissue_id}") from last_error


def validate_candidates(rows: list[dict[str, str]], expected_category: str) -> None:
    required = {
        "tissue",
        "tissue_key",
        "ensembl_gene_id",
        "selection_category",
    }
    missing = required.difference(rows[0] if rows else {})
    if missing:
        raise RuntimeError(f"Input is missing columns: {sorted(missing)}")
    wrong = Counter(row["selection_category"] for row in rows)
    if set(wrong) != {expected_category}:
        raise RuntimeError(
            f"Expected only {expected_category} rows; observed {dict(wrong)}"
        )
    keys = [
        (row["tissue_key"], stable_gene_id(row["ensembl_gene_id"]))
        for row in rows
    ]
    duplicates = [key for key, count in Counter(keys).items() if count > 1]
    if duplicates:
        raise RuntimeError(
            f"Duplicated tissue-gene pairs in {expected_category}: {duplicates[:5]}"
        )


def audit_rows(
    rows: list[dict[str, str]],
    method_specific_group: str,
    api_by_tissue: dict[str, dict[str, dict]],
) -> list[dict]:
    audited: list[dict] = []
    for original in rows:
        gene_id = stable_gene_id(original["ensembl_gene_id"])
        match = api_by_tissue[original["tissue_key"]].get(gene_id)
        row = dict(original)
        row.update(
            {
                "method_specific_group": method_specific_group,
                "gtex_v8_standard_egene": "TRUE" if match else "FALSE",
                "gtex_v8_gencode_id": match["gencodeId"] if match else "",
                "gtex_v8_gene_symbol": match["geneSymbol"] if match else "",
                "gtex_v8_q_value": match["qValue"] if match else "",
                "gtex_v8_empirical_p_value": (
                    match["empiricalPValue"] if match else ""
                ),
                "gtex_v8_lead_nominal_p_value": match["pValue"] if match else "",
                "gtex_v8_p_value_threshold": (
                    match["pValueThreshold"] if match else ""
                ),
                "gtex_v8_log2_allelic_fold_change": (
                    match["log2AllelicFoldChange"] if match else ""
                ),
                "evidence_scope": (
                    "same-cohort GTEx v8 FastQTL"
                    if match
                    else "not a GTEx v8 standard eGene"
                ),
            }
        )
        audited.append(row)
    return audited


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    annokn_only = read_tsv(args.annokn_only)
    knockoff_only = read_tsv(args.knockoff_only)
    validate_candidates(annokn_only, "AnnoKn_only")
    validate_candidates(knockoff_only, "Knockoff_only")

    annokn_keys = {
        (row["tissue_key"], stable_gene_id(row["ensembl_gene_id"]))
        for row in annokn_only
    }
    knockoff_keys = {
        (row["tissue_key"], stable_gene_id(row["ensembl_gene_id"]))
        for row in knockoff_only
    }
    if annokn_keys.intersection(knockoff_keys):
        raise RuntimeError("AnnoKn-only and Knockoff-only inputs are not disjoint")

    all_candidates = annokn_only + knockoff_only
    tissue_ids = sorted({row["tissue_key"] for row in all_candidates})
    api_by_tissue: dict[str, dict[str, dict]] = {}
    api_inventory: list[dict] = []

    for tissue_id in tissue_ids:
        records = request_egenes(tissue_id)
        stable_ids = [stable_gene_id(row["gencodeId"]) for row in records]
        duplicates = [
            gene_id for gene_id, count in Counter(stable_ids).items() if count > 1
        ]
        if duplicates:
            raise RuntimeError(
                f"GTEx API returned duplicated stable gene IDs for {tissue_id}: "
                + ", ".join(duplicates[:5])
            )
        api_by_tissue[tissue_id] = {
            stable_gene_id(row["gencodeId"]): row for row in records
        }
        api_inventory.append(
            {
                "tissue_key": tissue_id,
                "n_official_gtex_v8_egenes": len(records),
            }
        )
        print(f"Downloaded {len(records):,} GTEx v8 eGenes for {tissue_id}")

    audited_annokn = audit_rows(annokn_only, "AnnoKn-only", api_by_tissue)
    audited_knockoff = audit_rows(
        knockoff_only, "Vanilla-Knockoff-only", api_by_tissue
    )
    audited = audited_annokn + audited_knockoff

    counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    for row in audited:
        counts[(row["method_specific_group"], row["tissue"])][
            row["gtex_v8_standard_egene"]
        ] += 1

    by_tissue: list[dict] = []
    for (group, tissue), count in sorted(counts.items()):
        supported = count["TRUE"]
        unsupported = count["FALSE"]
        total = supported + unsupported
        by_tissue.append(
            {
                "method_specific_group": group,
                "tissue": tissue,
                "n_pairs": total,
                "gtex_v8_standard_egene": supported,
                "not_gtex_v8_standard_egene": unsupported,
                "support_fraction": supported / total if total else "",
            }
        )

    overall: list[dict] = []
    for group in ("AnnoKn-only", "Vanilla-Knockoff-only"):
        group_rows = [row for row in audited if row["method_specific_group"] == group]
        supported = sum(row["gtex_v8_standard_egene"] == "TRUE" for row in group_rows)
        overall.append(
            {
                "method_specific_group": group,
                "n_pairs": len(group_rows),
                "gtex_v8_standard_egene": supported,
                "not_gtex_v8_standard_egene": len(group_rows) - supported,
                "support_fraction": supported / len(group_rows),
            }
        )

    detail_fields = list(audited[0])
    write_tsv(
        args.output_dir / "annokn_only_egenes_with_gtex_v8_support.tsv",
        audited_annokn,
        detail_fields,
    )
    write_tsv(
        args.output_dir / "vanilla_knockoff_only_egenes_with_gtex_v8_support.tsv",
        audited_knockoff,
        detail_fields,
    )
    write_tsv(
        args.output_dir / "method_specific_egenes_with_gtex_v8_support.tsv",
        audited,
        detail_fields,
    )
    write_tsv(
        args.output_dir / "gtex_v8_overlap_by_method.tsv",
        overall,
        list(overall[0]),
    )
    write_tsv(
        args.output_dir / "gtex_v8_overlap_by_method_and_tissue.tsv",
        by_tissue,
        list(by_tissue[0]),
    )
    write_tsv(
        args.output_dir / "gtex_v8_api_inventory.tsv",
        api_inventory,
        list(api_inventory[0]),
    )

    for row in overall:
        print(
            f"{row['method_specific_group']}: "
            f"{row['gtex_v8_standard_egene']}/{row['n_pairs']} "
            f"({100 * row['support_fraction']:.1f}%) GTEx v8 eGenes"
        )


if __name__ == "__main__":
    main()
