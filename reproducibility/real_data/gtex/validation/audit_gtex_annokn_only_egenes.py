#!/usr/bin/env python3
"""Cross-reference AnnoKn-only eGenes with the official GTEx v8 eGene API."""

from __future__ import annotations

import os
os.chdir(repro_path("gtex_validation_workspace"))

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_INPUT = Path(
    "gtex_annokn_only_eGenes_q01_initial_seed_PET/annokn_only_egenes.tsv"
)
DEFAULT_OUTPUT_DIR = Path("gtex_annokn_only_eGenes_q01_initial_seed_PET")
API_URL = "https://gtexportal.org/api/v2/association/egene"


def stable_gene_id(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


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
        headers={"User-Agent": "AnnoKn-GTEx-eGene-audit/1.0"},
    )
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.load(response)
            data = payload.get("data", [])
            paging = payload.get("paging_info", {})
            expected = paging.get("totalNumberOfItems", len(data))
            if len(data) != expected:
                raise RuntimeError(
                    f"Incomplete API response for {tissue_id}: "
                    f"received {len(data)} of {expected} records"
                )
            return data
        except Exception as error:  # retry transient network/API failures
            last_error = error
            time.sleep(2**attempt)
    raise RuntimeError(f"GTEx API failed for {tissue_id}") from last_error


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


def main() -> None:
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INPUT
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    with input_path.open(newline="", encoding="utf-8") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))

    tissue_ids = sorted({row["tissue_key"] for row in candidates})
    api_by_tissue: dict[str, dict[str, dict]] = {}
    api_inventory: list[dict] = []

    for tissue_id in tissue_ids:
        records = request_egenes(tissue_id)
        indexed: dict[str, dict] = {}
        duplicates = Counter(stable_gene_id(row["gencodeId"]) for row in records)
        duplicate_ids = [gene_id for gene_id, count in duplicates.items() if count > 1]
        if duplicate_ids:
            raise RuntimeError(
                f"GTEx API returned duplicated stable gene IDs for {tissue_id}: "
                + ", ".join(duplicate_ids[:5])
            )
        for record in records:
            indexed[stable_gene_id(record["gencodeId"])] = record
        api_by_tissue[tissue_id] = indexed
        api_inventory.append(
            {
                "tissue_key": tissue_id,
                "n_official_gtex_v8_egenes": len(records),
            }
        )
        print(f"Downloaded {len(records):,} GTEx v8 eGenes for {tissue_id}")

    audited: list[dict] = []
    for candidate in candidates:
        gene_id = stable_gene_id(candidate["ensembl_gene_id"])
        match = api_by_tissue[candidate["tissue_key"]].get(gene_id)
        candidate = dict(candidate)
        candidate.update(
            {
                "gtex_v8_standard_egene": "TRUE" if match else "FALSE",
                "gtex_v8_gencode_id": match["gencodeId"] if match else "",
                "gtex_v8_gene_symbol": match["geneSymbol"] if match else "",
                "gtex_v8_q_value": match["qValue"] if match else "",
                "gtex_v8_empirical_p_value": (
                    match["empiricalPValue"] if match else ""
                ),
                "gtex_v8_lead_nominal_p_value": (
                    match["pValue"] if match else ""
                ),
                "gtex_v8_p_value_threshold": (
                    match["pValueThreshold"] if match else ""
                ),
                "gtex_v8_log2_allelic_fold_change": (
                    match["log2AllelicFoldChange"] if match else ""
                ),
                "evidence_scope": (
                    "same-cohort GTEx v8 FastQTL" if match else "not a GTEx v8 standard eGene"
                ),
            }
        )
        audited.append(candidate)

    summary_counts: dict[str, Counter] = defaultdict(Counter)
    for row in audited:
        summary_counts[row["tissue"]][row["gtex_v8_standard_egene"]] += 1

    summaries: list[dict] = []
    for tissue in sorted(summary_counts):
        supported = summary_counts[tissue]["TRUE"]
        unsupported = summary_counts[tissue]["FALSE"]
        total = supported + unsupported
        summaries.append(
            {
                "tissue": tissue,
                "annokn_only_pairs": total,
                "gtex_v8_standard_egene": supported,
                "not_gtex_v8_standard_egene": unsupported,
                "support_fraction": supported / total if total else "",
            }
        )

    output_fields = list(audited[0])
    write_tsv(
        output_dir / "annokn_only_egenes_with_gtex_v8_support.tsv",
        audited,
        output_fields,
    )
    write_tsv(
        output_dir / "gtex_v8_support_summary_by_tissue.tsv",
        summaries,
        list(summaries[0]),
    )
    write_tsv(
        output_dir / "gtex_v8_api_inventory.tsv",
        api_inventory,
        list(api_inventory[0]),
    )

    supported_total = sum(row["gtex_v8_standard_egene"] == "TRUE" for row in audited)
    print(
        f"Matched {supported_total:,}/{len(audited):,} AnnoKn-only pairs "
        "to official GTEx v8 eGenes."
    )


if __name__ == "__main__":
    main()
