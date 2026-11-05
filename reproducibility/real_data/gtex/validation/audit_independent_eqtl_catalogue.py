#!/usr/bin/env python3
"""Audit AnnoKn-only eGenes in independent eQTL Catalogue datasets."""

from __future__ import annotations

import os
os.chdir(repro_path("gtex_validation_workspace"))

import csv
import gzip
import io
import math
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path


DEFAULT_INPUT = Path(
    "gtex_annokn_only_eGenes_q01_initial_seed_PET/"
    "annokn_only_egenes_with_gtex_v8_support.tsv"
)
DEFAULT_OUTPUT_DIR = Path("gtex_annokn_only_eGenes_q01_initial_seed_PET")


DATASETS = [
    {
        "study": "FUSION",
        "dataset_id": "QTD000090",
        "study_id": "QTS000011",
        "source_tissue": "adipose",
        "scope": "same organ; adipose depot not resolved",
        "target_tissues": ["Adipose Subcutaneous"],
    },
    {
        "study": "TwinsUK",
        "dataset_id": "QTD000534",
        "study_id": "QTS000029",
        "source_tissue": "adipose",
        "scope": "same organ; adipose depot not resolved",
        "target_tissues": ["Adipose Subcutaneous"],
    },
    {
        "study": "FUSION",
        "dataset_id": "QTD000085",
        "study_id": "QTS000011",
        "source_tissue": "skeletal muscle",
        "scope": "same tissue",
        "target_tissues": ["Muscle Skeletal"],
    },
    {
        "study": "TwinsUK",
        "dataset_id": "QTD000544",
        "study_id": "QTS000029",
        "source_tissue": "skin",
        "scope": "same organ; sun-exposure site not resolved",
        "target_tissues": [
            "Skin Not Sun Exposed Suprapubic",
            "Skin Sun Exposed Lower leg",
        ],
    },
    {
        "study": "Lepik_2017",
        "dataset_id": "QTD000373",
        "study_id": "QTS000019",
        "source_tissue": "whole blood",
        "scope": "same tissue",
        "target_tissues": ["Whole Blood"],
    },
    {
        "study": "TwinsUK",
        "dataset_id": "QTD000549",
        "study_id": "QTS000029",
        "source_tissue": "whole blood",
        "scope": "same tissue",
        "target_tissues": ["Whole Blood"],
    },
    {
        "study": "van_de_Bunt_2015",
        "dataset_id": "QTD000554",
        "study_id": "QTS000030",
        "source_tissue": "pancreatic islet",
        "scope": "related pancreatic compartment; not whole pancreas",
        "target_tissues": ["Pancreas"],
    },
    {
        "study": "PISA",
        "dataset_id": "QTD000574",
        "study_id": "QTS000033",
        "source_tissue": "pancreatic islet",
        "scope": "related pancreatic compartment; not whole pancreas",
        "target_tissues": ["Pancreas"],
    },
]


def stable_gene_id(gene_id: str) -> str:
    return gene_id.split(".", 1)[0]


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=p_values.__getitem__)
    adjusted = [math.nan] * n
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = n - reverse_rank + 1
        running = min(running, p_values[index] * n / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def dataset_url(dataset: dict) -> str:
    return (
        "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
        f"{dataset['study_id']}/{dataset['dataset_id']}/"
        f"{dataset['dataset_id']}.permuted.tsv.gz"
    )


def load_permuted_results(dataset: dict) -> dict[str, dict]:
    request = urllib.request.Request(
        dataset_url(dataset),
        headers={"User-Agent": "AnnoKn-independent-eQTL-audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        compressed = response.read()
    text = gzip.decompress(compressed).decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))

    valid_indices: list[int] = []
    valid_p_values: list[float] = []
    for index, row in enumerate(rows):
        value = row.get("p_beta", "")
        try:
            p_value = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(p_value) and 0 <= p_value <= 1:
            valid_indices.append(index)
            valid_p_values.append(p_value)

    adjusted = bh_adjust(valid_p_values)
    for index, q_value in zip(valid_indices, adjusted):
        rows[index]["bh_q_value"] = q_value

    indexed: dict[str, dict] = {}
    for row in rows:
        gene_id = stable_gene_id(row["molecular_trait_id"])
        existing = indexed.get(gene_id)
        row_q = float(row.get("bh_q_value", "inf"))
        existing_q = float(existing.get("bh_q_value", "inf")) if existing else math.inf
        if existing is None or row_q < existing_q:
            indexed[gene_id] = row
    print(
        f"Loaded {len(rows):,} gene-level tests from "
        f"{dataset['study']} {dataset['dataset_id']} ({dataset['source_tissue']})"
    )
    return indexed


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

    results_by_dataset = {
        dataset["dataset_id"]: load_permuted_results(dataset)
        for dataset in DATASETS
    }

    audit_rows: list[dict] = []
    for candidate in candidates:
        gene_id = stable_gene_id(candidate["ensembl_gene_id"])
        applicable = [
            dataset
            for dataset in DATASETS
            if candidate["tissue"] in dataset["target_tissues"]
        ]
        supported: list[tuple[dict, dict]] = []
        tested: list[tuple[dict, dict]] = []
        for dataset in applicable:
            match = results_by_dataset[dataset["dataset_id"]].get(gene_id)
            if match is None:
                continue
            tested.append((dataset, match))
            q_value = float(match.get("bh_q_value", "inf"))
            if q_value <= 0.05:
                supported.append((dataset, match))

        best = min(
            supported,
            key=lambda item: float(item[1]["bh_q_value"]),
            default=None,
        )
        row = dict(candidate)
        row.update(
            {
                "independent_datasets_available": len(applicable),
                "independent_datasets_gene_tested": len(tested),
                "independent_eqtl_support": "TRUE" if supported else "FALSE",
                "independent_supporting_datasets": ";".join(
                    item[0]["dataset_id"] for item in supported
                ),
                "independent_supporting_studies": ";".join(
                    item[0]["study"] for item in supported
                ),
                "independent_best_bh_q_value": (
                    best[1]["bh_q_value"] if best else ""
                ),
                "independent_best_p_beta": best[1]["p_beta"] if best else "",
                "independent_best_nominal_p": best[1]["pvalue"] if best else "",
                "independent_best_lead_variant": best[1]["variant"] if best else "",
                "independent_best_lead_beta": best[1]["beta"] if best else "",
                "independent_best_scope": best[0]["scope"] if best else "",
            }
        )
        audit_rows.append(row)

    summary_accumulator: dict[str, dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    for row in audit_rows:
        summary = summary_accumulator[row["tissue"]]
        summary["annokn_only_pairs"] += 1
        if int(row["independent_datasets_available"]) > 0:
            summary["pairs_with_relevant_independent_dataset"] += 1
        if int(row["independent_datasets_gene_tested"]) > 0:
            summary["genes_tested_in_at_least_one_independent_dataset"] += 1
        if row["independent_eqtl_support"] == "TRUE":
            summary["independent_eqtl_supported"] += 1

    summary_rows: list[dict] = []
    for tissue in sorted(summary_accumulator):
        summary = summary_accumulator[tissue]
        tested = summary["genes_tested_in_at_least_one_independent_dataset"]
        supported = summary["independent_eqtl_supported"]
        summary_rows.append(
            {
                "tissue": tissue,
                "annokn_only_pairs": summary["annokn_only_pairs"],
                "pairs_with_relevant_independent_dataset": summary[
                    "pairs_with_relevant_independent_dataset"
                ],
                "genes_tested_in_at_least_one_independent_dataset": tested,
                "independent_eqtl_supported": supported,
                "support_fraction_among_tested": supported / tested if tested else "",
            }
        )

    write_tsv(
        output_dir / "annokn_only_egenes_with_independent_eqtl_support.tsv",
        audit_rows,
        list(audit_rows[0]),
    )
    write_tsv(
        output_dir / "independent_eqtl_support_summary_by_tissue.tsv",
        summary_rows,
        list(summary_rows[0]),
    )
    write_tsv(
        output_dir / "independent_eqtl_dataset_inventory.tsv",
        [
            {
                **dataset,
                "target_tissues": ";".join(dataset["target_tissues"]),
                "url": dataset_url(dataset),
            }
            for dataset in DATASETS
        ],
        [
            "study",
            "study_id",
            "dataset_id",
            "source_tissue",
            "scope",
            "target_tissues",
            "url",
        ],
    )

    total_supported = sum(
        row["independent_eqtl_support"] == "TRUE" for row in audit_rows
    )
    print(f"Independent eQTL support found for {total_supported} candidate pairs.")


if __name__ == "__main__":
    main()
