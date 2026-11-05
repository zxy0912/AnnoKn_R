#!/usr/bin/env python3
"""Track old-seed method-specific eGenes into new-seed selections."""

from __future__ import annotations

import os
os.chdir(repro_path("gtex_validation_workspace"))

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


OLD_DIR = Path("gtex_annokn_only_eGenes_q01_initial_seed_PET")
NEW_DIR = Path("gtex_method_specific_eGene_candidates_q01_new_seed_PET")
OUTPUT_DIR = Path("gtex_method_specific_eGene_seed_retention_q01_old_to_new_PET")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--old-annokn-only",
        type=Path,
        default=OLD_DIR / "annokn_only_egenes.tsv",
    )
    parser.add_argument(
        "--old-vanilla-only",
        type=Path,
        default=OLD_DIR / "knockoff_only_egenes.tsv",
    )
    parser.add_argument(
        "--new-all-selected",
        type=Path,
        default=NEW_DIR / "all_selected_egenes.tsv.gz",
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def stable_id(value: str) -> str:
    return value.split(".", 1)[0]


def key(row: dict[str, str]) -> tuple[str, str]:
    return row["source_file"], stable_id(row["ensembl_gene_id"])


def validate_group(rows: list[dict[str, str]], category: str) -> None:
    observed = Counter(row["selection_category"] for row in rows)
    if set(observed) != {category}:
        raise RuntimeError(f"Expected only {category}; observed {dict(observed)}")
    keys = [key(row) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError(f"Duplicated old-seed keys in {category}")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    old_annokn = read_tsv(args.old_annokn_only)
    old_vanilla = read_tsv(args.old_vanilla_only)
    new_selected = read_tsv(args.new_all_selected)
    validate_group(old_annokn, "AnnoKn_only")
    validate_group(old_vanilla, "Knockoff_only")

    new_index = {key(row): row for row in new_selected}
    if len(new_index) != len(new_selected):
        raise RuntimeError("Duplicated tissue-gene keys in new-seed selections")

    audited: list[dict] = []
    groups = (
        ("AnnoKn-only", "AnnoKn_only", old_annokn),
        ("Vanilla-Knockoff-only", "Knockoff_only", old_vanilla),
    )
    for group, same_only_category, candidates in groups:
        for candidate in candidates:
            new_row = new_index.get(key(candidate))
            new_category = new_row["selection_category"] if new_row else "Neither"
            new_annokn = new_category in {"AnnoKn_only", "Both"}
            new_vanilla = new_category in {"Knockoff_only", "Both"}
            corresponding = new_annokn if group == "AnnoKn-only" else new_vanilla
            row = dict(candidate)
            row.update(
                {
                    "old_seed_method_specific_group": group,
                    "new_seed_category": new_category,
                    "new_seed_annokn_egene": str(new_annokn).upper(),
                    "new_seed_vanilla_egene": str(new_vanilla).upper(),
                    "new_seed_any_method_egene": str(
                        new_category != "Neither"
                    ).upper(),
                    "retained_by_corresponding_method": str(corresponding).upper(),
                    "remains_same_only_category": str(
                        new_category == same_only_category
                    ).upper(),
                }
            )
            audited.append(row)

    category_levels = ("AnnoKn_only", "Both", "Knockoff_only", "Neither")
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in audited:
        grouped[(row["old_seed_method_specific_group"], row["tissue"])].append(row)
        grouped[(row["old_seed_method_specific_group"], "Overall")].append(row)

    summaries: list[dict] = []
    for (group, tissue), rows in sorted(grouped.items()):
        counts = Counter(row["new_seed_category"] for row in rows)
        corresponding = sum(
            row["retained_by_corresponding_method"] == "TRUE" for row in rows
        )
        any_method = sum(row["new_seed_any_method_egene"] == "TRUE" for row in rows)
        same_only = sum(row["remains_same_only_category"] == "TRUE" for row in rows)
        summary = {
            "old_seed_method_specific_group": group,
            "tissue": tissue,
            "old_seed_total": len(rows),
            "retained_by_corresponding_method": corresponding,
            "corresponding_method_retention_fraction": corresponding / len(rows),
            "retained_by_any_method": any_method,
            "any_method_retention_fraction": any_method / len(rows),
            "remains_same_only_category": same_only,
            "same_only_retention_fraction": same_only / len(rows),
        }
        summary.update(
            {f"new_{category}": counts[category] for category in category_levels}
        )
        summaries.append(summary)

    detail_fields = list(audited[0])
    summary_fields = list(summaries[0])
    write_tsv(
        args.output_dir / "old_method_specific_new_seed_status.tsv",
        audited,
        detail_fields,
    )
    write_tsv(
        args.output_dir / "old_method_specific_new_seed_summary.tsv",
        summaries,
        summary_fields,
    )

    for row in summaries:
        if row["tissue"] == "Overall":
            print(
                f"{row['old_seed_method_specific_group']}: "
                f"corresponding method "
                f"{row['retained_by_corresponding_method']}/"
                f"{row['old_seed_total']}; any method "
                f"{row['retained_by_any_method']}/{row['old_seed_total']}; "
                f"same-only {row['remains_same_only_category']}/"
                f"{row['old_seed_total']}"
            )


if __name__ == "__main__":
    main()
