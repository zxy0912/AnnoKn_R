#!/usr/bin/env python3
"""Audit Wood et al. (2014) Table 16 support flags against Table 13.

Uses the Python standard library to read the supplied single-sheet XLSX
files. Does not force the extracted set to match the article's count of 202.
"""

import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile


BASE = Path(repro_path("height_validation_workspace"))
REFERENCE = BASE / "reference"
OUTPUT = REFERENCE / "ng_height_gene_evidence_audit_2026-09-07"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
OFFICIAL_URLS = {
    13: "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fng.3097/MediaObjects/41588_2014_BFng3097_MOESM9_ESM.xlsx",
    16: "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fng.3097/MediaObjects/41588_2014_BFng3097_MOESM11_ESM.xlsx",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_xlsx(path):
    with ZipFile(path) as archive:
        strings = [
            "".join(t.text or "" for t in si.findall(".//m:t", NS))
            for si in ET.fromstring(archive.read("xl/sharedStrings.xml"))
            .findall("m:si", NS)
        ]
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        rows = {}
        for row in root.findall("m:sheetData/m:row", NS):
            require(row.get("hidden") != "1", "Unexpected hidden source row")
            cells = {}
            for cell in row.findall("m:c", NS):
                require(cell.find("m:f", NS) is None, "Unexpected source formula")
                value_node = cell.find("m:v", NS)
                value = value_node.text if value_node is not None else ""
                if cell.get("t") == "s" and value:
                    value = strings[int(value)]
                elif cell.get("t") == "inlineStr":
                    value = "".join(t.text or "" for t in cell.findall(".//m:t", NS))
                column = "".join(c for c in cell.get("r") if c.isalpha())
                cells[column] = (value or "").strip()
            rows[int(row.get("r"))] = cells
        return rows


def write_text(path, content):
    if path.exists():
        require(path.read_text() == content, f"Existing output differs: {path}")
    else:
        path.write_text(content, encoding="utf-8")


def write_table(name, rows, fields, delimiter="\t"):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter=delimiter,
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    write_text(OUTPUT / name, buffer.getvalue())


def main():
    tables = {n: read_xlsx(REFERENCE / f"supplementary{n}.xlsx") for n in (13, 16)}
    require(tables[13][2]["B"] == "Prioritized Gene", "Table 13 gene header changed")
    require(tables[16][3]["F"] == "Ensembl ID", "Table 16 ID header changed")
    require(tables[16][3]["K"] == "OMIM human stature gene", "OMIM header changed")
    require(tables[16][3]["L"] == "Nominally significant in GRAIL", "GRAIL header changed")
    table13 = [(n, r) for n, r in tables[13].items() if n >= 3 and r.get("B")]
    table16 = [(n, r) for n, r in tables[16].items() if n >= 4 and r.get("F")]
    require(len(table13) == 339, "Expected 339 Table 13 genes")
    require(len(table16) == 649, "Expected 649 Table 16 genes")
    require(len({r["B"] for _, r in table13}) == 339, "Duplicate Table 13 symbol")
    require(len({r["F"] for _, r in table16}) == 649, "Duplicate Table 16 Ensembl ID")
    require(all(0 <= float(r["C"]) < 0.01 for _, r in table13), "Unexpected Table 13 P-value")
    require(all(r["J"] == "Yes" for _, r in table16), "Table 16 includes a non-significant DEPICT row")
    require(all(r[c] in {"Yes", "No"} for _, r in table16 for c in ("K", "L")),
            "Unexpected OMIM/GRAIL flag")
    grail13 = {r["B"]: (n, r) for n, r in table13}

    audit = []
    for row_number, row in table16:
        matched_row, matched = grail13.get(row["G"], ("", {}))
        omim = row["K"] == "Yes"
        grail = row["L"] == "Yes"
        support = ("Both" if omim and grail else "OMIM_only" if omim else
                   "GRAIL_only" if grail else "Neither")
        audit.append({
            "supplementary16_excel_row": row_number,
            "ensembl_gene_id": row["F"],
            "gene_symbol": row["G"],
            "depict_locus": row["A"],
            "chromosome": row["B"],
            "depict_locus_start": row["C"],
            "depict_locus_end": row["D"],
            "depict_prioritization_pvalue": row["I"],
            "depict_fdr_lt_0p05": row["J"],
            "table16_omim": row["K"],
            "table16_grail_nominal": row["L"],
            "table16_support_class": support,
            "selected_by_table16_omim_or_grail": "Yes" if omim or grail else "No",
            "table13_exact_symbol_match": "Yes" if matched else "No",
            "table13_excel_row": matched_row,
            "table13_grail_locus": matched.get("A", ""),
            "table13_pvalue_697_snps": matched.get("C", ""),
            "table13_pvalue_previous_giant_snps": matched.get("D", ""),
            "selected_by_table13_or_table16_omim": "Yes" if matched or omim else "No",
        })

    selected = [r for r in audit if r["selected_by_table16_omim_or_grail"] == "Yes"]
    literal13 = [r for r in audit if r["selected_by_table13_or_table16_omim"] == "Yes"]
    nominal_not13 = [r for r in audit if r["table16_grail_nominal"] == "Yes"
                     and r["table13_exact_symbol_match"] == "No"]
    require(len({r["gene_symbol"] for r in selected}) == len(selected),
            "Selected set has duplicate symbols")
    require(all(r["gene_symbol"] not in {"", "-"} for r in selected),
            "Selected set has missing symbols")
    require(all(r["table16_grail_nominal"] == "Yes" for r in literal13
                if r["table13_exact_symbol_match"] == "Yes"),
            "Table 13 exact matches disagree with Table 16 GRAIL flags")

    OUTPUT.mkdir(exist_ok=True)
    fields = list(audit[0])
    write_table("supplementary16_all649_evidence_audit.tsv", audit, fields)
    write_table("ng_height_204_genes_from_table16_flags.tsv", selected, fields)
    write_table("ng_height_204_genes_from_table16_flags.csv", selected, fields, ",")
    write_table("ng_height_164_genes_table13_or_omim.tsv", literal13, fields)
    write_table("table16_grail_yes_absent_from_table13_40_genes.tsv", nominal_not13, fields)
    write_text(OUTPUT / "ng_height_204_gene_symbols.txt",
               "".join(r["gene_symbol"] + "\n" for r in sorted(selected, key=lambda r: r["gene_symbol"])))

    classes = Counter(r["table16_support_class"] for r in audit)
    counts = {
        "article_reported_supported_genes": 202,
        "table16_depict_genes": len(audit),
        "table16_grail_yes": classes["Both"] + classes["GRAIL_only"],
        "table16_omim_yes": classes["Both"] + classes["OMIM_only"],
        "table16_both_yes": classes["Both"],
        "table16_grail_only": classes["GRAIL_only"],
        "table16_omim_only": classes["OMIM_only"],
        "table16_grail_or_omim": len(selected),
        "table13_genes": len(table13),
        "table13_exact_symbol_matches_in_table16": sum(r["table13_exact_symbol_match"] == "Yes" for r in audit),
        "table13_or_table16_omim_within_depict": len(literal13),
        "table16_grail_yes_absent_from_table13": len(nominal_not13),
    }
    sources = {}
    for n in (13, 16):
        source = REFERENCE / f"supplementary{n}.xlsx"
        official = REFERENCE / f"ng3097_official_supplementary{n}.xlsx"
        sources[source.name] = {
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "official_url": OFFICIAL_URLS[n],
            "official_copy_byte_identical": source.read_bytes() == official.read_bytes()
            if official.is_file() else None,
            "title": tables[n][1]["A"],
            "columns": tables[n][2 if n == 13 else 3],
        }
    metadata = {
        "doi": "10.1038/ng.3097",
        "audit_date": "2026-09-07",
        "counts": counts,
        "sources": sources,
        "rule": "Table 16 J == Yes AND (K == Yes OR L == Yes)",
        "table13_matching_rule": "Exact original gene symbol; no alias substitution",
        "status": "204 genes extracted; exact article set of 202 unresolved",
        "coordinate_note": "Table 16 locus start/end describe DEPICT loci, not gene-body boundaries",
    }
    write_text(OUTPUT / "extraction_metadata.json", json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(counts, indent=2))
    print("Outputs:", OUTPUT)


if __name__ == "__main__":
    main()
