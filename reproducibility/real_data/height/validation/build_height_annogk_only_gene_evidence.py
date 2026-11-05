#!/usr/bin/env python3
"""Validate curated literature examples against frozen region/gene inputs."""
import csv
import hashlib
import json
from pathlib import Path

BASE = Path(repro_path("height_validation_workspace"))
DATA = BASE / "height_ten_seed_region_stability_q0p1"
CURATED = BASE / "height_annogk_only_gene_evidence_2026-09-14.json"
GENES = BASE / "reference/genecards_height_top1000_hg19_gene_positions.csv"
OUT = BASE / "height_annogk_only_gene_evidence_2026-09-14"


def read(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def tex(value):
    replacements = {"&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#"}
    return "".join(replacements.get(c, c) for c in str(value))


def main():
    curated = json.loads(CURATED.read_text())
    genes = read(GENES, ",")
    inputs = [CURATED, GENES]
    regions = {}
    for ancestry in ["EUR", "AFR"]:
        path = DATA / f"{ancestry}_all10_union_regions.tsv"
        inputs.append(path)
        regions[ancestry] = read(path)
        expected = 105 if ancestry == "EUR" else 12
        assert sum(r["all_10_class"] == "AnnoGK_only" for r in regions[ancestry]) == expected
    rows = []
    for item in curated:
        for ancestry in item["ancestries"]:
            matches = []
            for g in genes:
                if g["genecards_symbol"] != item["gene"] or g["autosomal_1_22"] != "TRUE":
                    continue
                assert g["coordinate_system"] == "GRCh37/hg19"
                assert g["interval_convention"] == "1-based closed"
                for r in regions[ancestry]:
                    if r["all_10_class"] != "AnnoGK_only":
                        continue
                    if (g["chromosome"] == r["chromosome"] and
                            int(r["start"]) <= int(g["start_1based"]) and
                            int(g["end_1based"]) <= int(r["end"])):
                        assert int(r["annogk_n_seeds"]) == 10 and int(r["ghost_n_seeds"]) < 10
                        assert float(r["annogk_worst_detected_min_q"]) <= .1
                        matches.append(dict(
                            ancestry=ancestry, region_id=r["region_id"], gene=item["gene"],
                            annogk_n_seeds=int(r["annogk_n_seeds"]),
                            ghost_n_seeds=int(r["ghost_n_seeds"]),
                            ensembl_gene_id=g["ensembl_gene_id"], genecards_rank=g["genecards_rank"],
                            gene_start=g["start_1based"], gene_end=g["end_1based"],
                            evidence_type=item["evidence_type"], study=item["study"],
                            title=item["title"], url=item["url"],
                            evidence_zh=item["evidence_zh"], function_zh=item["function_zh"],
                            evidence_en=item["evidence_en"], function_en=item["function_en"],
                            function_url=item["function_url"]))
            assert len(matches) == 1, (item["gene"], ancestry, len(matches))
            rows.extend(matches)
    assert len(rows) == 9 and len({r['gene'] for r in rows}) == 8
    OUT.mkdir(exist_ok=True)
    with (OUT / "evidence_table.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    intro = ["# EUR/AFR AnnoGK 特有稳定 regions 中的身高相关基因：代表性文献证据", "",
             "核查日期：2026-09-14。", "",
             "范围：M=1，seeds=100,200,...,1000，q≤0.1；AnnoGK 10/10 seeds 检出，GhostKnockoff <10/10。使用原有 EUR 定义的 outer risk regions（AFR 也沿用这些区间），以及 GENCODE v19 GRCh37/hg19 完整基因体包含规则。", "",
             "这是从冻结 GeneCards Top-1000 覆盖基因中核查的 8 个代表例子（9 个 ancestry-region-gene 条目），不是所有基因的系统文献穷尽或显著性检验。每个例子的完整基因体均直接重新匹配到相应区间。文献给出基因/位点与身高或生长障碍的既有联系，并不表示本次选中的 SNP 已被验证通过该基因起作用。", "",
             "| 人群 | Region（hg19） | Gene | Seeds：AnnoGK / Ghost | 研究及证据 | 功能 |",
             "|---|---|---|---|---|---|"]
    for r in rows:
        study = f"[{r['study']}]({r['url']})：{r['evidence_zh']}"
        function = r['function_zh']
        if r['function_url'] != r['url']:
            function += f" [功能实验]({r['function_url']})"
        intro.append(f"| {r['ancestry']} | {r['region_id']} | **{r['gene']}** | {r['annogk_n_seeds']}/10 / {r['ghost_n_seeds']}/10 | {study} | {function} |")
    intro += ["", "解释与写作建议：", "",
              "- IGF1R 所在 EUR region 为 10/10 对 0/10，可以明确描述为在这十个 seeds 中 GhostKnockoff 从未检出；其他例子应描述为未达到十个 seeds 均检出的稳定性标准。",
              "- GNAS 的同一区间在 EUR、AFR 中均为 AnnoGK 特有稳定 region，但 GhostKnockoff 分别为 9/10 和 8/10。",
              "- HHIP、EFEMP1 的证据是成人身高 GWAS 位点关联；其候选效应基因身份不能仅由该 GWAS 确定。其他条目主要是人类罕见变异/生长障碍证据，不能与一般人群身高 GWAS 证据混称。",
              "- 外部研究的受试人群未必是 AFR。位于 AFR 检出 region 内不构成该文献在 AFR 中独立复现的证明。",
              "- 宽区间内含多个基因；这里没有进一步确认 selected SNP 是否在该基因体内、与文献变异同一 LD 信号，或作用于该基因。因此建议称为生物学合理性支持，不称为新因果基因发现或 FDR 证明。",
              "", "Suggested manuscript text:", "",
              "> Regions reproducibly detected by AnnoGK in all ten knockoff realizations, but not meeting this criterion for GhostKnockoff, contained genes with established roles in human growth. Examples included IGF1R, FGFR3, PAPPA2 and IHH in EUR, and DMP1 and GNAS in AFR. The HHIP and EFEMP1 regions also overlapped previously reported adult-height GWAS loci. These observations support the biological plausibility of the additional stable regions, without assigning causality to individual genes or selected variants.", "",
              "Verification: nine unique ancestry-region-gene mappings; full containment reconstructed from frozen coordinate inputs; original seed recurrence and q-threshold fields verified. The pre-existing reference-panel build audit is not newly resolved by this literature lookup.", "",
              "Reproduce from repository root: `Rscript reproducibility/run.R real_data/height/validation/build_height_annogk_only_gene_evidence.py`.", ""]
    (OUT / "report.md").write_text("\n".join(intro))

    lines = [r"% Requires longtable, booktabs, hyperref; compile with XeLaTeX/LuaLaTeX for Unicode names.",
             r"\begin{longtable}{p{0.055\textwidth}p{0.19\textwidth}p{0.065\textwidth}p{0.07\textwidth}p{0.32\textwidth}p{0.17\textwidth}}",
             r"\caption{Representative growth-related genes fully contained in AnnoGK-specific stable regions.}\label{tab:height-known-genes}\\",
             r"\toprule", r"Ancestry & Region (hg19) & Gene & Seeds A/G & Published evidence & Function \\",
             r"\midrule", r"\endfirsthead", r"\toprule",
             r"Ancestry & Region (hg19) & Gene & Seeds A/G & Published evidence & Function \\",
             r"\midrule", r"\endhead"]
    for r in rows:
        fields = [r['ancestry'], tex(r['region_id']).replace('-', r'\allowbreak-'), tex(r['gene']),
                  f"10/{r['ghost_n_seeds']}",
                  r"\href{" + r['url'] + "}{" + tex(r['study']) + "}. " + tex(r['evidence_en']),
                  tex(r['function_en'])]
        if r['function_url'] != r['url']:
            fields[-1] += r" \href{" + r['function_url'] + "}{Functional study}."
        lines.append(" & ".join(fields) + r" \\")
    lines += [r"\bottomrule", r"\end{longtable}",
              r"\noindent A/G denotes the number of seeds with region detection for AnnoGK/GhostKnockoff, each out of ten ($M=1$, $q\leq0.1$). AnnoGK-specific stable regions satisfy A=10 and G<10. Genes are fully contained using GENCODE v19 coordinates. Human growth-disorder evidence and common-variant GWAS locus associations are distinguished; these examples do not establish causal genes for the selected variants.", ""]
    (OUT / "table.tex").write_text("\n".join(lines))
    metadata = dict(date="2026-09-14", examples=len(rows), distinct_genes=8,
                    selection="Representative literature-supported genes, not exhaustive",
                    sha256={str(p.relative_to(BASE)) if p.is_relative_to(BASE) else p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in inputs + [Path(__file__)]})
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Validated {len(rows)} ancestry-region-gene entries. Wrote {OUT}.")


if __name__ == "__main__":
    main()
