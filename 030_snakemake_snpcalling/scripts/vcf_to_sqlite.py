#!/usr/bin/env python
"""
Load annotated VCF (with SnpEff ANN field) into a SQLite database
w three tables: SNP, Effect, Call.
>    python vcf_to_sqlite.py <input.vcf> <output.sqlite>
"""

import argparse
import os
import sqlite3

import pandas as pd
import vcfpy


ANN_FIELDS = [
    "allele", "effect", "impact", "gene_name", "gene_id",
    "feature_type", "feature_id", "biotype", "rank",
    "hgvs_c", "hgvs_p", "cdna_pos", "cds_pos", "aa_pos",
    "distance", "errors",
]


def parse_vcf(vcf_path):
    """Read the VCF and return (snp_df, effect_df, call_df)."""
    snps, effects, calls = [], [], []

    reader = vcfpy.Reader.from_path(vcf_path)
    for record in reader:
        chrom = record.CHROM
        pos   = record.POS
        ref   = record.REF
        alt   = record.ALT[0].value          # biallelic after vt decompose
        snp_id = f"{chrom}_{pos}_{ref}_{alt}"

        snps.append({
            "snp_id": snp_id,
            "chrom":  chrom,
            "pos":    pos,
            "ref":    ref,
            "alt":    alt,
            "qual":   record.QUAL,
            "filter": ";".join(record.FILTER) if record.FILTER else "",
        })

        for call in record.calls:
            calls.append({
                "snp_id":   snp_id,
                "sample":   call.sample,
                "genotype": call.data.get("GT"),
            })

        for ann_str in record.INFO.get("ANN", []):
            parts = ann_str.split("|")
            parts += [""] * (16 - len(parts))   # defensive pad
            row = {"snp_id": snp_id}
            row.update(dict(zip(ANN_FIELDS, parts)))
            effects.append(row)

    reader.close()
    return pd.DataFrame(snps), pd.DataFrame(effects), pd.DataFrame(calls)


def write_db(snp_df, effect_df, call_df, db_path):
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    snp_df.to_sql("SNP",    conn, if_exists="replace", index=False)
    effect_df.to_sql("Effect", conn, if_exists="replace", index=False)
    call_df.to_sql("Call",   conn, if_exists="replace", index=False)
    conn.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="Input annotated VCF file")
    parser.add_argument("db",  help="Output SQLite database path")
    args = parser.parse_args()

    print(f"[vcf_to_sqlite] Reading {args.vcf}")
    snp_df, effect_df, call_df = parse_vcf(args.vcf)
    print(f"[vcf_to_sqlite] Parsed {len(snp_df)} SNPs, "
          f"{len(effect_df)} effects, {len(call_df)} calls")

    print(f"[vcf_to_sqlite] Writing {args.db}")
    write_db(snp_df, effect_df, call_df, args.db)
    print("[vcf_to_sqlite] Done")


if __name__ == "__main__":
    main()