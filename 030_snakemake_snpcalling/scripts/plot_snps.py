#!/usr/bin/env python
"""
Generate 2 figures from the annotated SNP database.
>    python plot_snps.py <input.sqlite> <fig1_out.svg> <fig2_out.svg>
"""

import argparse
import os
import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


IMPACT_ORDER = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

# Colorblind-safe palette (Wong / IBM): blue + orange are distinguishable
# under deuteranopia, protanopia and tritanopia.
SAMPLE_PALETTE = {"TLE66_N": "#0173B2", "TLE66_T": "#DE8F05"}


def fig1_impact_per_sample(conn, out_path):
    """Grouped bar chart: distinct SNP count per (sample, impact category)."""
    df = pd.read_sql(
        """
        SELECT c.sample, e.impact, COUNT(DISTINCT s.snp_id) AS n
        FROM SNP s
        JOIN Effect e ON s.snp_id = e.snp_id
        JOIN Call   c ON s.snp_id = c.snp_id
        WHERE c.genotype IN ('0/1', '1/1', '1/2', '2/2')
        GROUP BY c.sample, e.impact
        """,
        conn,
    )

    df["impact"] = pd.Categorical(df["impact"], categories=IMPACT_ORDER, ordered=True)
    df = df.sort_values(["impact", "sample"])

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(
        data=df, x="impact", y="n", hue="sample",
        order=IMPACT_ORDER, palette=SAMPLE_PALETTE, ax=ax,
    )
    ax.set_xlabel("Impact category", fontsize=12)
    ax.set_ylabel("Number of distinct SNPs", fontsize=12)
    ax.set_title("SNP impact severity per sample", fontsize=14)
    ax.set_ylim(bottom=0)
    ax.legend(title="Sample", fontsize=11)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def fig2_top_genes(conn, out_path):
    # Horizontal bar of top 15 genes by HIGH/MODERATE variant count:
    # ranking is the natural lens for "which genes deserve a closer look?".
    df = pd.read_sql(
        """
        SELECT gene_name, COUNT(DISTINCT snp_id) AS n_variants
        FROM Effect
        WHERE impact IN ('HIGH', 'MODERATE')
          AND gene_name != ''
        GROUP BY gene_name
        ORDER BY n_variants DESC
        LIMIT 15
        """,
        conn,
    )

    fig, ax = plt.subplots(figsize=(8, 7))
    sns.barplot(
        data=df, y="gene_name", x="n_variants",
        color="#0173B2", ax=ax,
    )
    ax.set_xlabel("Number of HIGH/MODERATE impact variants", fontsize=12)
    ax.set_ylabel("Gene", fontsize=12)
    ax.set_title("Top 15 genes by HIGH/MODERATE impact variants", fontsize=14)
    ax.set_xlim(left=0)
    ax.tick_params(labelsize=11)
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("db",   help="Input SQLite database")
    parser.add_argument("fig1", help="Output path for figure 1 (SVG)")
    parser.add_argument("fig2", help="Output path for figure 2 (SVG)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.fig1), exist_ok=True)
    os.makedirs(os.path.dirname(args.fig2), exist_ok=True)

    conn = sqlite3.connect(args.db)

    print(f"[plot_snps] Generating fig1 → {args.fig1}")
    fig1_impact_per_sample(conn, args.fig1)

    print(f"[plot_snps] Generating fig2 → {args.fig2}")
    fig2_top_genes(conn, args.fig2)

    conn.close()
    print("[plot_snps] Done")


if __name__ == "__main__":
    main()