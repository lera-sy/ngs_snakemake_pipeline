# NGS Snakemake Pipeline

End-to-end variant calling pipeline for paired tumor/normal WGS data, built on the VSC HPC cluster.

## Features
- FastQC quality control
- BWA-MEM alignment + samtools sorting/indexing
- Joint variant calling with bcftools mpileup
- SnpEff variant annotation
- SQLite database storage of SNPs and effects
- Enformer-based variant scoring (async parallel)
- Summary figures

Enformer scoring parallelised with asyncio + semaphore across 868 variants — reduced wall time from ~3 hours to ~4 minutes.

## Stack
Python · Snakemake · bcftools · SnpEff · SQLite · asyncio
