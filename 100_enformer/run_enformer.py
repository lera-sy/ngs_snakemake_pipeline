#
import asyncio
import csv
import sys
import time

import vcfpy
from fake_enformer import async_predict


def extract_coordinates(vcf_path):
    """Pull true SNP coordinates in fake_enformer's hg38:chrN:pos:REF:ALT format."""
    coords = []
    for record in vcfpy.Reader.from_path(vcf_path):
        ref = record.REF
        # CHROM should be 'chr9'; prepend 'chr' defensively in case it's bare '9'
        chrom = record.CHROM if record.CHROM.startswith("chr") else f"chr{record.CHROM}"
        for alt in record.ALT:
            a = alt.value
            # keep only single-base ACGT substitutions; indels would be rejected (HTTP 500)
            if len(ref) == 1 and len(a) == 1 and ref in "ACGT" and a in "ACGT":
                coords.append(f"hg38:{chrom}:{record.POS}:{ref}:{a}")
    return coords


async def score_all(coords, concurrency):
    sem = asyncio.Semaphore(concurrency)
    results = {}

    async def worker(coord):
        async with sem:                       # at most `concurrency` calls in flight
            try:
                results[coord] = await async_predict(coord)
            except Exception as e:             # one bad coord shouldn't abort the rest
                print(f"FAILED {coord}: {e}", file=sys.stderr)

    await asyncio.gather(*(worker(c) for c in coords))
    return results


def main():
    vcf_path = sys.argv[1]
    concurrency = int(sys.argv[2]) if len(sys.argv) > 2 else 20  # optional 2nd arg to experiment

    coords = extract_coordinates(vcf_path)
    unique = list(dict.fromkeys(coords))       # dedupe, preserve order = the caching win
    print(f"{len(coords)} SNPs, {len(unique)} unique coords, concurrency={concurrency}",
          file=sys.stderr)

    start = time.time()
    results = asyncio.run(score_all(unique, concurrency))
    elapsed = time.time() - start
    print(f"Scored {len(results)}/{len(unique)} coords in {elapsed:.1f}s "
          f"(concurrency={concurrency})", file=sys.stderr)

    with open("enformer_scores.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["coordinate", "score"])
        for c in unique:
            if c in results:
                w.writerow([c, results[c]])


if __name__ == "__main__":
    main()