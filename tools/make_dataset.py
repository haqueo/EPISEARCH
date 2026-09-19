#!/usr/bin/env python3
"""Generate a synthetic genotype dataset in the format EPISEARCH expects.

Rows are samples, columns are SNPs coded {0,1,2}, and the last column is the
binary phenotype class. One triple of SNPs is given a genuine 3-way epistatic
effect on the class: none of the three sites is individually predictive, only
the combination is. That triple is your ground truth - a correct search should
rank it at or near the top, and a broken one will not.

  ./tools/make_dataset.py --nsamples 400 --nsnps 60 --out data.txt --seed 7

Prints the planted triple and the exact CLI arguments to stderr.
"""

import argparse
import random
import sys


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsamples", type=int, default=400)
    ap.add_argument("--nsnps", type=int, default=60)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--noise", type=float, default=0.05,
                    help="probability a sample's class label is flipped")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    if args.nsnps < 4:
        print("need at least 4 SNPs", file=sys.stderr)
        return 1

    # Pick the causal triple somewhere in the middle of the matrix.
    triple = sorted(rng.sample(range(args.nsnps), 3))

    rows = []
    for _ in range(args.nsamples):
        snps = [rng.choices([0, 1, 2], weights=[0.25, 0.5, 0.25])[0]
                for _ in range(args.nsnps)]
        # 3-way parity: marginally uninformative, jointly deterministic.
        parity = (snps[triple[0]] + snps[triple[1]] + snps[triple[2]]) % 3
        label = 1 if parity == 0 else 0
        if rng.random() < args.noise:
            label = 1 - label
        rows.append(snps + [label])

    with open(args.out, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(" ".join(str(v) for v in row) + "\n")

    nvars = args.nsnps + 1
    print(f"planted triple (0-indexed SNP columns): {tuple(triple)}", file=sys.stderr)
    print(f"wrote {args.out}: {args.nsamples} samples x {nvars} columns "
          f"({args.nsnps} SNPs + class)", file=sys.stderr)
    print("", file=sys.stderr)
    print(f"  ./episearch -inputfile {args.out} -outputfile out.txt "
          f"-nsamples {args.nsamples} -nvars {nvars} -printall n "
          f"-assoclevel 0.05 -numThread 1", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
