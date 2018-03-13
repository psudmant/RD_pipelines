from __future__ import print_function
from __future__ import division

from Bio import SeqIO

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")

    args = parser.parse_args()

with open(args.outfile, "w") as outfile:
    for record in SeqIO.parse(args.infile, "fasta"):
        if len(record.seq) - record.seq.count("N") > 10000:
            SeqIO.write(record, outfile, "fasta")
