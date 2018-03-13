from __future__ import print_function
from __future__ import division

from Bio import SeqIO

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("contigs_file", help="Fasta index of contigs to include")
    parser.add_argument("outfile")

    args = parser.parse_args()

    contigs = []
    with open(args.contigs_file, "r") as infile:
        for line in infile:
            contigs.append(line.rstrip().split()[0])

    with open(args.outfile, "w") as outfile:
        for record in SeqIO.parse(args.infile, "fasta"):
            if record.id in contigs:
                SeqIO.write(record, outfile, "fasta")
