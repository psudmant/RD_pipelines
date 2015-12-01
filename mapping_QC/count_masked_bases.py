import argparse



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", help="Input fasta file")

    args = parser.parse_args()

    countN = {}
    with open(args.fasta, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                contig = line.rstrip()[1:]
                countN[contig] = 0
            else:
                countN[contig] += line.upper().count("N")    

    for contig in sorted(countN):
        print "%s\t%d" % (contig, countN[contig])
