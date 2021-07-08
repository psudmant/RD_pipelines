import os
import sys
import argparse
import tables
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import seaborn as sns
sns.set(style="whitegrid")

import matplotlib.pyplot as plt

MASKED_COUNT_HG19 = 1683340365

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count bases with 0 wssd read depth (ed 0 only)")
    parser.add_argument("--wssd_files", nargs="+", help="List of wssd_out_files to check")
    parser.add_argument("--stat_files", nargs="+", help="List of tab-delimited files with sample, contig, and 0-depth bases")
    parser.add_argument("--contigs", nargs="+", help="List of contigs to check")
    parser.add_argument("--masked_count", default="hg19_masked_count.tab", help="Tab delimited file with contig and masked base count (Default: %(default)s)")
    parser.add_argument("--out_df", help="Write to tab-delimited file")
    parser.add_argument("--out_pdf", help="Plot pdf")
    parser.add_argument("--corrected", action="store_true")

    args = parser.parse_args()

    if args.wssd_files is None and args.stat_files is None:
        sys.exit("No input files specified. Use --wssd_files or --stat_files.")

    masked_dict = {}
    with open(args.masked_count, "r") as reader:
        for line in reader:
            contig, count = line.rstrip().split()
            masked_dict[contig] = int(count)

    n_masked = 0

    for contig in args.contigs:
        if contig in masked_dict:
            n_masked += masked_dict[contig]
        else:
            print "Contig %s not found in masked contig file %s" % (contig, args.masked_count)

    if args.corrected:
        root = "depthAndStarts_wssd.combined_corrected"
    else:
        root = "depthAndStarts_wssd"

    n_wssd = len(args.wssd_files) if args.wssd_files is not None else 0
    n_stat = len(args.stat_files) if args.stat_files is not None else 0
    entries = range((n_wssd + n_stat) * (len(args.contigs) + 1))
    wssd_dat = pd.DataFrame(index = entries, 
                            columns = ["SAMPLE", "CONTIG", "NZERO"])

    i = 0
    if args.wssd_files is not None:
        for fn in args.wssd_files:
            wssd = tables.File(fn, "r")
            sn = os.path.basename(os.path.dirname(fn))
            print "%s\t%s" % (sn, fn)
            z_sum = 0
            for contig in args.contigs:
                print "\t%s" % contig
                contig_rd = wssd.getNode("/%s/%s" % (root, contig))
                contig_rd.shape[0]
                sum = contig_rd.shape[0] - np.count_nonzero(contig_rd[:, 0, 0])
                wssd_dat.loc[i] = [sn, contig, sum]
                z_sum += sum
                i += 1
            wssd_dat.loc[i] = [sn, "all", z_sum]
            i += 1
            wssd.close()

    stat_dict = {}

    if args.stat_files is not None:
        for fn in args.stat_files:
            with open(fn, "r") as reader:
                for line in reader:
                    sn, contig, nzero = line.rstrip().split()
                    if sn not in stat_dict:
                        stat_dict[sn] = {}
                    stat_dict[sn][contig] = int(nzero)

        for sn in stat_dict.keys():
            tot = 0
            for contig in args.contigs:
                wssd_dat.loc[i] = [sn, contig, stat_dict[sn][contig]]
                tot += stat_dict[sn][contig]
                i += 1
            wssd_dat.loc[i] = [sn, "all", tot]
            i += 1

    if args.out_df is not None:
        wssd_dat.ix[wssd_dat.CONTIG != "all",].to_csv(args.out_df, sep="\t", header=False, index=False)

    if args.out_pdf is not None:
        f, ax = plt.subplots(figsize=(6, 0.5*len(wssd_dat.SAMPLE.unique())+1))
        sns.set_color_codes("pastel")
        wssd_dat_summary = wssd_dat.ix[wssd_dat.CONTIG == "all",]
        wssd_dat_summary.NZERO -= n_masked
        sns.barplot(x="NZERO", y="SAMPLE", data = wssd_dat_summary.sort("NZERO", ascending=False))
        ax.set_title("Zero-depth unmasked base count QC")
        ticks = ax.get_xticks().tolist()
        tick_labels = map(lambda x: round(x/1000000., 1), ticks)
        plt.xticks(ticks, tick_labels)

        ax.set(ylabel = "", xlabel="Count of 0-depth bases (Mb)")
        #plt.axvline(x=n_masked, linewidth=2, color="k", label="hg19 masked bases")
        #sns.despine(left=True)
        plt.tight_layout()
        plt.savefig(args.out_pdf)
