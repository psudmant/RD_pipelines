import os
import argparse
import tables
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import seaborn as sns
sns.set(style="whitegrid")

import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count bases with 0 wssd read depth (ed 0 only)")
    parser.add_argument("--wssd_files", nargs="+", help="List of wssd_out_files to check")
    parser.add_argument("--contigs", nargs="+", help="List of contigs to check")
    parser.add_argument("outfile")
    parser.add_argument("--corrected", action="store_true")

    args = parser.parse_args()

    if args.corrected:
        root = "depthAndStarts_wssd.combined_corrected"
    else:
        root = "depthAndStarts_wssd"

    entries = range(len(args.wssd_files) * (len(args.contigs) + 1))
    wssd_dat = pd.DataFrame(index = entries, 
                            columns = ["SAMPLE", "CONTIG", "NZERO"])

    i = 0
    for fn in args.wssd_files:
        wssd = tables.File(fn, "r")
        sn = os.path.basename(os.path.dirname(fn))
        z_sum = 0
        for contig in args.contigs:
            contig_rd = wssd.getNode("/%s/%s" % (root, contig))
            contig_rd.shape[0]
            sum = contig_rd.shape[0] - np.count_nonzero(contig_rd[:, 0, 0])
            wssd_dat.loc[i] = [sn, contig, sum]
            z_sum += sum
            i += 1
        wssd_dat.loc[i] = [sn, "all", z_sum]
        i += 1
        wssd.close()

    f, ax = plt.subplots(figsize=(6, 0.5*len(args.wssd_files)+1))
    sns.set_color_codes("pastel")
    wssd_dat_summary = wssd_dat.ix[wssd_dat.CONTIG == "all",]
    sns.barplot(x="NZERO", y="SAMPLE", data = wssd_dat_summary.sort("NZERO"))
    ax.set_title("Zero-depth base count QC")
    ticks = ax.get_xticks().tolist()
    tick_labels = map(lambda x: round(x/1000000., 2), ticks)
    plt.xticks(ticks, tick_labels)

    ax.set(ylabel = "", xlabel="Count of 0-depth bases (Mb)")
    #sns.despine(left=True)
    plt.tight_layout()
    plt.savefig(args.outfile)
