#Usage : Rscript plot_GC.R combined_output/wssd_GC_cp.summary GC_summary.pdf
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

input.file <- args[1]
output.file <- args[2]

gc <- read.table(input.file, header=FALSE, sep=" ")
colnames(gc) <- c("indiv", "bases_corrected", "bases_total", "frac_cp2", "type")

# Order samples by GC descending.
gc <- gc[with(gc, order(frac_cp2)),]

# Keep order of sample names based on GC ordering.
gc$indiv <- factor(gc$indiv, levels=unique(gc$indiv))

pdf(output.file, width=6, height=8)
ggplot(gc, aes(x=indiv, y=frac_cp2)) + geom_point(size=1) + coord_flip() + theme_bw() + xlab("Sample") + ylab("Proportion of genome properly called copy 2")+ theme(axis.text.y = element_text(size=6,hjust=1,vjust=0)) + ylim(c(0.75, 1))
dev.off()
