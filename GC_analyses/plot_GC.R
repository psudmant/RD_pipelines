library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
input.file <- args[1]
output.file <- args[2]

data <- read.table(input.file, sep=" ", header=TRUE)

pdf(output.file, width=8, height=6)
ggplot(data, aes(x=GC, y=frac_cp2, colour=indiv)) + geom_line() + theme_bw()
dev.off()
