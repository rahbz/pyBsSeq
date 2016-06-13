args <- commandArgs(TRUE)

library(data.table)
library("qvalue")
library(futile.logger)

inFile=args[1]
outFile=args[2]


getErrorMeth <- function (bs.chr, bs.depth, bs.perMeth, chrName){
  bs.bed.chrName <- which(bs.chr == chrName)
  bs.methPer <- sum(bs.depth[bs.bed.chrName] * bs.perMeth[bs.bed.chrName], na.rm = T)
  bs.totDepth <- sum(bs.depth[bs.bed.chrName], na.rm = T)
  return(bs.methPer/bs.totDepth)
}

flog.info("Reading the input bed file")
bs.bed <- read.table(inFile, as.is = T)
bs.depth <- as.numeric(sapply(bs.bed$V4, function(x){unlist(strsplit(x, split = ":"))[2]}))
flog.info("Calculating the error rates for the chloroplast genome")
error.rate <- getErrorMeth(bs.bed$V1, bs.depth, bs.bed$V5, "chrC")
flog.info("%s is the error rate", error.rate)


flog.info("Getting the p-values for each site separately")
bs.pval <- sapply(seq(1, length(bs.depth)), function(x) {if (bs.depth[x]) {binom.test(as.integer(bs.depth[x] * bs.bed$V5[x]), bs.depth[x], p = error.rate, alternative = "greater")$p.value} else {1}})
bs.qval <- qvalue(p = bs.pval)$qvalues
flog.info("Corrected p-values based on storey method")

flog.info("Writing into BED file")
bs.out <- data.frame(chr = bs.bed$V1, start = bs.bed$V2, end = bs.bed$V3, name = paste(bs.bed$V4, bs.bed$V5, sep = ":"), qval = bs.qval, strand = bs.bed$V6)
write.table(bs.out, file = outFile, quote = F, sep = "\t", row.names = F, col.names = F)



