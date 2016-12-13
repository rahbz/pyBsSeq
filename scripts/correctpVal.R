library(data.table)
library("qvalue")
library(futile.logger)
library("optparse")


option_list = list(
  make_option(c("-i", "--inFile"), type="character", default=NULL, help="Input BED file", metavar="character"),
  make_option(c("-o", "--outFile"), type="character", default="pval-corrected.bed", help="Output file name [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



flog.info("Reading BED file")
bs.bed <- read.table(opt$inFile, as.is = T)


flog.info("Correcting the p-value based on the storey method")
bs.qval <- qvalue(p = bs.bed$V5)

bs.out <- data.frame(chr = bs.bed$V1, start = bs.bed$V2, end = bs.bed$V3, name = bs.bed$V4, qval = bs.qval$qvalues, strand = bs.bed$V6)
write.table(bs.out, file = outFile, quote = F, sep = "\t", row.names = F, col.names = F)
flog.info("Written corrected p-value in the output file")





