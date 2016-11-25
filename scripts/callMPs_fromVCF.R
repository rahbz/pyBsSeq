## This script has been written in python and it is much efficient

library("optparse")
library("qvalue")
library(futile.logger)
library(VariantAnnotation)
library(parallel)

option_list = list(
  make_option(c("-i", "--inFile"), type="character", default=NULL, help="Input VCF file name", metavar="character"),
  make_option(c("-o", "--outFile"), type="character", default="out.txt", help="output file name [default= %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


flog.info("Reading VCF file into VariantAnnotation frame")
bs.vcf <- readVcf(opt$inFile, genome="tair10")
flog.info("Done reading VCF file")

flog.info("Calculating the conversion rate at chloroplast genome")
bs.control.ind <- which(seqnames(bs.vcf) == "ChrC")
bs.control.depth <- geno(bs.vcf)$CV[bs.control.ind]
bs.control.perMeth <- geno(bs.vcf)$BT[bs.control.ind]
conv.rate <- 1 - sum(bs.control.perMeth * bs.control.depth, na.rm = T)/ sum(bs.control.depth, na.rm = T)
flog.info(paste("Conversion rate:", conv.rate * 100))


flog.info("Getting MPs")

getMPind <- function(ev, bs.vcf, conv.rate){
  ev.ref <- as.character(ref(bs.vcf)[ev])
  ev.pVal <- numeric()
  if(grepl("[CG]", ev.ref, ignore.case = T, perl = T) & !is.na(geno(bs.vcf)$BT[ev]) & geno(bs.vcf)$CV[ev] > 0 & grepl("[12345]", sub("chr", "", as.character(seqnames(bs.vcf)[ev]), ignore.case = T), perl = T)){
    if (grepl("C", ev.ref, ignore.case = T)){
      ev.strand <- "+"
    } else if (grepl("G", ev.ref, ignore.case = T)){
      ev.strand <- "-"
    }
    ev.depth <- geno(bs.vcf)$CV[ev]
    ev.perMeth <- geno(bs.vcf)$BT[ev]
    ev.pVal <- binom.test(as.integer(ev.perMeth * ev.depth), ev.depth, alternative = "greater", p = conv.rate)$p.value
  }
  return(ev.pVal)
}

vgetMPind <- Vectorize(getMPind)

bs.pVal.all <- numeric()
temp.pVal <- numeric()
chunk.size <- 10000
for (i in seq(1, length(ref(bs.vcf)), chunk.size)){
  start <- i
  end <- i + chunk.size - 1
  temp.pVal <- mclapply(X = start:end, FUN = function(x) {getMPind(x, bs.vcf, conv.rate)}, mc.cores = 16)
  bs.pVal.all <- c(bs.pVal.all, as.numeric(temp.pVal))
  flog.info("Progress: %s", i)
}
#bs.pVal.all <- mclapply(X = 1:length(ref(bs.vcf)), FUN = function(x) {getMPind(x, bs.vcf, conv.rate)}, mc.cores = 16)
#bs.pVal.all <- as.numeric(bs.pVal.all)
bs.pVal.ind <- which(!is.na(bs.pVal.all))
bs.pVal <- as.numeric(bs.pVal.all)[bs.pVal.ind]
bs.refs <- as.character(ref(bs.vcf)[bs.pVal.ind])
bs.strand <- sub("C", "+", sub("G", "-", bs.refs, ignore.case = T), ignore.case = T)

flog.info("Calculating q-values")
bs.qVal <- qvalue(p = bs.pVal)$qvalues
num_meth <- length(which(bs.qVal < 0.01))
flog.info(paste("Number of methylation sites:", num_meth, "out of", length(ref(bs.vcf)), "cytosine sites", sep = " "))

bs.info <- paste(info(bs.vcf)$CX[bs.pVal.ind], geno(bs.vcf)$CV[bs.pVal.ind], geno(bs.vcf)$BT[bs.pVal.ind], sep = ":")
bs.out <- cbind(chr = as.character(seqnames(bs.vcf)[bs.pVal.ind]), start <- start(ranges(bs.vcf)[bs.pVal.ind]) - 1, end = start(ranges(bs.vcf)[bs.pVal.ind]), info = bs.info, score = bs.qVal, strand = bs.strand)


flog.info("Writing into output file")
write.table(bs.out, file=opt$outFile, sep = "\t", quote = F, row.names = F, col.names = F)




warnings()





