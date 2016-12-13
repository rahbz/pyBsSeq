#make genome matrix from all the BED files generated  
# Get the list of all the BED files present 

library("optparse")


option_list = list(
  make_option(c("-i", "--inDIR"), type="character", default=NULL, help="Path to the directory with BED files", metavar="character"),
  make_option(c("-o", "--outFile"), type="character", default="pval-corrected.bed", help="Output file name [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



allScoreFiles <- list.files(opt$inDir, pattern = "[.]pVal.bed$")

sed 's/:/\t/g' $1 | awk '{gsub("+", "C", $8);gsub("-", "G", $8);print $1 "\t" $2 "\t" $4 "\t" $8 "\t" "40/" $6 "/" int($5 * $6) "/" $5 - int($6 * $5) }' > $2



