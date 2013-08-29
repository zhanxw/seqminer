suppressPackageStartupMessages(library(seqminer))
fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")

cat("--------------- test readVCFToMatrixByRange ---------------\n") 
try(cfh <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Nonsynonymous"))
try(cfh.nonsyn <- readVCFToMatrixByRange(fileName, "1:196621007-196642234,1:196642235-196716634", "Nonsynonymous"))
try(cfh.syn <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Synonymous"))

cat("--------------- test readVCFToMatrixByGene ---------------\n")
## Extract genotypes
try(cfh <- readVCFToMatrixByGene(fileName, geneFile, "CFH", "Synonymous"))
try(apoe <- readVCFToMatrixByGene(fileName, geneFile, "APOE", ""))
try(ssss <- readVCFToMatrixByGene(fileName, geneFile, "ssss", ""))
## print (ssss)

## Another way to extract from VCF File
cat("--------------- test VCFToListByGene ---------------\n")
try(t <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM", "POS", "ID"), "", ""))
##print (t)

cat("--------------- test VCFToListByGene ---------------\n")
try(t <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM","ID"), c("AC","AN"), ""))
##print (t)

cat("--------------- test VCFToListByGene ---------------\n")
try(t <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM","ID"), c("AC","AN"), c("GT","GQ")))
##print (t)

cat("--------------- test VCFToListByRange ---------------\n")
try(t <- readVCFToListByRange(fileName, "1:196621007-196716634", "Synonymous", c("CHROM","ID", "POS"), c("AC","AN"), c("GT","GQ")))
## print(t)


