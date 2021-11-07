a <- proc.time()

library(seqminer)
vcf <- "ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
gene <- "knownGene.txt.gz" 
ref <- "human.g1k.v37.fa"

param <- list(reference = ref,
              geneFile = gene,
              geneFileFormat = "knownGene"
              )
param <- makeAnnotationParameter(param)
inVcf <- vcf
outVcf <- paste0(getwd(), "/", "out.vcf.gz")
annotateVcf (inVcf, outVcf, param)


b <- proc.time()

print(c("Elapsed total time: ", sum( (b - a)[1:2]), " seconds\n"))
print("Total details")
print (b - a)
