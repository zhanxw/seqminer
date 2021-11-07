## NOTE: Need data integrated VCF file (processed by annotateVcf())

a <- proc.time()
library(seqminer)
vcfFile = "chr1.anno.vcf.gz"
geneFile = "../../data/refFlat_hg19.txt.gz"
geneName = ""
annoType = "Nonsynonymous"

fn <- "../../data/rand.chr1.gene.100"
genes = as.matrix(read.table(fn, stringsAsFactors = FALSE))[,1]
for (g in genes) {
  ##print(g)
  mat <- readVCFToMatrixByGene(vcfFile, geneFile, g, annoType)
}

b <- proc.time()

cat("read ", length(genes), " genes.\n")
cat("Elapsed time: ", b - a, " seconds\n")
print(b-a)

