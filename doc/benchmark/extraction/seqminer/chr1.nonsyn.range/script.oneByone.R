a <- proc.time()
library(seqminer)
vcfFile = "chr1.anno.vcf.gz"
geneFile = "../../data/refFlat_hg19.txt.gz"
geneName = ""
annoType = "Nonsynonymous"

fn <- "../../data/rand.chr1.range.100"
ranges = read.table(fn, stringsAsFactors = FALSE)[,1]
for (r in ranges) {
  ##print(g)
  mat <- readVCFToMatrixByRange(vcfFile, r, annoType)
}

b <- proc.time()

cat("read ", length(genes), " genes.\n")
cat("Elapsed time: ", b - a, " seconds\n")
print(b-a)

