##http://www.bioconductor.org/help/workflows/variants/
##http://www.bioconductor.org/packages/2.12/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

a <- proc.time()

library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

a.doneLib <- proc.time()

print(a.doneLib-a)
##vcfFile = "vcfRefFlat_index.vcf.gz"
#vcfFile = "/net/fantasia/home/zhanxw/anno/paper/benchmark.vcf/vcfRefFlat_index.gz"
vcfFile <- "ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"

ranges = read.table("../../data/rand.chr1.range.100", stringsAsFactors = FALSE)[,1]

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## ##
## ## some gene may have two ranges, so the returned gene
## getRangeFromName <- function(genesym) {
##     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
##     library(org.Hs.eg.db)

##     ## get GRangesList of all genes
##     txbygene = transcriptsBy(txdb, "gene")

##     # get entrezid from given gene symbol
##     geneid <- select(org.Hs.eg.db, keys=genesym, keytype="SYMBOL", columns="ENTREZID")
##     # get range by matching on the EntrezID
##     rngs <- lapply(geneid$ENTREZID, function(id) range(txbygene[names(txbygene) %in% id]))

##     gnrng <- unlist(do.call(c, rngs), use.names=TRUE)
##     gnrng <- renameSeqlevels(gnrng, c(chr1="1", chr2="2", chr3="3", chr4="4", chr5="5", chr6="6", chr7="7", chr8="8", chr9="9", chr10="10", chr11="11", chr12="12", chr13="13", chr14="14", chr15="15", chr16="16", chr17="17", chr18="18", chr19="19", chr20="20", chr21="21", chr22="22"))

##     names(gnrng) <- sapply(names(gnrng), function(id) {geneid$SYMBOL[ id == geneid$ENTREZID ]} )
##     gnrng
## }

for (j in ranges) {
    cat("Process", j, "\n")
    tmp <- strsplit(j, ":")[[1]]
    chrom <- tmp[1]
    if (grepl("^chr", chrom)) {
        ## strip out "chr" prefix
        chrom <- sub("^chr", "", chrom)
    }
    tmp <- strsplit(tmp[2], "-")[[1]]
    start <- as.integer(tmp[1])
    end <- as.integer(tmp[2])
    
    gnrng <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
    ## gnrng <- getRangeFromName(genes[j])
    if(length(gnrng) == 0) {next}
    param <- ScanVcfParam(which = gnrng, geno = c("GT"))
    ## print('okay till here')
    vcf <- try({
        vcf <- readVcf(vcfFile, "hg19", param)
        vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf))) ## add "chr" prefix
        
        ## find which are nonsynonymous
        nVariant <- dim(vcf)[1]
        print("step 1 - get coding annotation")
        # step 1. get rough annotation category
        # locate by CodingVariants() is not fast...
        loc_all <- locateVariants(vcf, txdb, CodingVariants())
        vcf.coding <- vcf[sort(unique(loc_all$QUERYID)),]
        ## anno <- list()
        ## annoFull <- list()
        ## options(stringsAsFactors = FALSE)
        ## d <- data.frame(anno = as.character(loc_all$LOCATION), idx = loc_all$QUERYID)
        ## annoLevel <- c(spliceSite = 6,
        ##                intron     = 3,
        ##                fiveUTR    = 4,
        ##                threeUTR   = 5,
        ##                coding     = 7,
        ##                intergenic = 1,
        ##                promoter   = 2)
        ## d$level <- factor(d$anno, levels = names(sort(annoLevel)))
        #head(d)
        ## d <- ddply(d, .(idx), function(x) {data.frame(anno = max(as.integer(x$level)), full = paste(unique(x$anno), collapse = ","))})
        #head(d)
        #str(d)
        ## d$anno <- names(annoLevel[match(d$anno, annoLevel)])
        ## #table(d$anno, d$full)
        ## d <- d[order(d$idx), ]
        ## tmp <- data.frame(idx = seq(nVariant),
        ##                   anno = rep(NA, nVariant),
        ##                   full = rep(NA, nVariant))
        ## tmp[d$idx,] <- d
        ## d <- tmp
        ## stopifnot(dim(d)[1] == dim(vcf)[1])

        # step 2. find finer coding annotation
        print("step 2 - fine annotation")
        # library(BSgenome.Hsapiens.UCSC.hg19)
        coding_all <- predictCoding(vcf.coding, txdb, Hsapiens)
        idx <- coding_all$QUERYID
        txt <- coding_all$CONSEQUENCE
        idx <- idx[txt == "nonsynonymous"]

        vcf.nonsyn <- vcf.coding[sort(unique(idx)),]
        geno <- genotypeToSnpMatrix(vcf.nonsyn)
        as(geno$genotypes, "integer")
    },silent=TRUE)
    if(class(vcf)=='try-error') {
        print("something wrong!") ;
        print(vcf);
        print(str(vcf));
        ##break;
    }
}

b <- proc.time()

cat("read ", length(genes), " genes.\n")
cat("Elapsed computing time: ", (b - a.doneLib)["elapsed"], " seconds\n")
cat("Elapsed total time: ", (b - a)["elapsed"], " seconds\n")
print(c("read ", length(genes), " genes.\n"))

print(c("Elapsed computing time: ", sum( (b - a.doneLib)[1:2]), " seconds\n"))
print(c("Elapsed total time: ", sum( (b - a)[1:2]), " seconds\n"))
print("Computation details")
print (b - a.doneLib)
print("Total details")
print (b - a)
