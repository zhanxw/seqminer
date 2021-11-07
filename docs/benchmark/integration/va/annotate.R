default.args <- c(system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation"), "tmp.vcf", "1000")
#default.args <- c("ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz", "anno.chr1.vcf", "1000")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    args <- default.args
}


##http://www.bioconductor.org/help/workflows/variants/
##http://www.bioconductor.org/packages/2.12/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf


## NOTE: VariantAnnotation cannot read too many chromosomes. Need to be iterative
## E.g. R will consume 32Gb and behave strangely (report alloc problems)

a <- proc.time()

setwd("~/anno/paper/benchmark.va/annotate")
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#library(SNPlocs.Hsapiens.dbSNP.20120608) ## dbSNP 137
library(SNPlocs.Hsapiens.dbSNP.20100427) ## dbSNP 131
source("loc2rsid.R", echo = TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)

library(PolyPhen.Hsapiens.dbSNP131)

library(plyr)

a.doneLib <- proc.time()

packageVersion("VariantAnnotation")
fl <- args[1]
## fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
print(fl)
## # read all at once
## system.time(vcf <- readVcf(fl, "hg19"))
## print(vcf)
## print(dim(vcf))

# read by chunks
yieldSize <- as.integer(args[3])
print(yieldSize)
tab <- TabixFile(fl, yieldSize = yieldSize)
open(tab)
out.fn <- args[2]
print(out.fn)
if (file.exists(out.fn)) {
    file.remove(out.fn)
}
con.out <- file(out.fn, open="a")
total <- 0
if (FALSE) {
    while(TRUE) {
        chunk <- readVcf(tab, "hg19")
        total <- total + dim(chunk)[1]
        if (nrow(chunk) == 0) {
            break
        } else {
            vcf <- chunk
        }
    }
}

while(TRUE) {
    chunk <- readVcf(tab, "hg19")
    total <- total + dim(chunk)[1]
    cat("read ", total, "records\n")
    ## if (total != 6000) next
    if (nrow(chunk) == 0) {
        break
    } else {
        vcf <- chunk

        ## append "chr" prefix
        ## head(seqlevels(vcf))
        ## head(seqlevels(txdb))
        ## intersect(seqlevels(vcf), seqlevels(txdb))
        vcf.nochr <- vcf
        vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf))) ## add "chr" prefix

        ## annotate variants
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

        getAnnotation <- function(vcf) {
            nVariant <- dim(vcf)[1]
            print("step 1 - rough annotation")
            # step 1. get rough annotation category
            # locate by AllVariants() is not fast...
            loc_all <- locateVariants(vcf, txdb, AllVariants())
            anno <- list()
            annoFull <- list()
            options(stringsAsFactors = FALSE)
            d <- data.frame(anno = as.character(loc_all$LOCATION), idx = loc_all$QUERYID)
            annoLevel <- c(spliceSite = 6,
                           intron     = 3,
                           fiveUTR    = 4,
                           threeUTR   = 5,
                           coding     = 7,
                           intergenic = 1,
                           promoter   = 2)
            d$level <- factor(d$anno, levels = names(sort(annoLevel)))
            #head(d)
            d <- ddply(d, .(idx), function(x) {data.frame(anno = max(as.integer(x$level)), full = paste(unique(x$anno), collapse = ","))})
            #head(d)
            #str(d)
            d$anno <- names(annoLevel[match(d$anno, annoLevel)])
            #table(d$anno, d$full)
            d <- d[order(d$idx), ]
            tmp <- data.frame(idx = seq(nVariant),
                              anno = rep(NA, nVariant),
                              full = rep(NA, nVariant))
            tmp[d$idx,] <- d
            d <- tmp
            stopifnot(dim(d)[1] == dim(vcf)[1])
            
            # step 2. find finer coding annotation
            print("step 2 - fine annotation")
            # library(BSgenome.Hsapiens.UCSC.hg19)
            coding_all <- predictCoding(vcf, txdb, Hsapiens)
            idx <- coding_all$QUERYID
            txt <- coding_all$CONSEQUENCE
            for (i in seq_len(length(coding_all))) {
                d$anno[idx[i]] <- sub("coding", txt[i], d$anno[idx[i]])
                d$full[idx[i]] <- sub("coding", txt[i], d$full[idx[i]])
            }
            head(d)

            print("step 3 - polyphen")
            # save(list = ls(), file = sprintf("poly.%d.Rdata", yieldSize))
            # step 3. add polyphen score
            ## as Polyphen can only be queried by rsid, need to get rsid first
            locs <- rowData(vcf) #GRanges(chrs, IRanges(as.integer(posi[,2]), width=1))
            idx <- which(width(locs) == 1)
            snp.locs <- locs[idx,]
            snp.locs <- renameSeqlevels(snp.locs, sub("chr", "ch", seqlevels(snp.locs))) ## change "chr" prefix to "ch"
            rsid <- loc2rsid(snp.locs)
            rsid <- do.call(c, lapply(rsid, function(x){
                if (length(x) == 0) {return(NA)}
                return(x[1])}))
            d$rsid <- NA
            d$rsid[idx] <- rsid
            ## also get REF, and first ALT
            d$ref <- as.character(locs$REF)
            alt <- do.call(c, lapply(locs$ALT, function(x){
                return(x[1])}))
            d$alt <- as.character(alt)
            d$key <- paste(d$rsid, d$ref, d$alt, sep = "_")
            ## get Polyphen
            pp <- select(PolyPhen.Hsapiens.dbSNP131, keys=d$rsid, 
                         columns=c("NT1", "NT2", "TRAININGSET", "PREDICTION", "PPH2PROB"))
            pp$key <- paste(pp$RSID, pp$NT1, pp$NT2, sep = "_")
            idx <- match(d$rsid, pp$RSID)
            d$polyphen <- NA
            d$polyphen <- pp$PREDICTION[idx]
            
            d
        }
        anno <- getAnnotation(vcf)

        ### try to write vcf back
        info.tag <- rownames(info(header(vcf)))
        info.tag <- c(info.tag, "ANNO", "ANNOFULL", "Polyphen")
        info <- rbind(info(header(vcf)),
                      DataFrame(Number = "1", Type = "String", Description = "Annotation"),
                      DataFrame(Number = "1", Type = "String", Description = "Fulll Annotation"),
                      DataFrame(Number = "1", Type = "String", Description = "Polyphen"))
        rownames(info) <- info.tag

        old.hdr <- header(vcf)
        old.hdr.info <- info(old.hdr)
        old.hdr@header@listData$INFO <- info
        new.hdr <- old.hdr

        info(vcf)$ANNO <- anno$anno
        info(vcf)$ANNOFULL <- anno$full
        vcf@exptData@listData$header <- new.hdr

        ## strip "chr" prefix
        vcf <- renameSeqlevels(vcf, seqlevels(vcf.nochr)) ## add "chr" prefix

        writeVcf(vcf, con.out)
    }
}
close(tab)
close(con.out)

b <- proc.time()

print(c("Elapsed computing time: ", sum( (b - a.doneLib)[1:2]), " seconds\n"))
print(c("Elapsed total time: ", sum( (b - a)[1:2]), " seconds\n"))
print("Computation details")
print (b - a.doneLib)
print("Total details")
print (b - a)
