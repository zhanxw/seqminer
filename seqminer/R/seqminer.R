#' Efficiently Read Sequencing Data (VCF format, METAL format) into R
#'
#' SeqMiner provides functions to easily load Variant Call Format (VCF) or METAL format into R
#'
#' The aim of this package is to save your time parsing large text file. 
#' That means data processing time can be saved for other researches.
#' This packages requires Bgzip compressed and Tabix indexed files as input.
#' If input files contaians annotation by TabAnno (),
#' it is possible to extract information at the unit of genes.
#' 
#' @docType package
#' @name SeqMiner
#' @useDynLib seqminer
#' @importFrom utils
NULL


#' Check input file has tabix index
#' 
#' @param fileName character vector, representing file names an input file name
#' @return TRUE if an index file with ".tbi" exists
#' @keywords internal
hasIndex <- function(fileName) {
  if (file.exists(paste(fileName, ".tbi", sep = ""))) {
    return (TRUE)
  } else {
    print("Cannot find index file (you can create it using: [ tabix -p vcf FILE_NAME.vcf.gz ])")  
    return (FALSE)
  }
}

#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param range character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @return genotype matrix
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' cfh <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Nonsynonymous")
readVCFToMatrixByRange <- function(fileName, range, annoType) {
  stopifnot(file.exists(fileName), length(fileName) == 1, hasIndex(fileName))
  storage.mode(fileName) <- "character"
  storage.mode(range)    <- "character"
  storage.mode(annoType) <- "character"
  .Call("readVCFToMatrixByRange", fileName, range, annoType, PACKAGE="seqminer");
};

#' Read a gene from VCF file and return a genotypes matrix
#'
#' @param fileName charactr, represents an input VCF file (Bgzipped, with Tabix index)
#' @param geneFile character, a text file listing all genes in refFlat format
#' @param geneName character vector, which gene(s) to be extracted
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @return genotype matrix
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- readVCFToMatrixByGene(fileName, geneFile, "CFH", "Synonymous")
readVCFToMatrixByGene <- function(fileName, geneFile, geneName, annoType) {
  stopifnot(file.exists(fileName), length(fileName) == 1)
  stopifnot(file.exists(geneFile), length(geneFile) == 1)

  storage.mode(fileName) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  storage.mode(annoType) <- "character"
  .Call("readVCFToMatrixByGene", fileName, geneFile, geneName, annoType, PACKAGE="seqminer");
};

#' Read information from VCF file in a given range and return a list
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param range character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @param vcfColumn character vector, which vcf columns to extract. It can be chosen from CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT and etc.
#' @param vcfInfo character vector, which should be tags in the INFO columns to extarct. Common choices include: DP, AC, AF, NS 
#' @param vcfIndv character vector, which values to extract at individual level. Common choices are: GT, GQ, GD
#' @return a list of genes, and each elements has specified vcfColumn, vcfinfo, vcfIndv
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' cfh <- readVCFToListByRange(fileName, "1:196621007-196716634", "Nonsynonymous", 
#'                             c("CHROM", "POS"), c("AF", "AC"), c("GT") )
readVCFToListByRange <- function(fileName, range, annoType, vcfColumn, vcfInfo, vcfIndv) {
  stopifnot(file.exists(fileName), length(fileName) == 1, hasIndex(fileName))
  storage.mode(fileName) <- "character"
  storage.mode(range)    <- "character"
  storage.mode(annoType) <- "character"
  storage.mode(vcfColumn)<- "character"
  storage.mode(vcfInfo)  <- "character"
  storage.mode(vcfIndv)  <- "character"
  .Call("readVCFToListByRange", fileName, range, annoType, vcfColumn, vcfInfo, vcfIndv, PACKAGE="seqminer");
};

#' Read information from VCF file in a given range and return a list
#'
#' @param fileName character, represents an input VCF file (Bgzipped, with Tabix index)
#' @param geneFile character, a text file listing all genes in refFlat format
#' @param geneName character vector, which gene(s) to be extracted
#' @param annoType character, annotated types you would like to extract, such as "Nonsynonymous", "Synonymous". This can be left empty.
#' @param vcfColumn character vector, which vcf columns to extract. It can be chosen from CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT and etc.
#' @param vcfInfo character vector, which should be tags in the INFO columns to extarct. Common choices include: DP, AC, AF, NS 
#' @param vcfIndv character vector, which values to extract at individual level. Common choices are: GT, GQ, GD
#' @return a list of genes, and each elements has specified vcfColumn, vcfinfo, vcfIndv
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", 
#'                            c("CHROM", "POS"), c("AF", "AC"), c("GT") )
readVCFToListByGene <- function(fileName, geneFile, geneName, annoType, vcfColumn, vcfInfo, vcfIndv) {
  stopifnot(file.exists(fileName), length(fileName) == 1, hasIndex(fileName))
  stopifnot(file.exists(geneFile), length(geneFile) == 1)
  storage.mode(fileName) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  storage.mode(annoType) <- "character"
  storage.mode(vcfColumn)<- "character"
  storage.mode(vcfInfo)  <- "character"
  storage.mode(vcfIndv)  <- "character"
  .Call("readVCFToListByGene", fileName, geneFile, geneName, annoType, vcfColumn, vcfInfo, vcfIndv, PACKAGE="seqminer");
};

#' Read association statistics by gene from METAL-format files. Both score statistics and covariance statistics will be extracted.
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @param covFiles character vector, covaraite files (rvtests outputs using --meta cov)
#' @param geneFile character, a text file listing all genes in refFlat format
#' @param geneName character vector, which gene(s) to be extracted
#' @return a list of statistics including chromosome, position, allele frequency, score statistics, covariance and annotation(if input files are annotated).
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- rvmeta.readDataByGene(scoreFileName, covFileName, geneFile, "CFH")
rvmeta.readDataByGene <- function(scoreTestFiles, covFiles, geneFile, geneName) {
  stopifnot(file.exists(scoreTestFiles))
  stopifnot(is.null(covFiles) || (file.exists(covFiles) && length(covFiles) == length(scoreTestFiles)))
  stopifnot(file.exists(geneFile), length(geneFile) == 1)
  storage.mode(scoreTestFiles) <- "character"
  storage.mode(covFiles) <- "character"
  storage.mode(geneFile) <- "character"
  storage.mode(geneName) <- "character"
  if (is.null(covFiles)) {
    .Call("rvMetaReadDataByGene", scoreTestFiles, "", geneFile, geneName, PACKAGE="seqminer");
  } else {
    .Call("rvMetaReadDataByGene", scoreTestFiles, covFiles, geneFile, geneName, PACKAGE="seqminer");
  }
};

#' Read association statistics by range from METAL-format files. Both score statistics and covariance statistics will be extracted.
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @param covFiles character vector, covaraite files (rvtests outputs using --meta cov)
#' @param ranges character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return a list of statistics including chromosome, position, allele frequency, score statistics, covariance and annotation(if input files are annotated).
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
#' cfh <- rvmeta.readDataByRange(scoreFileName, covFileName, "1:196621007-196716634")
rvmeta.readDataByRange <- function (scoreTestFiles, covFiles, ranges) {
  stopifnot(file.exists(scoreTestFiles))
  stopifnot(is.null(covFiles) || (file.exists(covFiles) && length(covFiles) == length(scoreTestFiles)))
  storage.mode(scoreTestFiles) <- "character"
  storage.mode(covFiles) <- "character"
  storage.mode(ranges) <- "character"
  if (is.null(covFiles)) {
    .Call("rvMetaReadDataByRange", scoreTestFiles, "",
          ranges, PACKAGE = "seqminer")
  } else {
    .Call("rvMetaReadDataByRange", scoreTestFiles, covFiles,
          ranges, PACKAGE = "seqminer")
  }
}

#' Read covariance by range from METAL-format files. 
#'
#' @param covFile character, a covariance file (rvtests outputs using --meta cov)
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return a matrix of covariance within given range
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
#' cfh <- rvmeta.readCovByRange(covFileName, "1:196621007-196716634")
rvmeta.readCovByRange <- function(covFile, tabixRange) {
  stopifnot(file.exists(covFile))
  storage.mode(covFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readCovByRange", covFile, tabixRange, PACKAGE="seqminer");
};

#' Read score test statistics by range from METAL-format files. 
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return score test statistics within given range
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
#' cfh <- rvmeta.readScoreByRange(scoreFileName, "1:196621007-196716634")
rvmeta.readScoreByRange <- function(scoreTestFiles, tabixRange) {
  stopifnot(file.exists(scoreTestFiles))
  storage.mode(scoreTestFiles) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readScoreByRange", scoreTestFiles, tabixRange, PACKAGE="seqminer");
};

#' Read tabix file, similar to running tabix in command line.
#'
#' @param tabixFile character, an tabix indexed file
#' @param tabixRange  character, a text indicating which range in the VCF file to extract. e.g. 1:100-200
#' @return character vector, each elements is an individual line
#' @export
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
#' snp <- tabix(fileName, "1:196621007-196621007")
tabix <- function(tabixFile, tabixRange) {
  stopifnot(file.exists(tabixFile))
  storage.mode(tabixFile) <- "character"
  storage.mode(tabixRange) <- "character"
  .Call("readTabixByRange", tabixFile, tabixRange, PACKAGE="seqminer");
}

.onAttach <- function(libname, pkgname){
  newVersionLink = "http://zhanxw.com:8080/seqminer/version"
  conn <- url(newVersionLink)
  ret <- tryCatch(readLines(conn, n = 2), error = function(e) {NULL})
  close(conn)

  if (!is.null(ret) && length(ret) == 2) {
    version <- ret[1]
    if (utils::packageVersion("seqminer") < version) {
      if (length(ret) > 1) {
        packageStartupMessage(ret[2])
      } else {
        packageStartupMessage("Found new version of seqminer: ", ret)
      }
    }
  }
  ## packageStartupMessage("seqminer loaded ...")
}

#' Read null model statistics
#'
#' @param scoreTestFiles character vector, score test output files (rvtests outputs using --meta score)
#' @return a list of statistics fitted under the null mode (without genetic effects)
#' @export
#' @import stringr
#' @seealso http://zhanxw.com/seqminer/ for online manual and examples
#' @examples
#' scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
rvmeta.readNullModel <- function(scoreTestFiles) {
  stopifnot(file.exists(scoreTestFiles))
  read.null.model <- function(fn) {
      ## cat("gzFile, fn = ", fn, "\n")
      f <- gzfile(fn, "rb")
      ret = list()
      beginStore = FALSE
      while(TRUE) {
          suppressWarnings(p <- readLines(con = f, n = 1, warn = FALSE))
          if (substr(p, 1, 1) != "#") {
              break
          }
          #print(p)
          if (!beginStore) {
              if ( p[1] == "##NullModelEstimates") {
                  beginStore = TRUE
                  next
              }
          } else {
              ret[length(ret) + 1] <- p
          }
      }
      close(f)
      ret
  }
  model.to.matrix <- function(ret) {
      if (is.null(ret) || length(ret) < 1) {
          return(numeric(0))
      }
      ## library(stringr)
      ## ret <- read.null.model(ret)
      ret <- lapply(ret, function(x) { str_split(str_replace(x, "## - ", ""), "\t")[[1]]})
      cnames <- ret[[1]]
      ret[[1]] <- NULL
      rnames <- lapply(ret, function(x) {x[1]})
      rnames <- do.call(c, rnames)
      data <- lapply(ret, function(x) {x[-1]})
      mat <- do.call(rbind, data)
      colnames(mat) <- cnames[-1]
      rownames(mat) <- rnames
      suppressWarnings(class(mat) <- "numeric")
      mat
  }
  ret <- lapply(scoreTestFiles, function(x) { model.to.matrix(read.null.model(x))})
  ret
};

