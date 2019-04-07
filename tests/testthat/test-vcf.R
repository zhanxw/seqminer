context("Test read from VCF file")
outFile <- system.file("tests/correct.Rdata", package = "seqminer")

sysname <- Sys.info()[['sysname']]
systemTestable <- FALSE
if (! sysname %in%  c("Linux", "Windows", "Darwin")) {
  cat("Skip unit-testing tabix, probably your system ", sysname, " is big endian!\n")
} else {
  systemTestable <- TRUE
}

if (file.exists(outFile) && systemTestable) {
  load(outFile, verbose = TRUE)

  test_that("readVCFToMatrixByRange", {
    cat("--------------- test readVCFToMatrixByRange ---------------\n")
    fileName <- system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")

    test.cfh <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Nonsynonymous")
    expect_equal(test.cfh, cfh.range)

    test.cfh.nonsyn <- readVCFToMatrixByRange(fileName, "1:196621007-196642234,1:196642235-196716634", "Nonsynonymous")
    expect_equal(test.cfh.nonsyn, cfh.nonsyn)

    test.cfh.nonsyn.2 <- readVCFToMatrixByRange(fileName, c("1:196621007-196642234", "1:196642235-196716634"), "Nonsynonymous")
    expect_equal(test.cfh.nonsyn, cfh.nonsyn)

    test.cfh.syn <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Synonymous")
    expect_equal(test.cfh.syn, cfh.syn)
  })

  test_that("test readVCFToMatrixByGene", {
    ## Extract genotypes
    fileName <- system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    geneFile <- system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")

    test.cfh.syn.2 <- readVCFToMatrixByGene(fileName, geneFile, "CFH", "Synonymous")
    expect_equal(test.cfh.syn.2, cfh.syn.2)

    test.apoe <- readVCFToMatrixByGene(fileName, geneFile, "APOE", "")
    expect_equal(test.apoe, apoe)

    expect_warning(ssss <- readVCFToMatrixByGene(fileName, geneFile, "ssss", ""))
  })

  test_that("test VCFToListByGene", {
    fileName <- system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    geneFile <- system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
    test.out.gene.1 <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM", "POS", "ID"), "", "")
    expect_equal(test.out.gene.1, out.gene.1)

    test.out.gene.2 <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM","ID"), c("AC","AN"), "")
    expect_equal(test.out.gene.2, out.gene.2)

    test.out.gene.3 <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM","ID"), c("AC","AN"), c("GT","GQ"))
    expect_equal(test.out.gene.3, out.gene.3)
  })

  test_that("test VCFToListByRange", {
    fileName <- system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    test.out.range.1 <- readVCFToListByRange(fileName, "1:196621007-196716634", "Synonymous", c("CHROM","ID", "POS"), c("AC","AN"), c("GT","GQ"))
    expect_equal(test.out.range.1, out.range.1)
  })

  ## test bcf
  test_that("readVCFToMatrixByRange.BCF", {
    ## cat("--------------- test readVCFToMatrixByRange ---------------\n")
    fileName <- system.file("vcf/all.anno.filtered.extract.bcf.gz", package = "seqminer")
    test.cfh <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Nonsynonymous")
    expect_equal(test.cfh, cfh)
  })

  ## test tabix functions
  test_that("tabix", {
    fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    snp <- tabix.read(fileName, "1:196623337-196632470")
    expect_equal(class(snp), "character")
    expect_equal(length(snp), 3)
  })

  test_that("tabix", {
    fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    snp <- tabix.read(fileName, "1:196633607-196633607")
    expect_equal(class(snp), "character")
    expect_equal(length(snp), 0)
  })
  test_that("tabix.read.table", {
    fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    snp <- tabix.read.table(fileName, "1:196623337-196632470")
    expect_equal(class(snp), "data.frame")
    expect_equal(nrow(snp), 3)
    expect_equal(ncol(snp), 12)
  })
  test_that("tabix.read.table", {
    fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    snp <- tabix.read.table(fileName, "1:196633607-196633607")
    expect_equal(class(snp), "data.frame")
    expect_equal(nrow(snp), 0)
    expect_equal(ncol(snp), 0)
  })
  test_that("readBGENToMatrixByRange", {
    fileName = system.file("bgen/all.anno.filtered.extract.bgen", package = "seqminer")
    cfh <- readBGENToMatrixByRange(fileName, "1:196621007-196716634")
    expect_length(cfh , 1)
    expect_equivalent(cfh , cfh.all)   ## compare value not dimnames
  })
  test_that("readBGENToMatrixByGene", {
    fileName = system.file("bgen/all.anno.filtered.extract.bgen", package = "seqminer")
    geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
    cfh <- readBGENToMatrixByGene(fileName, geneFile, "CFH")
    expect_length(cfh , 1)
    expect_equivalent(cfh , cfh.all)   ## compare value not dimnames
  })

  test_that("readBGENToListByRange", {
    fileName = system.file("bgen/all.anno.filtered.extract.bgen", package = "seqminer")
    cfh <- readBGENToListByRange(fileName, "1:196621007-196716634")
    expect_length(cfh, 8)
  })

  test_that("readBGENToListByGene", {
    fileName = system.file("bgen/all.anno.filtered.extract.bgen", package = "seqminer")
    geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
    cfh <- readBGENToListByGene(fileName, geneFile, "CFH")
    expect_length(cfh, 8)
  })

  test_that("readPlinkToMatrixByIndex", {
    fileName = system.file("plink/all.anno.filtered.extract.bed", package = "seqminer")
    fileName = sub(fileName, pattern = ".bed", replacement = "")
    sampleIndex = seq(3)
    markerIndex =c(14, 36) ## position 1:196642233 and 1:196659237
    cfh <- readPlinkToMatrixByIndex(fileName, sampleIndex, markerIndex)
    expect_equal(t(cfh), cfh.range[[1]])
  })

  test_that("openPlink", {
    cat("--------------- test openPlink ---------------\n")
    fileName = system.file("plink/all.anno.filtered.extract.bed", package = "seqminer")
    fileName = sub(fileName, pattern = ".bed", replacement = "")
    plinkObj <- openPlink(fileName)
    expect_length(plinkObj, 3)
  })

  test_that("readPlinkBySubscript", {
    cat("--------------- test readPlinkBySubscript ---------------\n")
    fileName = system.file("plink/all.anno.filtered.extract.bed", package = "seqminer")
    filePrefix = sub(fileName, pattern = ".bed", replacement = "")
    plinkObj = openPlink(filePrefix)
    sampleIndex = seq(3)
    markerIndex =c(14, 36)
    cfh <- plinkObj[sampleIndex, markerIndex]
    expect_equal(t(cfh), cfh.range[[1]])
  })

  test_that("readSingleChromosomeVCFToMatrixByRange", {
    cat("--------------- test readSingleChromosomeVCFToMatrixByRange ---------------\n")
    fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    cfh <- readSingleChromosomeVCFToMatrixByRange(fileName, "1:196621007-196716634")
    expect_equivalent(cfh[[1]], t(cfh.all[[1]]))
  })

  test_that("createSingleChromosomeVCFIndex", {
    cat("--------------- test createSingleChromosomeVCFIndex ---------------\n")
    fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
    cfh <- createSingleChromosomeVCFIndex(fileName, indexFileName = tempfile())
    expect(nchar(cfh) > 0)
  })

  test_that("readSingleChromosomeBCFToMatrixByRange", {
    cat("--------------- test readSingleChromosomeBCFToMatrixByRange ---------------\n")
    fileName = system.file("vcf/all.anno.filtered.extract.headerFixed.bcf.gz", package = "seqminer")
    cfh <- readSingleChromosomeBCFToMatrixByRange(fileName, "1:196621007-196716634")
    expect_equivalent(cfh[[1]], t(cfh.all[[1]]))
  })

  test_that("createSingleChromosomeBCFIndex", {
    cat("--------------- test createSingleChromosomeBCFIndex ---------------\n")
    fileName = system.file("vcf/all.anno.filtered.extract.headerFixed.bcf.gz", package = "seqminer")
    cfh <- createSingleChromosomeBCFIndex(fileName, indexFileName = tempfile())
    expect(nchar(cfh) > 0)
  })
}

if (FALSE) {
  ##  code that generate the correct output in a .Rdata file
  rm(list=ls())

  setwd("/net/fantasia/home/zhanxw/mycode/seqminer/seqminer/inst/tests")
  outFile <- "correct.Rdata"

  suppressPackageStartupMessages(library(seqminer))
  fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
  geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")

  cat("--------------- test readVCFToMatrixByRange ---------------\n")
  try(cfh <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Nonsynonymous"))
  try(cfh.nonsyn <- readVCFToMatrixByRange(fileName, "1:196621007-196642234,1:196642235-196716634", "Nonsynonymous"))
  try(cfh.nonsyn.2 <- readVCFToMatrixByRange(fileName, c("1:196621007-196642234", "1:196642235-196716634"), "Nonsynonymous"))
  try(cfh.syn <- readVCFToMatrixByRange(fileName, "1:196621007-196716634", "Synonymous"))

  cat("--------------- test readVCFToMatrixByGene ---------------\n")
  ## Extract genotypes
  try(cfh.syn.2 <- readVCFToMatrixByGene(fileName, geneFile, "CFH", "Synonymous"))
  try(apoe <- readVCFToMatrixByGene(fileName, geneFile, "APOE", ""))
  try(ssss <- readVCFToMatrixByGene(fileName, geneFile, "ssss", ""))

  ## Another way to extract from VCF File
  cat("--------------- test VCFToListByGene ---------------\n")
  try(out.gene.1 <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM", "POS", "ID"), "", ""))
  try(out.gene.2 <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM","ID"), c("AC","AN"), ""))
  try(out.gene.3 <- readVCFToListByGene(fileName, geneFile, "CFH", "Synonymous", c("CHROM","ID"), c("AC","AN"), c("GT","GQ")))
  ##print (t)

  cat("--------------- test VCFToListByRange ---------------\n")
  try(out.range.1 <- readVCFToListByRange(fileName, "1:196621007-196716634", "Synonymous", c("CHROM","ID", "POS"), c("AC","AN"), c("GT","GQ")))
  ## print(t)

  save.image(file = outFile)
  ## load(outFile)
}
