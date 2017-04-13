library(testthat)
library(seqminer)
context("Test bi-allelic and tri-allelic")
## setwd("~/mycode/seqminer/debug/triallelic")

getFile <- function(x) {
  fn <- system.file(sprintf("test-triallelic/%s" , x), package = "seqminer")
  if (!file.exists(fn)) {
    warning(message("File ", fn, " does not exist!"))
  }
  fn
}

test_that("Read bi-allelic", {
  res <- res.bi <- rvmeta.readDataByRange(
           getFile("bi.MetaScore.assoc.gz"),
           getFile("bi.MetaCov.assoc.gz"),
           "1:1-3")
  ##print(str(res))
  print(res[[1]]$pos)
  print(res[[1]]$alt)
  print(res[[1]]$cov)
  expect_equal(res[[1]]$pos, sprintf("1:%d", seq(3)))
  expect_equal(res[[1]]$alt[[1]], rep("G", 3))
  expect_equal(dim(res[[1]]$cov[[1]]),  c(3, 3))
  expect_equal(diag(res[[1]]$cov[[1]]),  c(0.5126420, 0.1993610, 0.1139210))
  expect_equal(res[[1]]$cov[[1]], t(res[[1]]$cov[[1]]))
})

test_that("Read bi-allelic", {
  res <- res.tri.as.bi <- (rvmeta.readDataByRange(
    getFile('tri.MetaScore.assoc.gz'),
    getFile("tri.MetaCov.assoc.gz"),
    "1:1-3"))
  ##print(str(res))
  print(res[[1]]$pos)
  print(res[[1]]$alt)
  print(res[[1]]$cov)
  expect_equal(res[[1]]$pos, sprintf("1:%d", c(1, 3)))
  expect_equal(res[[1]]$alt[[1]], c("T", "G"))
  expect_equal(dim(res[[1]]$cov[[1]]),  c(2, 2))
  expect_equal(diag(res[[1]]$cov[[1]]),  c(0.1993610, 0.1139210))
  expect_equal(res[[1]]$cov[[1]], t(res[[1]]$cov[[1]]))
})

test_that("Read multi-allelic", {
  res <- res.tri.as.tri <- rvmeta.readDataByRange(
           getFile('tri.MetaScore.assoc.gz'),
           getFile("tri.MetaCov.assoc.gz"),
           "1:1-3", multiAllelic = TRUE)
  ##print(str(res))
  print(res[[1]]$pos)
  print(res[[1]]$alt)
  print(res[[1]]$cov)
  expect_equal(res[[1]]$pos, c("1:1_A/G", "1:1_A/T", "1:3_A/G"))
  expect_equal(res[[1]]$alt[[1]], c("G", "T", "G"))
  expect_equal(dim(res[[1]]$cov[[1]]),  c(3, 3))
  expect_equal(diag(res[[1]]$cov[[1]]),  c(0.5126420, 0.1993610, 0.1139210))
  expect_equal(res[[1]]$cov[[1]], t(res[[1]]$cov[[1]]))
})


test_that("Read covariate", {
  print("=============")
  res <- rvmeta.readCovByRange(getFile("bi.MetaCov.assoc.gz"), "1:1-3")
  print(res)
  expect_equal(colnames(res), rownames(res))
  expect_equal(colnames(res), c("1:1", "1:2", "1:3"))
  expect_equal((res), t(res))
  expect_equivalent(diag(res), c(0.5126420, 0.1993610, 0.1139210)) ## do not compare labels
  print(res)

  res <- rvmeta.readCovByRange(getFile("tri.MetaCov.assoc.gz"), "1:1-3")
  expect_equal(colnames(res), rownames(res))
  expect_equal(colnames(res), c("1:1", "1:1/1", "1:3"))
  expect_equal((res), t(res))
  expect_equivalent(diag(res), c(0.5126420, 0.1993610, 0.1139210)) ## do not compare labels
  print(res)
})

test_that("Read score statistics", {
  print("=============")
  res <- rvmeta.readScoreByRange(getFile("bi.MetaScore.assoc.gz"), "1:1-3")
  print(res)
  expect_equal(res$pos, c(1, 2, 3))
  expect_equal(res$pVal, c(0.894958, 1.000000, 0.747845))

  res <- rvmeta.readScoreByRange(getFile("tri.MetaScore.assoc.gz"), "1:1-3")
  print(res)
  expect_equal(res$pos, c(1, 1, 3))
  expect_equal(res$pVal, c(0.894958, 1.000000, 0.747845))
})

test_that("Read from multiple-studies", {
  score.files <- c("test2.study1.MetaScore.assoc.gz", "test2.study2.MetaScore.assoc.gz")
  cov.files <- c("test2.study1.MetaCov.assoc.gz", "test2.study2.MetaCov.assoc.gz")
  score.files <- sapply(score.files, getFile)
  cov.files <- sapply(cov.files, getFile)
  res <- rvmeta.readDataByRange(score.files, cov.files, "1:1-3")
  expect_equal(res[[1]]$alt[[1]], c("T", "G"))
  expect_equal(res[[1]]$alt[[2]], c("C", "G"))
  expect_equal(res[[1]]$pVal[[1]], c(1.000000, 0.747845))
  expect_equal(res[[1]]$pVal[[2]], c(0.894958, 0.747845))
  expect_equal(res[[1]]$cov[[1]], matrix(c(0.199361, -0.0284801,-0.0284801,0.113921), 2, 2))
  expect_equal(res[[1]]$cov[[2]], matrix(c(0.512642, 0.0427202, 0.0427202, 0.113921), 2, 2))

  res <- rvmeta.readDataByRange(score.files, cov.files, "1:1-3", multiAllelic = TRUE)
  expect_equal(res[[1]]$alt[[1]], c("C", "T", "G"))
  expect_equal(res[[1]]$alt[[2]], c("C", NA, "G"))
  expect_equal(res[[1]]$pVal[[1]], c(0.894958, 1.000000, 0.747845))
  expect_equal(res[[1]]$pVal[[2]], c(0.894958, NA, 0.747845))
  expect_equal(diag(res[[1]]$cov[[1]]), c(0.5126420, 0.199361,0.113921))
  expect_equal(diag(res[[1]]$cov[[2]]), c(0.512642, NA, 0.113921))
})

