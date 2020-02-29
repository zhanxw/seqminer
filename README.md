SEQMINER
========

[![Build Status](https://travis-ci.org/zhanxw/seqminer.svg?branch=master)](https://travis-ci.org/zhanxw/seqminer)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/zhanxw/seqminer?branch=master&svg=true)](https://ci.appveyor.com/project/zhanxw/seqminer)
![](https://cranlogs.r-pkg.org/badges/seqminer)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/seqminer)](https://cran.r-project.org/package=seqminer)


**Table of Contents**

- [Introduction](#introduction)
- [Download](#download)
- [Showcase](#showcase)
  - [Index VCF/BCF files](#Index-VCF/BCF-files)
  - [Query VCF/BCF files](#Query-VCF/BCF-files)
  - [Query BGEN/PLINK files](#query-BGEN/PLINK-files)
  - [Command line linterface](#command-line-linterface)

# Introduction
Seqminer is a highly efficient R-package for retrieving sequence variants from biobank scale datasets of millions of individuals and billions of genetic variants. It supports all variant types, including multi-allelic variants and imputation dosages. It takes VCF/BCF/BGEN/PLINK format as input file, indexes, queries them based upon variant-based index and loads them as R data types such as list or matrix.

# Download
Install the development version ([devtools](https://github.com/hadley/devtools) package is required):

    devtools::install_github("zhanxw/seqminer")

# Showcase
Here are some examples of how to use seqminer to index and query files in real-life scenarios.

## Index VCF/BCF files

    library(seqminer)
    bcf.ref.file <- "input.bcf"
    bcf.idx.file <- "input.bcf.scIdx"
    out <- seqminer::createSingleChromosomeBCFIndex(bcf.ref.file, bcf.idx.file)

or

    vcf.ref.file <- "input.vcf.gz"
    vcf.idx.file <- "input.vcf.gz.scIdx"
    out <- seqminer::createSingleChromosomeVCFIndex(vcf.ref.file, vcf.idx.file)

This would generate variant-based index that works with commonly used sequence variant file format, such as VCF/BCF files.

## Query VCF/BCF files

Query VCF file:

    vcf.ref.file <-  "input.vcf.gz"
    vcf.idx.file <-  "input.vcf.gz.scIdx"
    tabix.range <- "1:123-1234"
    geno <- seqminer::readSingleChromosomeVCFToMatrixByRange(vcf.ref.file, tabix.range, vcf.idx.file)

Query BCF file:

    bcf.ref.file <- "input.bcf"
    bcf.idx.file <- "input.bcf.scIdx"
    tabix.range <- "1:123-1234"
    geno <- seqminer::readSingleChromosomeBCFToMatrixByRange(bcf.ref.file, tabix.range, bcf.idx.file)

Querying multiple regions is also doable, simply specify multiple regions and separte them by a comma, e.g. `"1:123-124,1:1234-1235"`.


## Query BGEN/PLINK files

Query BGEN file:

    bg.ref.file <- "input.bgen"
    bg.range <- "1:123-1234"
    geno.mat <- seqminer::readBGENToMatrixByRange(bg.ref.file, bg.range)
    geno.list <- seqminer::readBGENToListByRange(bg.ref.file, bg.range)

Query PLINK file:

    plink.ref.file <- "input"
    geno <- seqminer::readPlinkToMatrixByIndex(plink.ref.file, sampleIndex=1:20000, markerIndex=1:100)

## Command line linterface
We also developed a seqminer command line interface:

    ./queryVCFIndex.intel input.vcf.gz input.vcf.gz.scIdx 1:123-1234

Citation: 

[Zhan, X. and Liu, D. J. (2015), SEQMINER: An R-Package to Facilitate the Functional Interpretation of Sequence-Based Associations. Genet. Epidemiol., 39: 619–623. doi:10.1002/gepi.21918](http://onlinelibrary.wiley.com/doi/10.1002/gepi.21918/abstract)
