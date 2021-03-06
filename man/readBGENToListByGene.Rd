% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqminer.R
\name{readBGENToListByGene}
\alias{readBGENToListByGene}
\title{Read information from BGEN file in a given range and return a list}
\usage{
readBGENToListByGene(fileName, geneFile, geneName)
}
\arguments{
\item{fileName}{character, represents an input BGEN file (Bgzipped, with Tabix index)}

\item{geneFile}{character, a text file listing all genes in refFlat format}

\item{geneName}{character vector, which gene(s) to be extracted}
}
\value{
a list of chrom, pos, varid, rsid, alleles, isPhased, probability, sampleId
}
\description{
Read information from BGEN file in a given range and return a list
}
\examples{
fileName = system.file("bgen/all.anno.filtered.extract.bgen", package = "seqminer")
geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
cfh <- readBGENToListByGene(fileName, geneFile, "CFH")
}
\seealso{
http://zhanxw.com/seqminer/ for online manual and examples
}
