\name{readVCFToMatrixByGene}
\alias{readVCFToMatrixByGene}
\title{Read a gene from VCF file and return a genotypes matrix}
\usage{
  readVCFToMatrixByGene(fileName, geneFile, geneName,
    annoType)
}
\arguments{
  \item{fileName}{charactr, represents an input VCF file
  (Bgzipped, with Tabix index)}

  \item{geneFile}{character, a text file listing all genes
  in refFlat format}

  \item{geneName}{character vector, which gene(s) to be
  extracted}

  \item{annoType}{character, annotated types you would like
  to extract, such as "Nonsynonymous", "Synonymous". This
  can be left empty.}
}
\value{
  genotype matrix
}
\description{
  Read a gene from VCF file and return a genotypes matrix
}
\examples{
fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
cfh <- readVCFToMatrixByGene(fileName, geneFile, "CFH", "Synonymous")
}
\seealso{
  http://zhanxw.com/seqminer/ for online manual and
  examples
}

