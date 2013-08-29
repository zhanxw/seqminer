\name{rvmeta.readDataByGene}
\alias{rvmeta.readDataByGene}
\title{Read association statistics by gene from METAL-format files. Both score statistics and covariance statistics will be extracted.}
\usage{
  rvmeta.readDataByGene(scoreTestFiles, covFiles, geneFile,
    geneName)
}
\arguments{
  \item{scoreTestFiles}{character vector, score test output
  files (rvtests outputs using --meta score)}

  \item{covFiles}{character vector, covaraite files
  (rvtests outputs using --meta cov)}

  \item{geneFile}{character, a text file listing all genes
  in refFlat format}

  \item{geneName}{character vector, which gene(s) to be
  extracted}
}
\value{
  a list of statistics including chromosome, position,
  allele frequency, score statistics, covariance and
  annotation(if input files are annotated).
}
\description{
  Read association statistics by gene from METAL-format
  files. Both score statistics and covariance statistics
  will be extracted.
}
\examples{
scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
cfh <- rvmeta.readDataByGene(scoreFileName, covFileName, geneFile, "CFH")
}
\seealso{
  http://zhanxw.com/seqminer/ for online manual and
  examples
}

