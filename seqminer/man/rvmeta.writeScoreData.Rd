\name{rvmeta.writeScoreData}
\alias{rvmeta.writeScoreData}
\title{Write score-based association statistics files.}
\usage{
rvmeta.writeScoreData(rvmetaData, outName, createIndex = FALSE)
}
\arguments{
  \item{rvmetaData}{a list vector. It's usually read by
  rvmeta.readDataByRange or rvmeta.readDataByGene function}

  \item{outName}{character, a text indicating output file
  prefix}

  \item{createIndex}{boolean, (default FALSE), whether or
  not to create the index}
}
\value{
TRUE only if succeed
}
\description{
Write score-based association statistics files.
}
\examples{
scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
covFileName = system.file("rvtests/rvtest.MetaCov.assoc.gz", package = "seqminer")
geneFile = system.file("vcf/refFlat_hg19_6col.txt.gz", package = "seqminer")
cfh <- rvmeta.readDataByRange(scoreFileName, covFileName, "1:196621007-196716634")
rvmeta.writeScoreData(cfh, "cfh.MetaScore.assoc")
}
\seealso{
http://zhanxw.com/seqminer/ for online manual and examples
}

