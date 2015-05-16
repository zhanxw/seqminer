\name{tabix.createIndex.meta}
\alias{tabix.createIndex.meta}
\title{Create tabix index for bgzipped MetaScore/MetaCov file}
\usage{
tabix.createIndex.meta(bgzipFile)
}
\arguments{
  \item{bgzipFile}{character, input vcf file}
}
\description{
Create tabix index for bgzipped MetaScore/MetaCov file
}
\examples{
fileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
tabix.createIndex.meta(fileName)
}
\seealso{
http://zhanxw.com/seqminer/ for online manual and examples

http://zhanxw.github.io/rvtests/ for rvtests
}

