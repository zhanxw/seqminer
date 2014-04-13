\name{rvmeta.readNullModel}
\alias{rvmeta.readNullModel}
\title{Read null model statistics}
\usage{
rvmeta.readNullModel(scoreTestFiles)
}
\arguments{
  \item{scoreTestFiles}{character vector, score test output
  files (rvtests outputs using --meta score)}
}
\value{
a list of statistics fitted under the null mode (without
genetic effects)
}
\description{
Read null model statistics
}
\examples{
scoreFileName = system.file("rvtests/rvtest.MetaScore.assoc.anno.gz", package = "seqminer")
}
\seealso{
http://zhanxw.com/seqminer/ for online manual and examples
}

