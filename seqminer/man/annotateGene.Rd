\name{annotateGene}
\alias{annotateGene}
\title{Annotate a test variant}
\usage{
annotateGene(param, chrom, position, ref, alt)
}
\arguments{
  \item{param}{a list of annotation configuaration (e.g.
  reference file, gene definition)}

  \item{chrom}{a vector of chromosome names}

  \item{position}{a vector of chromosome positions}

  \item{ref}{a vector of reference alleles}

  \item{alt}{a vector of alternative alleles}
}
\value{
annotated results in a data frame structure
}
\description{
Annotate a test variant
}
\examples{
reference = "/net/fantasia/home/zhanxw/anno/resources/hs37d5.fa"
refBase <- getRefBase(reference, chrom = c(19, 19), pos = c(6718146, 6718387), len = c(1, 1))

param <- list(reference = "/net/fantasia/home/zhanxw/anno/resources/hs37d5.fa",
                  geneFile = "/net/fantasia/home/zhanxw/anno/resources/refFlat_hg19.txt.gz")
param <- makeAnnotationParameter(param)
print(param)
annotateGene(param, c("19", "19"), c(6718146, 6718387) , c("T", "G"), c("G", "C"))
}
\seealso{
makeAnnotationParameter
}

