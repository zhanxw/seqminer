\name{tabix.read.header}
\alias{tabix.read.header}
\title{Read tabix file, similar to running tabix in command line.}
\usage{
  tabix.read.header(tabixFile, skippedLine = FALSE)
}
\arguments{
  \item{tabixFile}{character, an tabix indexed file}

  \item{skippedLine}{logitcal, whether to read tabix
  skipped lines (when used 'tabix -S NUM')}
}
\value{
  a list
}
\description{
  Read tabix file, similar to running tabix in command
  line.
}
\examples{
fileName = system.file("vcf/all.anno.filtered.extract.vcf.gz", package = "seqminer")
snp <- tabix.read.header(fileName)
}
\seealso{
  http://zhanxw.com/seqminer/ for online manual and
  examples
}

