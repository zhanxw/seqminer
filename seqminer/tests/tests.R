## follow http://adv-r.had.co.nz/Testing.html
library(testthat)
library(seqminer)
## test code are under inst/tests
## test_package("seqminer", reporter="tap")
test_package("seqminer")
