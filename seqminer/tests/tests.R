## follow http://adv-r.had.co.nz/Testing.html
library(testthat)
library(seqminer)
test_package("seqminer", reporter="tap")
