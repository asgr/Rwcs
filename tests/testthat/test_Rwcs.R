context("Check Rwcs world to pixel (s2p) and pixel to world (p2s)")
library(Rwcs)
library(testthat)

testmat=cbind(RA=10:20, Dec=20:30)

expect_equal(testmat, Rwcs_p2s(Rwcs_s2p(testmat)))
