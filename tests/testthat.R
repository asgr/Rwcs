library(testthat)
library(Rwcs)

test_check("Rwcs")

testmat=cbind(RA=10:20, Dec=20:30)

expect_equal(testmat, Rwcs_p2s(Rwcs_s2p(testmat)))
