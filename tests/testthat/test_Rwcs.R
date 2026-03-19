context("Check Rwcs world to pixel (s2p) and pixel to world (p2s)")
library(Rwcs)
library(testthat)

testmat=cbind(RA=10:20, Dec=20:30)

expect_equal(testmat, Rwcs_p2s(Rwcs_s2p(testmat, inherit=FALSE), inherit=FALSE))

# Example WCS header for more comprehensive tests
header <- c("SIMPLE  =                    T / file does conform to FITS standard",
            "BITPIX  =                  -32 / number of bits per data pixel",
            "NAXIS   =                    2 / number of data axes",
            "NAXIS1  =                  356 / length of data axis 1",
            "NAXIS2  =                  356 / length of data axis 2",
            "EXTEND  =                    T / FITS dataset may contain extensions",
            "COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy",
            "COMMENT and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H",
            "EQUINOX =                 2000 / equinox of celestial coord. system",
            "EPOCH   =                 2000 / epoch of celestial coord. system",
            "WCSAXES =                    2 / Number of World Coordinate System axes",
            "CRPIX1  =                  178 / ref pixel x",
            "CRPIX2  =                  178 / ref pixel y",
            "CRVAL1  =          352.2914408 / ref pixel x value",
            "CRVAL2  =          -31.8223455 / ref pixel y value",
            "CTYPE1  = 'RA---TAN'           / the coordinate type for the first axis",
            "CTYPE2  = 'DEC--TAN'           / the coordinate type for the second axis",
            "CUNIT1  = 'deg     '           / Axis unit",
            "CUNIT2  = 'deg     '           / Axis unit",
            "CD1_1   = -0.00009416666295793 / pixel size x",
            "CD1_2   =                    0 / xy rotation",
            "CD2_1   =                    0 / yx rotation",
            "CD2_2   =  0.00009416666295793 / pixel size y",
            "OBJECT  = '352.2914_-31.8223'  / GAMA CATAID",
            "RADESYS = 'ICRS    '           / Astrometric system")

raw = paste(formatC(substr(header,1,79), width=80, flag='-'),sep='',collapse = '')

# Pixel to sky and sky to pixel round-trip
test_that("Pixel to sky and back round-trip works", {
  x <- c(100, 512.5, 900)
  y <- c(100, 512.5, 900)
  nkey <- length(header)
  sky <- Cwcs_head_p2s(x, y, raw, nkey)
  pix <- Cwcs_head_s2p(sky[,1], sky[,2], raw, nkey)
  expect_equal(x, pix[,1], tolerance = 1e-5)
  expect_equal(y, pix[,2], tolerance = 1e-5)
})

# Sky to pixel and pixel to sky round-trip
test_that("Sky to pixel and pixel to sky round-trip is accurate", {
  ra <- c(351.95, 352.0, 352.05)
  dec <- c(-31.95, -32.0, -32.05)
  nkey <- length(header)
  pix <- Cwcs_head_s2p(ra, dec, raw, nkey)
  sky <- Cwcs_head_p2s(pix[,1], pix[,2], raw, nkey)
  expect_equal(ra, sky[,1], tolerance = 1e-5)
  expect_equal(dec, sky[,2], tolerance = 1e-5)
})

# NA and out-of-bounds handling
test_that("NA and out-of-bounds handling", {
  nkey <- length(header)
  pix <- Cwcs_head_s2p(352, NA, raw, nkey)
  expect_true(is.na(pix[1,1]))
  expect_true(is.na(pix[1,2]))
})
