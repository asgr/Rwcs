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

# Headers for a cube and a 4D array. wcslib takes the number of coordinate axes
# from WCSAXES, or from NAXIS when WCSAXES is absent, so these parse to structs
# with three and four axes even though only RA and Dec are ever projected.
set_key <- function(lines, kw, val) {
  i <- grep(paste0("^", kw, "\\s*="), lines)
  if (!length(i)) stop("no such keyrecord: ", kw)
  lines[i] <- sprintf("%-8s= %s", kw, val)
  lines
}

as_raw <- function(lines) {
  raw <- paste(formatC(substr(lines, 1, 79), width = 80, flag = '-'),
               sep = '', collapse = '')
  list(raw = raw, nkey = nchar(raw) / 80)
}

nd_header <- function(ndim) {
  l <- header
  if (ndim >= 3) {
    l <- c(l,
      "NAXIS3  =                  100 / length of data axis 3",
      "CRPIX3  =                    1 / ref pixel z",
      "CRVAL3  =              1.4e-03 / ref pixel z value",
      "CTYPE3  = 'FREQ    '           / the coordinate type for the third axis",
      "CUNIT3  = 'Hz      '           / Axis unit",
      "CD1_3   =                    0 / xz rotation",
      "CD2_3   =                    0 / yz rotation",
      "CD3_1   =                    0 / zx rotation",
      "CD3_2   =                    0 / zy rotation",
      "CD3_3   =                    1 / pixel size z")
  }
  if (ndim >= 4) {
    l <- c(l,
      "NAXIS4  =                    8 / length of data axis 4",
      "CRPIX4  =                    1 / ref pixel w",
      "CRVAL4  =                    0 / ref pixel w value",
      "CTYPE4  = 'STOKES  '           / the coordinate type for the fourth axis",
      "CD4_4   =                    1 / pixel size w")
  }
  l <- set_key(l, "NAXIS", sprintf("%20d", ndim))
  # The bundled header declares WCSAXES = 2; raise it in step with NAXIS.
  set_key(l, "WCSAXES", sprintf("%20d", ndim))
}

test_that("a cube header projects RA/Dec identically to its 2D parent", {
  ref <- Cwcs_head_p2s(c(100, 512.5, 900), c(100, 512.5, 900), raw, length(header))
  for (ndim in 3:4) {
    h <- as_raw(nd_header(ndim))
    # Multiple points used to return NULL, which Rwcs_p2s turned into zeros.
    multi <- Cwcs_head_p2s(c(100, 512.5, 900), c(100, 512.5, 900), h$raw, h$nkey)
    expect_equal(dim(multi), c(3L, 2L))
    expect_equal(multi, ref, tolerance = 1e-10)
    # A single point used to skip wcslib's ncoord/nelem guard entirely and write
    # past the end of the image-coordinate buffer.
    single <- Cwcs_head_p2s(100, 100, h$raw, h$nkey)
    expect_equal(single[1, ], ref[1, ], tolerance = 1e-10)
  }
})

test_that("the sky-to-pixel direction round-trips through a cube header", {
  for (ndim in 2:4) {
    h <- if (ndim == 2) list(raw = raw, nkey = length(header)) else as_raw(nd_header(ndim))
    ra <- c(351.95, 352.0, 352.05)
    dec <- c(-31.95, -32.0, -32.05)
    pix <- Cwcs_head_s2p(ra, dec, h$raw, h$nkey)
    sky <- Cwcs_head_p2s(pix[, 1], pix[, 2], h$raw, h$nkey)
    expect_equal(sky[, 1], ra, tolerance = 1e-5)
    expect_equal(sky[, 2], dec, tolerance = 1e-5)
  }
})

test_that("repeated projection against a cube header stays stable", {
  # The out-of-bounds write scaled with the number of axes, so repeat enough
  # calls that corrupting adjacent heap would perturb the result.
  h <- as_raw(nd_header(3))
  x <- c(100, 512.5, 900)
  want <- Cwcs_head_p2s(x, x, h$raw, h$nkey)
  expect_equal(dim(want), c(3L, 2L))
  for (i in 1:500) {
    expect_identical(Cwcs_head_p2s(x, x, h$raw, h$nkey), want)
  }
})

test_that("an unusable WCS is reported rather than answered with zeros", {
  # A singular CD matrix leaves the two celestial axes undefined.
  sing <- set_key(header, "CD2_2", "                   0 / yx rotation")
  h <- as_raw(sing)
  expect_error(Cwcs_head_p2s(c(100, 512.5), c(100, 512.5), h$raw, h$nkey),
               "singular")
  expect_error(Cwcs_p2s(c(100, 512.5), c(100, 512.5), CTYPE1 = "RA---NOPE"),
               "wcsset")
  # An alternate label that is not present in the header.
  expect_error(Cwcs_head_p2s(1, 1, raw, length(header), WCSref = 20L), "alternate")
  # RA coupled to the spectral axis: RA/Dec are no longer separable from it.
  coupled <- as_raw(set_key(nd_header(3), "CD1_3", "            1e-07 / xz rotation"))
  expect_error(Cwcs_head_p2s(c(100, 512.5), c(100, 512.5),
                             coupled$raw, coupled$nkey),
               "coupled|separable")
})

test_that("parsing a header does not disturb later parses of it", {
  # wcspih() compacts the keyrecords it accepts in place, so it must be handed a
  # private copy.  Re-using and interleaving the same caller-held string is what
  # exposes any cross-call contamination.
  h2 <- as_raw(header)
  h3 <- as_raw(nd_header(3))
  x <- c(100, 512.5, 900)
  a2 <- Cwcs_head_p2s(x, x, h2$raw, h2$nkey)
  a3 <- Cwcs_head_p2s(x, x, h3$raw, h3$nkey)
  # B then A then A again must match A alone.
  invisible(Cwcs_head_p2s(x, x, h3$raw, h3$nkey))
  expect_identical(Cwcs_head_p2s(x, x, h2$raw, h2$nkey), a2)
  invisible(Cwcs_head_p2s(x, x, h2$raw, h2$nkey))
  expect_identical(Cwcs_head_p2s(x, x, h3$raw, h3$nkey), a3)
  # And the 2D and reduced 3D answers continue to agree.
  expect_equal(a3, a2, tolerance = 1e-10)
})

test_that("per-point failures still return NAs rather than erroring", {
  # A mix of valid and invalid sky positions: the good ones must still project.
  ra <- c(352.2914408, 172.29, 352.2914408)
  dec <- c(-31.8223455, 58.18, -31.8223455)
  out <- unname(Rwcs_s2p(cbind(ra, dec), header = raw, inherit = FALSE))
  expect_equal(out[1, 1], 178, tolerance = 1e-3)
  expect_true(is.na(out[2, 1]))
  expect_equal(out[3, 1], 178, tolerance = 1e-3)
})

