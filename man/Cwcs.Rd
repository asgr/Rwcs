\name{Cwcs}
\alias{Cwcs}
\alias{Cwcs_s2p}
\alias{Cwcs_p2s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mid Level WCSlib Interface
}
\description{
Offers basic mid level access to the C \code{wcslib} library for projecting and deprojecting [x,y] and [RA,Dec] image data using a target WCS. Faster than Rwcs functions, but less safety and utility. Advanced users can use these functions.
}
\usage{
Cwcs_s2p(RA, Dec, CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN", CRVAL1 = 0, CRVAL2 = 0,
  CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, PV1_1 = NA,
  PV1_2 = NA, PV2_1 = NA, PV2_2 = NA)
Cwcs_p2s(x, y, CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN", CRVAL1 = 0, CRVAL2 = 0,
  CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, PV1_1 = NA,
  PV1_2 = NA, PV2_1 = NA, PV2_2 = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RA}{
Numeric vector; target right ascension in degrees.
}
  \item{Dec}{
Numeric vector; target declination in degrees.
}
  \item{x}{
Numeric vector; target x-pixel.
}
  \item{y}{
Numeric vector; target y-pixel.
}
  \item{CTYPE1}{
Scalar character; the RA projection system type. Must be a WCS legal type and 8 characters long.
}
  \item{CTYPE2}{
Scalar character; the Dec projection system type. Must be a WCS legal type and 8 characters long.
}
  \item{CRVAL1}{
Numeric scalar; FITS header CRVAL1 for the CTYPE1 projection system. This is the RA in degrees at the location of CRPIX1.
}
  \item{CRVAL2}{
Numeric scalar; FITS header CRVAL2 for the CTYPE2 projection system. This is the Dec in degrees at the location of CRPIX2.
}
  \item{CRPIX1}{
Numeric scalar; FITS header CRPIX1 for the CTYPE1 projection system. This is the x pixel value at the location of CRVAL1.
}
  \item{CRPIX2}{
Numeric scalar; FITS header CRPIX2 for the CTYPE2 projection system. This is the y pixel value at the location of CRVAL2.
}
  \item{CD1_1}{
Numeric scalar; FITS header CD1_1 for the CTYPE1 projection system. Change in CTYPE1 in degrees along x-Axis.
}
  \item{CD1_2}{
Numeric scalar; FITS header CD1_2 for the CTYPE1 projection system. Change in CTYPE1 in degrees along y-Axis.
}
  \item{CD2_1}{
Numeric scalar; FITS header CD2_1 for the CTYPE2 projection system. Change in CTYPE2 in degrees along x-Axis.
}
  \item{CD2_2}{
Numeric scalar; FITS header CD2_2 for the CTYPE2 projection system. Change in CTYPE2 in degrees along y-Axis.
}
  \item{PV1_1}{
Numeric scalar; additional axis 1 projection parameter 1.
}
  \item{PV1_2}{
Numeric scalar; additional axis 1 projection parameter 2.
}
  \item{PV2_1}{
Numeric scalar; additional axis 1 projection parameter 1.
}
  \item{PV2_2}{
Numeric scalar; additional axis 1 projection parameter 2.
}
}
\details{
The function input and output is deliberately simple.
}
\value{
\code{Cwcs_s2p} is the spherical to pixel (hence s2p) projection function, i.e. RA/Dec to x/y. The output is a two column numeric matrix containing [x,y] data.

\code{Cwcs_p2s} is the pixel to spherical (hence p2s) projection function, i.e. x/y to RA/Dec. The output is a two column numeric matrix containing [RA,Dec] data.
}
\references{
Greisen, E.W. & Calabretta, M.R. (2002), Astronomy & Astrophysics, 395, 1061-1075
Calabretta, M.R., & Greisen, E.W., (2002), Astronomy & Astrophysics, 395, 1077-1122
}
\author{
Aaron Robotham
Rodrgio Tobar
}
\seealso{
\code{\link{Cwcs}}
}
\examples{
Cwcs_s2p(10:20, 20:30)
Cwcs_p2s(10:20, 20:30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
