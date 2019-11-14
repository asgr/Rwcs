\name{Rwcs_keypass}
\alias{Rwcs_keypass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WCS Header Creator
}
\description{
Ensures we make a WCS legal header.
}
\usage{
Rwcs_keypass(keyvalues = NULL, CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN", CRVAL1 = 0,
CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, PV1 = 0,
PV2 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keyvalues}{
List; header values to be used for the base WCS. If this is empty then a fresh header is created using the other arguments. Any specified argument (e.g. those below here) will take precedent other the same named list element in \option{keyvalues}. This means \code{Rwcs_keypass} is a simple and safe way to make small changes to a current legal FITS header object.
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
  \item{PV1}{
Numeric scalar; additional projection parameter 1.
}
  \item{PV2}{
Numeric scalar; additional projection parameter 2.
}
}
\details{
Simple function that safely creates a legal FITS \option{keyvalues} output, similar to the one made by \code{Rfits_read_image}.
}
\value{
\option{keyvalues} list output.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{Rwcs}}
}
\examples{
\dontrun{
library(Rfits)
file_image=system.file('extdata', 'image.fits', package = "Rfits")
temp_image=Rfits_read_image(file_image)
Rwcs_keypass(temp_image$keyvalues, CTYPE1 = 'RA---SIN', CTYPE2 = 'DEC--SIN')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
