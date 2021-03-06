\name{Rwcs_setscreen}
\alias{Rwcs_setscreen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a WCS Projection
}
\description{
This is useful for creating a simple reference WCS to plot various images on top of, e.g. if you have many images that overlap on the sky and they have a complex WCS geometry- it is usually best to SWarp to a common simple WCS.
}
\usage{
Rwcs_setscreen(CRVAL1 = 0, CRVAL2 = 0, pixscale = 1, NAXIS1 = 1000, NAXIS2 = 1000,
  CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN", CD1_1 = pixscale/3600, CD1_2 = 0, CD2_1 = 0,
  CD2_2 = pixscale/3600, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CRVAL1}{
Numeric scalar; FITS header CRVAL1 for the CTYPE1 projection system. This is the RA in degrees at the location of \option{NAXIS1}/2. This should be chosen to be near the centre of the target projection (i.e. the average RA of the frames to be overlaid).
}
  \item{CRVAL2}{
Numeric scalar; FITS header CRVAL2 for the CTYPE2 projection system. This is the Dec in degrees at the location of \option{NAXIS2}/2. This should be chosen to be near the centre of the target projection (i.e. the average Dec of the frames to be overlaid).
}
  \item{pixscale}{
Numeric scalar; the target pixel scale in arc seconds. By default \option{CD1_1} and \option{CD2_2} inherit this in degress (i.e. divided by 3600).
}
  \item{NAXIS1}{
Numeric scalar; the number of pixels to display on the x-axis.
}
  \item{NAXIS2}{
Numeric scalar; the number of pixels to display on the y-axis.
}
  \item{CTYPE1}{
Character scalar; the RA projection system type. Must be a WCS legal type (up to 4 letter axis followed by 3 letter projection padded as required with "-", e.g.: "RA---TAN") and 8 characters long.
}
  \item{CTYPE2}{
Character scalar; the Dec projection system type. Must be a WCS legal type (up to 4 letter axis followed by 3 letter projection padded as required with "-", e.g.: "DEC--TAN") and 8 characters long. 
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
  \item{\dots}{
Extra arguments to be passed to \code{\link{Rwcs_image}}. These will be the decoration arguments in practice.
}
}
\details{
This is mostly useful for creating a simple projection of multiple over-lapping images, which can each be added to the current view using \code{Rwcs_image(image, add=TRUE)}. See Examples.
}
\value{
Run for the side effect of creating a new WCS view, updating the global options for current_keyvalues in the process.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{Rwcs_image}}
}
\examples{
#Recommended to use Rfits (will be supported longer term):

library(Rfits)
file_image = system.file('extdata', 'image.fits', package = "Rfits")
temp_image = Rfits_read_image(file_image)

Rwcs_setscreen(352.2914, -31.82235, 1, 200, 200)
Rwcs_image(temp_image, add=TRUE)

Rwcs_setscreen(352.2914, -31.82235, 0.5, 400, 400, CD1_2 = 1/3600, CD2_1 = -1/3600)
Rwcs_image(temp_image, add=TRUE)
}

