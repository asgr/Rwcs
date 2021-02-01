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
  CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, RADESYS = "ICRS",
  EQUINOX = 2000L, PV1_1 = NA, PV1_2 = NA, PV1_3 = NA, PV2_1 = NA, PV2_2 = NA, PV2_3 = NA,
  PV2_4 = NA, PV2_5 = NA)
Cwcs_p2s(x, y, CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN", CRVAL1 = 0, CRVAL2 = 0,
  CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, RADESYS = "ICRS",
  EQUINOX = 2000L, PV1_1 = NA, PV1_2 = NA, PV1_3 = NA, PV2_1 = NA, PV2_2 = NA, PV2_3 = NA,
  PV2_4 = NA, PV2_5 = NA)
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
Character scalar; the RA projection system type. Must be a WCS legal type (up to 4 letter axis followed by 3 letter projection padded as required with "-", e.g.: "RA---TAN") and 8 characters long. See Details.
}
  \item{CTYPE2}{
Character scalar; the Dec projection system type. Must be a WCS legal type (up to 4 letter axis followed by 3 letter projection padded as required with "-", e.g.: "DEC--TAN") and 8 characters long. See Details.
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
  \item{RADESYS}{
Character scalar; the RA-DEC system to use. See details.
}
  \item{EQUINOX}{
Character scalar; the equinox for the system. Usually 2000 for ICRS and FK5, and 1950 for all others. Historically used "EPOCH", but this is deprecated now (although this function will search for EPOCH if EQUINOX is missing and use that in place).
}
  \item{PV1_1}{
Numeric scalar; additional axis 1 projection parameter 1.
}
  \item{PV1_2}{
Numeric scalar; additional axis 1 projection parameter 2.
}
  \item{PV1_3}{
Numeric scalar; additional axis 1 projection parameter 3.
}
  \item{PV2_1}{
Numeric scalar; additional axis 2 projection parameter 1.
}
  \item{PV2_2}{
Numeric scalar; additional axis 2 projection parameter 2.
}
  \item{PV2_3}{
Numeric scalar; additional axis 2 projection parameter 3.
}
  \item{PV2_4}{
Numeric scalar; additional axis 2 projection parameter 4.
}
  \item{PV2_5}{
Numeric scalar; additional axis 2 projection parameter 5.
}
}
\details{
The function input and output is deliberately simple.

The following axis types are supported:

\tabular{ll}{
RA   \tab right ascension \cr 
DEC  \tab declination \cr 
GLON \tab galactic longitude \cr 
GLAT \tab galactic latitude \cr 
ELON \tab ecliptic longitude \cr 
ELAT \tab ecliptic latitude \cr 
HLON \tab helioecliptic longitude \cr 
HLAT \tab helioecliptic latitude \cr 
SLON \tab supergalactic longitude \cr 
SLAT \tab supergalactic latitude
}

The following projection types are supported:

\tabular{ll}{
AZP  \tab zenithal/azimuthal perspective \cr
SZP  \tab slant zenithal perspective \cr
TAN  \tab gnomonic \cr
STG  \tab stereographic \cr
SIN  \tab orthographic/synthesis \cr
NCP  \tab unofficially supported SIN-like projectio \cr
ARC  \tab zenithal/azimuthal equidistant \cr
ZPN  \tab zenithal/azimuthal polynomial \cr
ZEA  \tab zenithal/azimuthal equal area \cr
AIR  \tab Airy’s projection \cr
CYP  \tab cylindrical perspective \cr
CEA  \tab cylindrical equal area \cr
CAR  \tab plate carrée \cr
MER  \tab Mercator’s projection \cr
COP  \tab conic perspective \cr
COE  \tab conic equal area \cr
COD  \tab conic equidistant \cr
COO  \tab conic orthomorphic \cr
SFL  \tab Sanson-Flamsteed (“global sinusoid”) \cr
PAR  \tab parabolic \cr
MOL  \tab Mollweide’s projection \cr
AIT  \tab Hammer-Aitoff \cr
BON  \tab Bonne’s projection \cr
PCO  \tab polyconic \cr
TSC  \tab tangential spherical cube \cr
CSC  \tab COBE quadrilateralized spherical cube \cr
QSC  \tab quadrilateralized spherical cube \cr
HPX  \tab HEALPix \cr
XPH  \tab HEALPix polar, aka “butterfly”
}

The following RA-DEC systems are supported:

\tabular{ll}{
ICRS  \tab  International Celestial Reference System \cr
FK5   \tab  mean place, new (IAU 1984) system \cr
FK4   \tab  mean place, old (Bessell-Newcomb) system \cr
FK4-NO-E  \tab  mean place, old system but without e-terms \cr
GAPPT \tab  geocentric apparent place,  IAU 1984 system
}
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
Aaron Robotham \cr
Rodrgio Tobar \cr
Mark Calabretta (main \code{wcslib} author, not responsible for mistakes in the \code{Rwcs} package)
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
