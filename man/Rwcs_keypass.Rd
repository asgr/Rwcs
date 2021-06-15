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
  CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1,
  RADESYS = 'ICRS', EQUINOX = 'infer', PV1_1 = NA, PV1_2 = NA, PV1_3 = NA, PV2_1 = NA,
  PV2_2 = NA, PV2_3 = NA, PV2_4 = NA, PV2_5 = NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keyvalues}{
List; header values to be used for the base WCS. If this is empty then a fresh header is created using the other arguments. Any specified argument (e.g. those below here) will take precedent other the same named list element in \option{keyvalues}. This means \code{Rwcs_keypass} is a simple and safe way to make small changes to a current legal FITS header object.
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
  \item{\dots}{
Does nothing, but means it will not fail if people pass in extra header info.  
}
}
\details{
Simple function that safely creates a legal FITS \option{keyvalues} output, similar to the one made by \code{Rfits_read_image}.

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
\option{keyvalues} list output.
}
\author{
Aaron Robotham
}
\note{
Note that SIP / TPV / DSS / WAT / TPD is only supported with the use of raw headers (see Examples).
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
