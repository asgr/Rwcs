\name{Rwcs}
\alias{Rwcs}
\alias{Rwcs_s2p}
\alias{Rwcs_p2s}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High Level WCSlib Interface
}
\description{
Offers basic high level access to the C \code{wcslib} library for projecting and deprojecting [x,y] and [RA,Dec] image data using a target WCS. Slower than Cwcs functions, but safer and with more utility. Normal users should use these functions.
}
\usage{
Rwcs_s2p(RA, Dec, keyvalues = NULL, pixcen = "FITS", loc.diff = c(0, 0),
  coord.type = "deg", sep = ":", ...)
Rwcs_p2s(x, y, keyvalues = NULL, pixcen = "FITS", loc.diff = c(0, 0),
  coord.type = "deg", sep = ":", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RA}{
Numeric vector/matrix; target right ascension in degrees. If matrix then the first column will be used as RA and the second column as Dec.
}
  \item{Dec}{
Numeric vector; target declination in degrees. Ignored if \option{RA} is a matrix.
}
  \item{x}{
Numeric vector/matrix; target x-pixel. If Matrix then the first column will be used as the x-axis and the second column as y-axis.
}
  \item{y}{
Numeric vector; target y-pixel. Ignored if \option{x} is a matrix.
}
  \item{keyvalues}{
List; header values to be used for the WCS. The \option{\dots} can be used to over-ride any values for convenience.
}
  \item{pixcen}{
Character scalar; either 'FITS' (pixel centres are integer, e.g. the centre of the bottom-left pixel is [1,1]) of 'R' (pixel centres are half integer, e.g. the centre of the bottom-left pixel is [0.5,0.5]). The default is 'FITS', to be consistent with other programs, but within \code{R} images are usually interpretted with a half pixel offset, so for converting object image [x,y] positions to true WCS positions option 'R' should be used.
}
  \item{loc.diff}{
Numeric vector; additional [x,y] pixel offset to apply (useful if using a cutout from a larger image when the WCS has not been fixed).
}
  \item{coord.type}{
Character scalar; the format of the input/output RA and Dec. Allowed options are 'deg' for degress (default) and 'sex' for sexigesimal (i.e. HMS for RA and DMS for Deg).
}
  \item{sep}{
Character scalar; when coord.type='sex', \option{sep} defines the type of separator used for the HMS and DMS strings (i.e. H:M:S and D:M:S would be \option{sep}=':', which is the default). See \code{\link{hms2deg}} and \code{\link{dms2deg}} for more details.
}
  \item{\dots}{
Optional keyvalues to overwrite the \option{keyvalues} entries with. If not present, the defaults are CRVAL1=0, CRVAL2=0, CRPIX1=0, CRPIX2=0, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1, PV1=0, PV2=0, CTYPE1='RA---TAN', CTYPE2='DEC--TAN'.
}
}
\details{
The function input and output is deliberately simple. The important thing is to pass in a senesible WCS to the \option{keyvalues} argument. This would be the same \option{keyvalues} output supplied by the \code{Rfits} package when reading in an image.

The following projections are supported:

AZP: zenithal/azimuthal perspective
SZP: slant zenithal perspective
TAN: gnomonic
STG: stereographic
SIN: orthographic/synthesis
NCP: unofficially supported SIN-like projection
ARC: zenithal/azimuthal equidistant
ZPN: zenithal/azimuthal polynomial
ZEA: zenithal/azimuthal equal area
AIR: Airy’s projection
CYP: cylindrical perspective
CEA: cylindrical equal area
CAR: plate carrée
MER: Mercator’s projection
COP: conic perspective
COE: conic equal area
COD: conic equidistant
COO: conic orthomorphic
SFL: Sanson-Flamsteed (“global sinusoid”)
PAR: parabolic
MOL: Mollweide’s projection
AIT: Hammer-Aitoff
BON: Bonne’s projection
PCO: polyconic
TSC: tangential spherical cube
CSC: COBE quadrilateralized spherical cube
QSC: quadrilateralized spherical cube
HPX: HEALPix
XPH: HEALPix polar, aka “butterfly”
}
\value{
\code{Rwcs_s2p} is the spherical to pixel (hence s2p) projection function, i.e. RA/Dec to x/y. The output is a two column numeric matrix containing [x,y] data.

\code{Rwcs_p2s} is the pixel to spherical (hence p2s) projection function, i.e. x/y to RA/Dec. The output is a two column numeric matrix containing [RA,Dec] data.
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
\code{\link{Rwcs}}
}
\examples{
Rwcs_s2p(10:20, 20:30)
Rwcs_p2s(Rwcs_s2p(10:20, 20:30))
Rwcs_p2s(Rwcs_s2p(10:20, 20:30, CTYPE1='RA---SIN', CTYPE2='DEC--SIN'), CTYPE1='RA---SIN', 
CTYPE2='DEC--SIN')

#A more complicated example, where we transform and rotate large amounts:

exdata_start=expand.grid(1:10,21:30)
plot(exdata_start)
exradec=Rwcs_s2p(exdata_start, CRVAL1=20, CRPIX1=100, CRVAL2=30, CRPIX2=130, CD1_1=0.1,
CD1_2=-0.05, CD2_1=0.05, CD2_2=0.1)
plot(exradec)
exdata_end=Rwcs_p2s(exradec, CRVAL1=20, CRPIX1=100, CRVAL2=30, CRPIX2=130, CD1_1=0.1,
CD1_2=-0.05, CD2_1=0.05, CD2_2=0.1)
plot(exdata_start,cex=2)
points(exdata_end,col='red')

#The residuals should be very small (in the noise of double precision arithmetic):

plot(density(exdata_start[,1]-exdata_end[,1]))
lines(density(exdata_start[,2]-exdata_end[,2]),col='red')

#Different world projections (ones not listed here do not work well across large angles):

library(magicaxis)

world.grid=expand.grid(seq(0,360,by=10), seq(-90,90,by=10))

magplot(Rwcs_s2p(world.grid, CTYPE1='RA---ARC', CTYPE2='DEC--ARC'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---ZEA', CTYPE2='DEC--ZEA'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CYP', CTYPE2='DEC--CYP'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CEA', CTYPE2='DEC--CEA'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CAR', CTYPE2='DEC--CAR'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---MER', CTYPE2='DEC--MER'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---SFL', CTYPE2='DEC--SFL'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---PAR', CTYPE2='DEC--PAR'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---PCO', CTYPE2='DEC--PCO'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---MOL', CTYPE2='DEC--MOL'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---AIT', CTYPE2='DEC--AIT'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---TSC', CTYPE2='DEC--TSC'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CSC', CTYPE2='DEC--CSC'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---QSC', CTYPE2='DEC--QSC'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---HPX', CTYPE2='DEC--HPX'), pch='.', cex=3)
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---XPH', CTYPE2='DEC--XPH'), pch='.', cex=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
