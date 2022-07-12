\name{Rwcs}
\alias{Rwcs}
\alias{Rwcs_s2p}
\alias{Rwcs_p2s}
\alias{Rwcs_pixscale}
\alias{Rwcs_in_image}
\alias{Rwcs_overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High Level WCSlib Interface
}
\description{
Offers basic high level access to the C \code{wcslib} library for projecting and deprojecting [x,y] and [RA,Dec] image data using a target WCS. Slower than Cwcs functions, but safer and with more utility. Normal users should use these functions.
}
\usage{
Rwcs_s2p(RA, Dec, keyvalues = NULL, pixcen = "FITS", loc.diff = c(0, 0),
  coord.type = "deg", sep = ":", header = NULL, inherit = TRUE, WCSref = NULL, ctrl = 2L,
  ...)
  
Rwcs_p2s(x, y, keyvalues = NULL, pixcen = "FITS", loc.diff = c(0, 0),
  coord.type = "deg", sep = ":", header = NULL, inherit = TRUE, WCSref = NULL, ctrl = 2L,
  ...)
  
Rwcs_pixscale(keyvalues = NULL, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, type = 'old',
  dim = NULL)
  
Rwcs_in_image(RA, Dec, xlim, ylim, buffer = 0, plot = FALSE, style = 'points', ...)

Rwcs_overlap(keyvalues_test, keyvalues_ref, buffer = 0.5, plot = FALSE)
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
Character scalar; either 'FITS' (pixel centres are integer, e.g. the centre of the bottom-left pixel is [1,1]) or 'R' (pixel centres are half integer, e.g. the centre of the bottom-left pixel is [0.5,0.5]). The default is 'FITS', to be consistent with other programs, but within \code{R} images are usually interpretted with a half pixel offset, so for converting object image [x,y] positions to true WCS positions option 'R' should be used.
}
  \item{loc.diff}{
Numeric vector; additional [x,y] pixel offset to apply. For \code{Rwcs_p2s} the input \option{x} is seen internally as x + \option{loc.diff}[1]. For \code{Rwcs_s2p} the output \option{x} is given as  x - \option{loc.diff}[1]. This is useful if using a cutout from a larger image when the WCS has not been fixed.
}
  \item{coord.type}{
Character scalar; the format of the input/output RA and Dec. Allowed options are 'deg' for degress (default) and 'sex' for sexigesimal (i.e. HMS for RA and DMS for Deg).
}
  \item{sep}{
Character scalar; when coord.type='sex', \option{sep} defines the type of separator used for the HMS and DMS strings (i.e. H:M:S and D:M:S would be \option{sep}=':', which is the default). See \code{\link{hms2deg}} and \code{\link{dms2deg}} for more details.
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
  \item{type}{
Character scalar; either 'old' style pixel scale computed from CD matrix, or 'new' that detects the pixel scale in the centre of the image.
}
  \item{dim}{
Numeric vector; the image dimensions. This is only relevant for \option{type} = 'new', and if not provided it is taken from NAXIS1 / NAXIS2 from the \option{keyvalues}. The only tie this tends to be wrong (when you shoudl enter it manually) is when the image has been compressed and these key values have different meanings (check for ZNAXIS1 and ZNAXIS2 in that case).
}
  \item{header}{
To aid backwards compatibility, you can also pass in a 1D character vector of header values, e.g. the \option{hdr} component of \option{readFITS} from the \code{FITSio} package. \option{keyvalues} must be left as NULL in this case. If you want to use all distortion terms you can also pass in the raw header output of \code{Rfits_read_header_raw}.
}
  \item{inherit}{
Logical; manual over-ride for inheriting the current WCS. This is usually only an issue if no options are passed in (see first Examples).
}
  \item{WCSref}{
Scalar; either NULL for base WCS (default), or a letter [a-z] for the relevant projection desired. Only used if using a raw \option{header}, and obviously this has to exist within the provided raw \option{header}.  
}
  \item{ctrl}{
Integer scalar; if using raw headers, this will specify the level of verbosity printed by wcslib when passing it. 0 = none, 1 = summary, 2 = fails, 3 = fails + non-WCS, 4 = all.
}
  \item{xlim}{
Numeric vector; x extent of frame domain to test. If missing and \option{keyvalues} is provided then default will be c(0, NAXIS1) extent. Note domain limits should be provided in the 'R' sense, so the bottom left pixel spans x=[0,1] and y=[0,1] rather than x=[0.5,1.5] and y=[0.5,1.5] in FITS.
}
  \item{ylim}{
Numeric vector; y extent of frame domain to test. If missing and \option{keyvalues} is provided then default will be c(0, NAXIS2) extent.
}
  \item{buffer}{
Numeric scalar; an additional buffer to apply to the 'in domain' test. You might want to set this to a low value (like 1e-9) to avoid spurious overlap error due to inperfections in forward and inverse WCS projections.
}
  \item{plot}{
Logical; should a plot be shown showing the location of test points on the WCS of interest?  
}
  \item{style}{
Character scalar; way to show overlap, either red 'points' or a red filled 'polygon'.
}
  \item{keyvalues_test}{
List; a WCS header to test. The question of interest is whether any of \option{keyvalues_test} overlaps with \option{keyvalues_ref}.
}
  \item{keyvalues_ref}{
List; a WCS header to reference. The question of interest is whether any of \option{keyvalues_test} overlaps with \option{keyvalues_ref}.
}
  \item{\dots}{
For \code{Rwcs_p2s} and \code{Rwcs_s2p} optional keyvalues to overwrite the \option{keyvalues} entries with. If not present, the defaults are CTYPE1 = "RA---TAN", CTYPE2 = "DEC--TAN", CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, RADESYS = "ICRS", EQUINOX = 2000L, PV1_1 = NA, PV1_2 = NA, PV2_1 = NA, PV2_2 = NA.

For \code{Rwcs_in_image} arguments to pass into \code{Rwcs_s2p} to project your \option{RA} and \option{Dec} points of interest into your test frame domain.
}
}
\details{
The function input and output is deliberately simple. The important thing is to pass in a senesible WCS to the \option{keyvalues} argument. This would be the same \option{keyvalues} output supplied by the \code{Rfits} package when reading in an image.

For convenience the last used \option{keyvalues} is written to options()$keyvalues. \code{Rwcs_s2p} and \code{Rwcs_p2s} will look for this object if neither \option{keyvalues} or any \dots are passed in (which in both cases will prevent this happening). The Examples shows an example of this being used in practice, where the \code{lines} functions after each all sky projection correctly pick up the last used WCS.

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
\code{Rwcs_s2p} is the spherical to pixel (hence s2p) projection function, i.e. RA/Dec to x/y. The output is a two column numeric matrix containing [x,y] data.

\code{Rwcs_p2s} is the pixel to spherical (hence p2s) projection function, i.e. x/y to RA/Dec. The output is a two column numeric matrix containing [RA,Dec] data.

\code{Rwcs_pixscale} is the scalar image pixel scale in asec/pixel (so typically a value of 0.1-0.5 for modern survey instruments).

\code{Rwcs_in_image} returns a logical vector for the positions given by \option{RA} and \option{Dec}, answering the questions whether in projection they exist within the domain of \option{xlim} and \option{ylim}. This might be used to find out whether sources fall onto a frame for instance.

\code{Rwcs_overlap} returns a logical scalar, answering the question of whether any of \option{keyvalues_test} overlaps with \option{keyvalues_ref}.
}
\references{
Greisen, E.W. & Calabretta, M.R. (2002), Astronomy & Astrophysics, 395, 1061-1075
Calabretta, M.R., & Greisen, E.W., (2002), Astronomy & Astrophysics, 395, 1077-1122
}
\author{
Aaron Robotham \cr
Rodrgio Tobar
}
\note{
By R convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y. This differs to the FITS convention of pixel mid points being integer values. As such the R [x,y] = FITS [x-0.5,y-0.5]. This rarely matters too much in practice, but for accurate overlays you will want to get it right.

Note that SIP / TPV / DSS / WAT / TPD is only supported with the use of raw headers (see Examples).
}
\seealso{
\code{\link{Cwcs}}, \code{\link{Rwcs_keypass}}
}
\examples{
Rwcs_s2p(10:20, 20:30, inherit=FALSE)
Rwcs_p2s(Rwcs_s2p(10:20, 20:30, inherit=FALSE), inherit=FALSE)
Rwcs_p2s(Rwcs_s2p(10:20, 20:30, CTYPE1='RA---SIN', CTYPE2='DEC--SIN'), CTYPE1='RA---SIN',
  CTYPE2='DEC--SIN')

#A more complicated example, where we transform and rotate large amounts:

exdata_start=expand.grid(1:10,21:30)
plot(exdata_start)
exradec=Rwcs_s2p(exdata_start, CRVAL1=20, CRPIX1=100, CRVAL2=30, CRPIX2=130, CD1_1=0.1,
CD1_2=-0.05, CD2_1=0.05, CD2_2=0.1, inherit=FALSE)
plot(exradec)
exdata_end=Rwcs_p2s(exradec, CRVAL1=20, CRPIX1=100, CRVAL2=30, CRPIX2=130, CD1_1=0.1,
CD1_2=-0.05, CD2_1=0.05, CD2_2=0.1, inherit=FALSE)
plot(exdata_start,cex=2)
points(exdata_end,col='red')

#The residuals should be very small (in the noise of double precision arithmetic):

plot(density(exdata_start[,1]-exdata_end[,1]))
lines(density(exdata_start[,2]-exdata_end[,2]),col='red')

#Different world projections (ones not listed here do not work well across large angles).
#Note that the lines functions are picking up the last used WCS.

library(magicaxis)

world.grid=expand.grid(seq(0,360,by=10), seq(-90,90,by=10))

magplot(Rwcs_s2p(world.grid, CTYPE1='RA---ARC', CTYPE2='DEC--ARC'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---ZEA', CTYPE2='DEC--ZEA'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CYP', CTYPE2='DEC--CYP'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CEA', CTYPE2='DEC--CEA'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CAR', CTYPE2='DEC--CAR'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---MER', CTYPE2='DEC--MER'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---SFL', CTYPE2='DEC--SFL'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---PAR', CTYPE2='DEC--PAR'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---PCO', CTYPE2='DEC--PCO'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---MOL', CTYPE2='DEC--MOL'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---AIT', CTYPE2='DEC--AIT'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---TSC', CTYPE2='DEC--TSC'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---CSC', CTYPE2='DEC--CSC'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---QSC', CTYPE2='DEC--QSC'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---HPX', CTYPE2='DEC--HPX'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')
magplot(Rwcs_s2p(world.grid, CTYPE1='RA---XPH', CTYPE2='DEC--XPH'), pch='.', cex=3)
lines(Rwcs_s2p(cbind(-50:50,-50:50)), col='red')

#check pixscale

\dontrun{
library(Rfits)
file_image=system.file('extdata', 'image.fits', package = "Rfits")
image=Rfits_read_image(file_image)
Rwcs_pixscale(image$keyvalues)

#Rwcs also supports more complex distortion headers (e.g. SIP)
#but this removes user control

#Read example image with SIP distortion terms
file_image = system.file('extdata', 'sipsample.fits', package = "Rwcs")
temp_image = Rfits_read_image(file_image)

#Using it this way the SIP terms are being ignored. If we pass the raw header
#these can be used:

longhead = Rfits_read_header_raw(file_image)

#The difference from the version where we ignore the terms is very small! It is more
#obvious when inspecting the same pixel:

Rwcs_p2s(1, 1, temp_image$keyvalues, coord.type='sex')
Rwcs_p2s(1, 1, header=longhead, coord.type='sex')
}

#test if position is in frame

\dontrun{
library(Rfits)
file_image = system.file('extdata', 'image.fits', package = "Rfits")
temp_image = Rfits_read_image(file_image)

Rwcs_in_image(352.2917, -31.825, keyvalues = temp_image$keyvalues)
Rwcs_in_image(352.2917 + 1, -31.825, keyvalues = temp_image$keyvalues)

#A single pixel of overlap:
Rwcs_overlap(temp_image[1:200,1:200]$keyvalues, temp_image[200:356,200:356]$keyvalues)
#No pixel of overlap:
Rwcs_overlap(temp_image[1:200,1:200]$keyvalues, temp_image[201:356,200:356]$keyvalues)

#Show the overlap plot:
Rwcs_overlap(temp_image[100:200,100:200]$keyvalues, temp_image$keyvalues, plot=TRUE)
#Swapping the order still works but changes the sense of the plot since ref is in test:
Rwcs_overlap(temp_image$keyvalues, temp_image[100:200,100:200]$keyvalues, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
