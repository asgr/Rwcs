\name{Rwcs_s2p}
\alias{Rwcs_s2p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High Level WCSlib Interface
}
\description{
Offers basic high level access to the C \code{wcslib} library for projecting and deprojecting [x,y] and [RA,Dec] image data using a target WCS.
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
List; header values to be used for the WCS. the \option{\dots} can be used to over-ride any values for convenience.
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
Optional keyvalies the overwrite the option{keyvalues} entries with. CRVAL1=0, CRVAL2=0, CRPIX1=0, CRPIX2=0, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1, PV1=0, PV2=0, CTYPE1='RA---TAN', CTYPE2='DEC--TAN'
}
}
\details{
The function input and output is deliberately simple. The important thing is to pass in a senesible WCS to the \option{keyvalues} argument. This would be the same \option{keyvalues} output supplied by the \code{Rfits} package when reading in an image.
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

#The residuals should be very small (in the noice of double precision arithmetic):

plot(density(exdata_start[,1]-exdata_end[,1]))
lines(density(exdata_start[,2]-exdata_end[,2]),col='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
