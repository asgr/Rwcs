\name{Rwcs_image}
\alias{Rwcs_image}
\alias{Rwcs_grid}
\alias{Rwcs_labels}
\alias{Rwcs_compass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WCS Image Decoration
}
\description{
These functions add various WCS information to R images. It is particularly designed to work in unison with \option{magimage} package. Can handle all projections supported by \code{wcslib}, although some might need more tweaking to display well. See \code{\link{Rwcs}}.
}
\usage{
Rwcs_image(image, keyvalues = NULL, n, grid.col = "grey", grid.lty = 2, grid.lwd = 0.5,
lab.col = "green", coord.type = "sex", margin = TRUE, loc.diff = c(0, 0),
xlab = "Right Ascension", ylab = "Declination", mgp = c(2, 0.5, 0), mtline = 2,
position = "topright", com.col = "green", com.length = 0.05, coord.axis = 'auto',
pretty = 'auto', header = NULL, add = FALSE, ...)

Rwcs_grid(keyvalues = NULL, n, grid.col = "grey", grid.lty = 2, grid.lwd = 0.5,
coord.type = "sex", loc.diff = c(0, 0), pretty= 'auto', header = NULL, ...)

Rwcs_labels(keyvalues = NULL, n, lab.col = "green", coord.type = "sex", margin = TRUE,
loc.diff = c(0, 0), xlab = "Right Ascension", ylab = "Declination", mgp = c(2, 0.5, 0),
mtline = 2, coord.axis='auto', pretty= 'auto', header = NULL, ...)

Rwcs_compass(keyvalues = NULL, position = "topright", com.col = "green",
com.length = 0.05, loc.diff = c(0, 0), header = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to decorate. If \option{image} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{image} and the correct header part is parsed to \option{keyvalues}.
}
  \item{keyvalues}{
List; header values to be used for the WCS.
}
  \item{n}{
The target number of major-axis sub-divisions. Will not necessarily be achieved.
}
  \item{grid.col}{
The colour of the overlaid grid lines.
}
  \item{grid.lty}{
The line type of the overlaid grid lines.
}
  \item{grid.lwd}{
The line width of the overlaid grid lines.
}
  \item{lab.col}{
The colour of the labels when \option{margin}=FALSE.
}
  \item{coord.type}{
Should the labels be drawn using degrees (deg) or colon delimited sexigesimal (sex).
}
  \item{margin}{
Should the labels be drawn in the outer margin region (default).
}
  \item{loc.diff}{
The pixel offset to apply. Only relevant if the image being plotted is a cutout from within a FITS legal image.
}
  \item{xlab}{
x axis name. If left as default either H:M:S or D:M:S (\option{coord.type}='sex') or deg (\option{coord.type}='deg') will be appended.
}
  \item{ylab}{
y axis name. If left as default either H:M:S or D:M:S (\option{coord.type}='sex') or deg (\option{coord.type}='deg') will be appended.
}
  \item{mgp}{
The margin line (in mex units) for the axis title, axis labels and axis line. This has different (i.e. prettier) defaults than R of c(2,0.5,0) rather than c(3,1,0). This pushes the numbers and labels nearer to the plot compared to the defaults. For \option{margin}=FALSE \option{mgp} = -\option{mgp} - 3, which has the effect of shifting the tick labels nicely inside the margin.
}
  \item{mtline}{
Number of lines separating axis name from axis. For \option{margin}=FALSE \option{mtline} = -\option{mtline}, which has the effect of shifting the axis labels nicely inside the margin.
}
\item{coord.axis}{
Integer vector; specifies whether the RA and Dec axes should be 1 or 2 (i.e. x or y axis). The default 'auto' tries to guess based on the header information (and this usually works okay). Otherwise \option{coord.axis}=c(1,2) would be sensible for a N/S vertical aligend frame and \option{coord.axis}=c(2,1) would be sensible for a E/W vertically aligned frame. When the orientation is between the two then it is not always obvious which will work better.
}
\item{pretty}{
If 'auto' then it will try to compute the prettiest grid and label scaling *usually this is pretty good). Otherwise set to 1 to be degrees major ticks, 60 for minutes (DMS or HMS) and 3600 for seconds (DMS or HMS).
}
  \item{position}{
Relative position of the compass bar. This argument is used like the 'legend' function. Specify one of 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'bottomright' and 'centre'. 
}
  \item{com.col}{
Colour of the compass.
}
  \item{com.length}{
Length of the edges of the arrow head.
}
  \item{header}{
To aid backwards compatibility, you can also pass in a 1D character vector of header values, e.g. the \option{hdr} component of \option{readFITS} from the \code{FITSio} package. \option{keyvalues} must be left as NULL in this case. If you want to use all distortion terms you can also pass in the raw header output of \code{Rfits_read_header_raw}.
}
  \item{add}{
Logical; should the target image be added to the current sky view?  
}
  \item{\dots}{
These are parsed to \code{\link{magimage}}, (\code{Rwcs_image}), \code{\link{lines}} (\code{Rwcs_grid}), \code{\link{axis}} (\code{Rwcs_labels}) or \code{\link{arrows}} (\code{Rwcs_compass}).
}
}
\details{
Most people will be content to use the higher level \code{Rwcs_image} function, which calls (in order) \code{magimage}, \code{Rwcs_grid}, \code{Rwcs_labels} and \code{Rwcs_compass}. For finer control you might want to execute these individual functions directly.
}
\value{
Outputs the final image list containing x,y and z (\code{Rwcs_image}). Generally run for the side effect of producing rapid and well-scaled image plots.
}
\author{
Aaron Robotham
}
\note{
By R convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y. This differs to the FITS convention of pixel mid points being integer values. As such the R [x,y] = FITS [x-0.5,y-0.5]. This rarely matters too much in practice, but for accurate overlays you will want to get it right.
}
\seealso{
\code{\link{magimage}}, \code{\link{Rwcs}}
}
\examples{
\dontrun{
#Recommended to use Rfits (will be supported longer term):

library(Rfits)
file_image = system.file('extdata', 'image.fits', package = "Rfits")
temp_image = Rfits_read_image(file_image)
Rwcs_image(temp_image)

#This will also work for the time being (perhaps not in the future):

library(FITSio)
temp_image = readFITS(file_image)
Rwcs_image(temp_image)

#we can see the curvature better if we ramp up the size of the image, but keep WCS:

temp_mat = list(x=seq(1,1e5,len=356), y=seq(1,1e5,len=356), z=temp_image$imDat)
Rwcs_image(temp_mat, temp_image$keyvalues)

#Now let's try some other projections:

tempkey = Rwcs_keypass(temp_image$keyvalues, CTYPE1 = 'RA---SIN', CTYPE2 = 'DEC--SIN')
Rwcs_image(temp_mat, keyvalues = tempkey)
tempkey = Rwcs_keypass(temp_image$keyvalues, CTYPE1 = 'RA---NCP', CTYPE2 = 'DEC--NCP')
Rwcs_image(temp_mat, keyvalues = tempkey)

#Rwcs also supports more complex distortion headers (e.g. SIP)
#but this removes user control

#Read example image with SIP distortion terms
file_image = system.file('extdata', 'sipsample.fits', package = "Rwcs")
temp_image = Rfits_read_image(file_image)
Rwcs_image(temp_image)

#Using it this way the SIP terms are being ignored. If we pass the raw header
#these can be used:

longhead = Rfits_read_header_raw(file_image)
Rwcs_image(temp_image$imDat, header=longhead)

#The difference from the version where we ignore the terms is very small! It is more
#obvious when inspecting the same pixel:

Rwcs_p2s(1, 1, temp_image$keyvalues, coord.type='sex')
Rwcs_p2s(1, 1, header=longhead, coord.type='sex')

#So we are looking at differences of about 1/10th of an arcsecond.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
