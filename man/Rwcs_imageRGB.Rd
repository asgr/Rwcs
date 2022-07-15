\name{Rwcs_imageRGB}
\alias{Rwcs_imageRGB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RBG Images with WCS
}
\description{
Create 3-band colour RGB images with pretty WCS projection labels.
}
\usage{
Rwcs_imageRGB(R, G, B, keyvalues_out = NULL, Rkeyvalues = NULL, Gkeyvalues = NULL,
  Bkeyvalues = NULL, dowarp = "auto", direction = "auto", boundary = "dirichlet",
  interpolation = "cubic", n, grid.col = "grey", grid.lty = 2, grid.lwd = 0.5,
  lab.col = "green", coord.type = "sex", margin = TRUE, loc.diff = c(0, 0),
  xlab = "Right Ascension", ylab = "Declination", mgp = c(2, 0.5, 0), mtline = 2,
  position = "topright", com.col = "green", com.length = 0.05, coord.axis = "auto",
  pretty = "auto", decorate = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
Numeric matrix; containing the red colour values to be plotted (NAs are allowed). If \option{R} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{R} and the correct header part is parsed to \option{Rkeyvalues}.
}
  \item{G}{
Numeric matrix; containing the red colour values to be plotted (NAs are allowed). If \option{G} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{G} and the correct header part is parsed to \option{Gkeyvalues}.
}
  \item{B}{
Numeric matrix; containing the red colour values to be plotted (NAs are allowed). If \option{B} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{B} and the correct header part is parsed to \option{Bkeyvalues}.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS. This is the target WCS projection that the input images will be mapped onto.
}
  \item{Rkeyvalues}{
List; input R header values to be used for the WCS. This should be the header WCS that matches \option{R}.
}
  \item{Gkeyvalues}{
List; input G header values to be used for the WCS. This should be the header WCS that matches \option{G}.
}
  \item{Bkeyvalues}{
List; input B header values to be used for the WCS. This should be the header WCS that matches \option{B}.
}
  \item{dowarp}{
Image warping flag, either TRUE, FALSE or 'auto' (default). If TRUE then the images will be warped onto the \option{keyvalues_out} WCS scheme, if FALSE they will not be (which will break if the WCSs are not truly the same and the iamges are pixel matched). 'auto' tries to detect if the headers are all the same, and if they are then it will set \option{dowarp}=FALSE, otherwise it will set to TRUE. This can be fooled by small (and unimportant) differences in the headers, e.g. they have different SWarp dates etc, since it just checks if the text differs anywhere. If in doubt set this manually.
}
  \item{direction}{
Only used if re-mapping via \code{\link{Rwcs_warp}}. "auto" (default), "forward" or "backward", see \code{imwarp}. Since it is usally better to go from the higher resolution image and map this onto the lower resolution grid, "auto" selects the better direction given the pixel scales recovered from the header information.
}
  \item{boundary}{
Only used if re-mapping via \code{\link{Rwcs_warp}}. Boundary conditions: "dirichlet", "neumann", "periodic" (default "dirichlet"), see \code{imwarp}.
}
  \item{interpolation}{
Only used if re-mapping via \code{\link{Rwcs_warp}}. "nearest", "linear", "cubic" (default "linear"), see \code{imwarp}.
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
  \item{position}{
Relative position of the compass bar. This argument is used like the 'legend' function. Specify one of 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'bottomright' and 'centre'. 
}
  \item{com.col}{
Colour of the compass.
}
  \item{com.length}{
Length of the edges of the arrow head.
}
  \item{coord.axis}{
Integer vector; specifies whether the RA and Dec axes should be 1 or 2 (i.e. x or y axis). The default 'auto' tries to guess based on the header information (and this usually works okay). Otherwise \option{coord.axis}=c(1,2) would be sensible for a N/S vertical aligend frame and \option{coord.axis}=c(2,1) would be sensible for a E/W vertically aligned frame. When the orientation is between the two then it is not always obvious which will work better.
}
  \item{pretty}{
If 'auto' then it will try to compute the prettiest grid and label scaling *usually this is pretty good). Otherwise set to 1 to be degrees major ticks, 60 for minutes (DMS or HMS) and 3600 for seconds (DMS or HMS).
}
  \item{decorate}{
Logical; should image decorations be shown (like axes labels, tick, grid line and a compass via \code{\link{Rwcs_image}}). Not showing these might seem odd (since you could use \code{\link{magimageRGB}} in most cases), but it does allow you to generate plain WCS matched images easily.  
}
  \item{\dots}{
\code{\dots} are passed to \code{\link{magimageRGB}}.
}
}
\details{
\code{Rwcs_imageRGB} can be used in a few different ways, but the recommended route is to supply combined image and header list objects to the \option{R}, \option{G} and \option{B} arguments. If \option{keyvalues_out} is not supplied then by default it will project the three images onto the first available WCS header it finds searching \option{R}, \option{G} and \option{B} in order. If you are happy to map the colour image into the \option{R} WCS then you do not need to supply \option{keyvalues_out} at all. Note if image remapping is required (the images are different sizes or the headers differ) then the \code{imager} library will need to be installed in order to use \code{\link{Rwcs_warp}}.
}
\value{
Outputs the final image list x and y vectors and R/G/B matrices. Generally run for the side effect of producing rapid and well-scaled image plots.
}
\author{
Aaron Robotham
}
\note{
By R convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y. This differs to the FITS convention of pixel mid points being integer values. As such the R [x,y] = FITS [x-0.5,y-0.5]. This rarely matters too much in practice, but for accurate overlays you will want to get it right.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Rwcs_image}}, \code{\link{Rwcs_warp}}, \code{\link{magimageRGB}}
}
\examples{
library(Rfits)

VISTA_K = Rfits_read_image(system.file("extdata", 'VISTA_K.fits', package="magicaxis"))
VST_r = Rfits_read_image(system.file("extdata", 'VST_r.fits', package="magicaxis"))
GALEX_NUV = Rfits_read_image(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))

Rwcs_imageRGB(VISTA_K, VST_r, GALEX_NUV)
}

\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
