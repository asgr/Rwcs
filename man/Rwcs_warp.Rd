\name{Rwcs_warp}
\alias{Rwcs_warp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remap Image WCS via Warping
}
\description{
Remaps an input projection system to a different target WCS.
}
\usage{
Rwcs_warp(image_in, keyvalues_out = NULL, keyvalues_in = NULL, dim_out = NULL,
direction = "auto", boundary = "dirichlet", interpolation = "cubic", doscale = TRUE,
plot = FALSE, header_out = NULL, header_in = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_in}{
Numeric matrix; required, the image we want to warp. If \option{image} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{image} and the correct header part is parsed to \option{keyvalues}.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS. This is the target WCS projection that \option{image_in} will be mapped onto.
}
  \item{keyvalues_in}{
List; input header values to be used for the WCS. This should be the header WCS that matches \option{image_in}.
}
  \item{dim_out}{
Integer vector; this defines the desired dimensions of the output image. If this is not provided then the output image is made to be the same size as the NAXIS1 and NAXIS2 arguments taken from \option{header_out} (which is usually what you will want TBH).
}
  \item{direction}{
"auto" (default), "forward" or "backward", see \code{imwarp}. Since it is usally better to go from the higher resolution image and map this onto the lower resolution grid, "auto" selects the better direction given the pixel scales recovered from the header information.
}
  \item{boundary}{
boundary conditions: "dirichlet", "neumann", "periodic" (default "dirichlet"), see \code{imwarp}
}
  \item{interpolation}{
"nearest", "linear", "cubic" (default "linear"), see \code{imwarp}
}
  \item{doscale}{
Logical; if TRUE (default) then the image is scaled by the relative change the change in flux per pixel found in a refernece normalisation image. This is usually what you want when converting images between different WCS (since we really want to make sure we conserve overall flux). If FALSE then nothing is done to rescale the image. This is useful when the input \option{image_in} is something like a segmentation map, where the actual values matter.
}
  \item{plot}{
Logical; should a \code{\link{Rwcs_image}} plot of the output be generated? 
}
  \item{header_out}{
To aid backwards compatibility, you can also pass in a 1D character vector of header values, e.g. the \option{hdr} component of \option{readFITS} from the \code{FITSio} package. \option{keyvalues_out} must be left as NULL in this case.
}
  \item{header_in}{
To aid backwards compatibility, you can also pass in a 1D character vector of header values, e.g. the \option{hdr} component of \option{readFITS} from the \code{FITSio} package. \option{keyvalues_in} must be left as NULL in this case.
}
  \item{\dots}{
Dots are parsed to \code{\link{Rwcs_image}} (only relevant if \option{plot}=TRUE).
}  
}
\details{
The function allows for arbitrary WCS remapping from any \code{wcslib} supported projection (see \code{\link{Rwcs}}). The process internally does the following:

\itemize{
  \item Rwcs_p2s; maps the input image to RA and Dec per pixel using the input header
  \item Rwcs_s2p; maps the pixel RA and Dec coordinates onto the desired output pixel grid using the output header
}  

There are a few different ways to consider doing pixel remapping (or warping). The main question is usually whether to operate on a forward or backwards manner (see \code{imwarp}). Backwards mapping (the default) finds the best interpolation of every output pixel in the input image, whereas forward mapping finds the best interpolation of every input pixel on the output image. Backwards mapping usually provides the fewest artefacts (hence it is the default), but in general if the input image is higher resolution then you might prefer forward mapping and vica-versa.

The actual warping is done using the \code{imwarp} function in the \code{imager} package, so this will need to be installed prior to use (it is available on CRAN).
}
\value{
A list containing:

  \item{image}{Numeric matrix; the remapped image using the target WCS.}
  \item{keyvalues}{The target \option{keyvalues_out}.}
  \item{header}{The target \option{header_out}.}
}
\author{
Aaron Robotham
}
\note{
This function uses a bi-cubic interpolation scheme by default. It should approximately conserve the flux in \option{image_in}, but this is not guaranteed. The 'linear' interpolation scheme may be closer to conserving flux in images with sharp features. The conservation is usually good to about 0.1\% (i.e. 0.01 mag). If you require better conservation than this then alternative schemes (e.g. SWarp) will be required.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Rwcs}}
}
\examples{
\dontrun{
library(Rfits)
GALEX_NUV=Rfits_read_image(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))
VST_r=Rfits_read_image(system.file("extdata", 'VST_r.fits', package="magicaxis"))
VISTA_K=Rfits_read_image(system.file("extdata", 'VISTA_K.fits', package="magicaxis"))

Rwcs_warp(VST_r, GALEX_NUV$keyvalues, plot=TRUE)
Rwcs_warp(VISTA_K, GALEX_NUV$keyvalues, plot=TRUE)

Rwcs_warp(GALEX_NUV, VST_r$keyvalues, plot=TRUE)
Rwcs_warp(VISTA_K, VST_r$keyvalues, plot=TRUE)

Rwcs_warp(GALEX_NUV, VISTA_K$keyvalues, plot=TRUE)
Rwcs_warp(VST_r, VISTA_K$keyvalues, plot=TRUE)

#Check we can warp forwards and backwards correctly:

plot(GALEX_NUV)
Rwcs_warp(GALEX_NUV, GALEX_NUV$keyvalues, plot=TRUE)
Rwcs_warp(Rwcs_warp(GALEX_NUV, VST_r$keyvalues, interpolation='nearest'),
        GALEX_NUV$keyvalues, interpolation='nearest', plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ warp}% use one of  RShowDoc("KEYWORDS")
