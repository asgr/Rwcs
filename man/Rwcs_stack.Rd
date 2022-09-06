\name{Rwcs_stack}
\alias{Rwcs_stack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Warping and Stacking
}
\description{
Given a list of images with WCS, this function will re-project and stack all un-masked pixels.
}
\usage{
Rwcs_stack(image_list = NULL, inVar_list = NULL, mask_list = NULL, magzero_in = 0,
  magzero_out = 0, keyvalues_out = NULL, dim_out = NULL, cores = 4, return_all = FALSE,
  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; list of the input \code{Rfits_image} or \code{Rfits_pointer} objects. Masked pixels can be flagged with NA in the imDat component.
}
  \item{inVar_list}{
Numeric Vector or list of numeric matrices; the inverse variances. Either one value per \option{image_list} entry, or one pixel matched matrix per \option{image_list} entry. This is used to achieve S/N weighted stacking.
}
  \item{mask_list}{
List of numeric matrices; pixel matched matrix per \option{image_list} entry. Anything non-zero is considered masked. This can be used in conjuction with setting pixels in \option{image_list} to NA.
}
  \item{magzero_in}{
Numeric vector; the input mag-zero points. If length 1 then it is assumed all input frames have the same mag-zero point.
}
  \item{magzero_out}{
Numeric scalar; the output mag-zero point desired. A nice \option{magzero_out} is perhaps 23.9, because this means all measured fluxed are in micro-Jansky (rather than arbitrary ADUs). \option{magzero_out} = 8.9 gives Jansky; 16.4 gives milli-Jansky and 31.4 gives nano-Jansky.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS. This is the target WCS projection that each image in \option{image_list} will be mapped onto.
}
  \item{dim_out}{
Integer vector; this defines the desired dimensions of the output image. If this is not provided then the output image is made to be the same size as the NAXIS1 and NAXIS2 arguments taken from \option{header_out} (which is usually what you will want TBH).
}
  \item{cores}{
Integer scalar; the number of cores to run on.
}
  \item{return_all}{
Logical; should the pre stack projected images and inverse variances be returns also? Note if this is TRUE then the stacking can no longer be computed in batch mode, so be sure that you have enough memory to hold all pre stack projected images and inverse variances in memory (where note these will each have the dimensions of the desired output stacked image, which can be very large).
}
  \item{\dots}{
Other arguments to pass into \code{\link{Rwcs_warp}}.
}
}
\details{
The warping code used here is effectively \code{\link{Rwcs_warp}}, and the stacking part is very similar to that used in \code{profoundMakeStack}.
}
\value{
The output is an object of class 'ProMo' containing:

\item{image}{The stacked image, of class Rfits_image.}
\item{weight}{The number of available un-masked pixels in every output pixel of \option{image}. Class Rfits_image.}
\item{inVar}{The stacked inverse variance of the output \option{image}. Class Rfits_image.}

And optionally (if \option{returnall} = TRUE):

\item{image_pre_stack}{The projected pre stack images.}
\item{inVar_pre_stack}{The projected pre stack inverse variance.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{Rwcs_warp}}
}
\examples{
# Nothing here yet!
}
\concept{ ~warp }
\concept{ ~stack }

