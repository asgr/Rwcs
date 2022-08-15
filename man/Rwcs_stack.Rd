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
  magzero_out = 0, keyvalues_out = NULL, dim_out = NULL, cores = 4, ...)
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
Numeric scalar; the output mag-zero point desired.
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
  \item{\dots}{
Other arguments to pass into \code{\link{Rwcs_warp}}.
}
}
\details{
The warping code used here is effectively \code{\link{Rwcs_warp}}, and the stacking part is very similar to that used in \code{profoundMakeStack}.
}
\value{
\item{image}{The stacked image, of class Rfits_image.}
\item{weight}{The number of available un-masked pixels in every output pixel of \option{image}. Class Rfits_image.}
\item{inVar}{The stacked inverse variance of the output \option{image}. Class Rfits_image.}
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

