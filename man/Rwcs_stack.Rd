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
Rwcs_stack(image_list = NULL, inVar_list = NULL, exp_list = NULL, mask_list = NULL,
  magzero_in = 0, magzero_out = 23.9, keyvalues_out = NULL, dim_out = NULL, cores = 4,
  Nbatch = cores, keep_extreme_pix = FALSE, doclip = FALSE, clip_tol = 100,
  clip_dilate = 0, clip_sigma = 5, return_all = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; required. List of the input \code{Rfits_image} or \code{Rfits_pointer} objects. Masked pixels can be flagged with NA in the imDat component.
}
  \item{inVar_list}{
Numeric vector or list of numeric matrices; optional. The inverse variances. Either one value per \option{image_list} entry, or one pixel matched matrix per \option{image_list} entry. This is used to achieve S/N weighted stacking. If not provided then images are stacked with equal weight.
}
  \item{exp_list}{
Numeric vector or list of numeric matrices; optional. The exposure time (in whatever units you like, but they should be consistent between input images). Either one value per \option{image_list} entry, or one pixel matched matrix per \option{image_list} entry. This is used to track maximal pixel exposure (so nothing is masked).  
}
  \item{mask_list}{
List of numeric matrices; optional. Pixel matched mask matrix per \option{image_list} entry. Anything non-zero is considered masked. This can be used in conjuction with setting pixels in \option{image_list} to NA.
}
  \item{magzero_in}{
Numeric vector; the input mag-zero points. If length 1 then it is assumed all input frames have the same mag-zero point.
}
  \item{magzero_out}{
Numeric scalar; the output mag-zero point desired. A nice \option{magzero_out} is perhaps 23.9 (the default), because this means all measured fluxed are in micro-Jansky (rather than arbitrary ADUs). \option{magzero_out} = 8.9 gives Jansky; 16.4 gives milli-Jansky and 31.4 gives nano-Jansky.
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
  \item{Nbatch}{
Integer scalar; the size of processing batches to be used internally. This can be larger than \option{cores}, but usually it makes more sense to set \option{Nbatch} = \option{cores} since that will minimise memory and maximise CPU usage. If \option{return_all} = TRUE or \option{doclip} = TRUE then \option{Nbatch} will be set to the number of images (since all pre-stack projected images will either be returned or need to be processed for clipped stacking). The fastest possible speed (but maximal memory impact) will be \option{Nbatch} = \option{cores} = length(\option{image_list}).
}
  \item{keep_extreme_pix}{
Logical; should a cold and hot pixel matrix be stored and returned? This can be useful for diagnosing and later masking bad cold and hot pixels from the stacked \option{image}. Potentially you can use these outputs to detect hot/cold biased pixels in the output \option{image} and replace them with the \option{cold}/\option{hot} pixel in the location, or alternatively just mask them entirely.
}
  \item{doclip}{
Logical; should the stack be cold/hot pixel clipped and restacked? Note \option{inVar_list} is required in this case, since the image variances are used to determine the clipping.  
}
  \item{clip_tol}{
Numeric scalar or vector; determines how should cold (\option{clip_tol}[1]) and hot (\option{clip_tol}[2]) pixels are clipped. The internal logic is (stack - cold).sqrt(inVar) > \option{clip_tol}[1] cold pixels and (hot - stack).sqrt(inVar) > \option{clip_tol}[2] hot pixels are clipped out for a subsequent re-stacking of the output \option{image}. If length is 1 then value is repeated for cold and hot pixels.
}
  \item{clip_dilate}{
Integer scalar; how much to dilate around cold/hot pixels you wish to mask out. This is only relevant when \option{doclip} is TRUE. The value of \option{clip_dilate} is the diameter of the circular kernal used, so it should be odd valued and >= 3 (which would be one pixel of dilation) to have any effect.
}
  \item{clip_sigma}{
Numeric scalar; cipping significance threshold for hot pixels. Only pixels where the \option{cold} pixel is below \option{clip_sigma} of the sky RMS will be flagged for hot pixel clipping. Generally you want this to be a value between 5-10, and it stops the cores of bright stars/galaxies being erroneously flagged as containing hot pixels. The trade off from setting this to higher values is you might start to lose the cores of bright objects, and the trade off setting this lower is you might not remove genuine hot pixels. In general this parameter should not need much fine-tuning.
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
\item{inVar}{If \option{inVar_list} is not NULL, the stacked inverse variance of the output \option{image}. Class Rfits_image when present.}
\item{exp}{If \option{exp_list} is not NULL, the total exposure time of the output \option{image}. Class Rfits_image when present. Note this will be unmasked, so values of 0 means no photon ever hit the output \option{image}, i.e. masked pixels can have postive exposures if they sit within the detector domain.}
\item{cold}{If \option{keep_extreme_pix} = TRUE, the extreme hot pixel value for the input \option{image_list} that formed the final stacked \option{image}. If \option{keep_extreme_pix} = FALSE this will be NULL. Class Rfits_image when present.}
\item{hot}{If \option{keep_extreme_pix} = TRUE, the extreme hot pixel value for the input \option{image_list} that formed the final stacked \option{image}. If \option{keep_extreme_pix} = FALSE this will be NULL. Class Rfits_image when present.}
\item{clip}{The number of frames clipped from the stacked \option{image}. This is only relevant when \option{doclip} is TRUE. Values larger than 1 can indicate that both cold and hot pixels have been clipped in this location, and also that the dilation process has masked a pixel multiple times.}


And optionally (if \option{returnall} = TRUE) these outputs are appended to the list:

\item{image_pre_stack}{The projected pre stack images.}
\item{inVar_pre_stack}{The projected pre stack inverse variance.}
\item{exp_pre_stack}{The projected pre stack exposure times.}
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

