\name{Rwcs_tweak}
\alias{Rwcs_tweak}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tweak the Translation Astrometry of an Image
}
\description{
Function to rapidly find translation tweaks directly in the image domain (i.e. it does not require a catalogue to be generated first). This is generally pretty fast, especially when running in integer shift mode.
}
\usage{
Rwcs_tweak(image_ref, image_pre_fix, delta_max = 3, quan_cut = 0.99, Nmeta = 3, cores = 4,
  shift_int = TRUE, return_image = FALSE, direction = "backward", final_centre = TRUE,
  verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_ref}{
Numeric matrix; static reference image (must be exactly the same dimensions as \option{image_pre_fix})
}
  \item{image_pre_fix}{
Numeric matrix; image to be tweaked (must be exactly the same dimensions as \option{image_ref})
}
  \item{delta_max}{
Numeric scalar; maximum number of pixels to explore in +/-x and +/- y.
}
  \item{quan_cut}{
Numeric scalar; bright pixel thesholding to apply to \option{image_ref} and \option{image_pre_fix}. Generally you want this quite high so that only cores of bright objects are compared.
}
  \item{Nmeta}{
Integer scalar; how many shifts of the search box can be applied. This is only relevant when \option{shift_int} is TRUE.
}
  \item{cores}{
Integer scalar; the number of cores to run on.
}
  \item{shift_int}{
Logical; should the shift applied be purely integer (TRUE) or fractional pixel (FALSE). Integer shifts are much faster to compute, and are better at conserving flux and not introducing pixel covariance, so think carefully before setting to FALSE.
}
  \item{return_image}{
Logical; should the tweaked image be returned as an output. You might only want the solution (for applying to another process), in which case set to FALSE for speed.
}
  \item{direction}{
Character scalar; "auto" (default), "forward" or "backward", see \code{imwarp}. Since it is usally better to go from the higher resolution image and map this onto the lower resolution grid, "auto" selects the better direction given the pixel scales recovered from the header information.
}
  \item{final_centre}{
Logical; should a final round of costs be computed once the best tweak is found. This is only relevant if \option{shift_int} = TRUE, and can be useful if you want to explore the output \option{cost_mat}.
}
  \item{verbose}{
Logical; should verbose outputs be printed as it runs?
}
}
\details{
Image tweaking is useful before stacking with e.g. \code{\link{Rwcs_stack}} if images are slightly misaligned in translation, but not rotated or distorted otherwise. This is quite common with HST and JWST images when guide starts are switched between observing blocks.
}
\value{
When \option{shift_int} is TRUE:

\item{optim_out}{The output of the actual optimisation (includes the best par).}
\item{image_post_fix}{The tweaked image (if \option{return_image} is TRUE).}
\item{time}{The total running time in seconds.}
\item{cost_mat}{The cost outputs for the final iteration. Thie has dimensions x = delta_x, y = delta_y, z = cost.}
\item{at_lim}{Logical; has the solution hit a limit. If so maybe re-run with larger \option{delta_max} and/or larger \option{Nmeta}.}

When \option{shift_int} is FALSE:

\item{optim_out}{The output of the actual optimisation (includes the best par).}
\item{image_post_fix}{The tweaked image (if \option{return_image} is TRUE).}
\item{time}{The total running time in seconds.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{Rwcs_tran}}
}
\examples{
library(Rfits)

image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"), header=FALSE)

tweak = Rwcs_tweak(image, Rwcs_tran(image, 6, 6), cores=1, verbose=FALSE)
tweak$optim_out$par
}

\concept{ ~tweak }

