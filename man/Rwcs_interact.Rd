\name{Rwcs_interact}
\alias{Rwcs_interact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Interaction
}
\description{
Allows dynamic zooming and scrolling on an image, and also the extraction of coordinate information.
}
\usage{
Rwcs_interact(mode = "zoom", coord.type = "sex", col = "green", point.pch = 4,
  point.cex = 5, ...)
}
\arguments{
  \item{mode}{
Character scalar; either 'zoom' (the default, which allows zooming and scrolling) or 'coord' (which extracts the coordinates of locations on an image).
}
  \item{coord.type}{
Character scalar; the format of the output RA and Dec. Allowed options are 'deg' for degress (default) and 'sex' for sexigesimal (i.e. HMS for RA and DMS for Deg).
}
  \item{col}{
Scalar; the colour of information plotted on top of the image 
}
  \item{point.pch}{
Integer scalar; the \option{pch} point type.
}
  \item{point.cex}{
Integer scalar; the \option{cex} point type.
}
  \item{\dots}{
Extra arguments passed to \code{\link{Rwcs_image}}.
}
}
\details{
For this function to work \code{\link{Rwcs_image}} (or the plot method) will need to have been run with \option{interactive} = TRUE set. This will save the reference image as a hidden object in the global environment, allowing us to rapidly interact with it.

When \option{mode} = 'zoom' the basic interaction is the following: 1 click and then hitting ESC will move you around the scene (centering on the clicked point) and print the coordinate of the click to the terminal; 2 clicks will zoom in (first click is left of the second) or zoom out (first click is right of the second); 0 click and then hitting ESC will exit interactive mode. 
}
\value{
If \option{mode} = 'coord' then the output will be a two column matrix (RA, Dec) with the format as specified by \option{coord.type}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{Rwcs_image}}
}
\examples{
\dontrun{
library(Rfits)

file_image = system.file('extdata', 'image.fits', package = "Rfits")
temp_image = Rfits_read_image(file_image)
Rwcs_image(temp_image, interactive=TRUE)
Rwcs_interact() #this will let you interact with the image until you hit ESC with 0 click
}
}

