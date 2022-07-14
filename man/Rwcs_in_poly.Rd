\name{Rwcs_in_poly}
\alias{Rwcs_in_poly}
\alias{Rwcs_overlap_poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polygon Intersection Tests
}
\description{
Tests whether points exist inside generic concave polygon, or whether two concave polygons overlap at all.
}
\usage{
Rwcs_in_poly(x, y, poly_x, poly_y)
Rwcs_overlap_poly(poly_x1, poly_y1, poly_x2, poly_y2, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; x points to test if in polygon.
}
  \item{y}{
Numeric vector; y points to test if in polygon.
}
  \item{poly_x}{
Numeric vector; x vertices of polygon. This need not be closing (so 4 points define square, not 5).
}
  \item{poly_y}{
Numeric vector; y vertices of polygon. This need not be closing (so 4 points define square, not 5).
}
  \item{poly_x1}{
Numeric vector; x vertices of polygon 1 to test for overlap with polygon 2. This need not be closing (so 4 points define square, not 5).
}
  \item{poly_y1}{
Numeric vector; y vertices of polygon 1 to test for overlap with polygon 2. This need not be closing (so 4 points define square, not 5).
}
  \item{poly_x2}{
Numeric vector; x vertices of polygon 2 to test for overlap with polygon 1. This need not be closing (so 4 points define square, not 5).
}
  \item{poly_y2}{
Numeric vector; y vertices of polygon 2 to test for overlap with polygon 1. This need not be closing (so 4 points define square, not 5).
}
  \item{plot}{
Logical; make a simple plot showing the polygons.  
}
}
\details{
Uses rapid method to test for polygon overlap. This should scale well to large number of points and vertices.

Note in most cases (bar machine precision) exactly touching means "inside" or "overlapping".
}
\value{
\code{Rwcs_in_poly}: returns logical vector of length \option{x}, where each \option{x}/\option{y} are flagged for being inside polygon.

\code{Rwcs_overlap_poly}; returns scalar logical specifying if polygon 1 and polygon 2 overlap at all.
}
\references{
Various geometry things I Googled.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{Rwcs_in_image}}, \code{\link{Rwcs_overlap}}
}
\examples{
#Define a concave star polygon:
star_x = c(0,0.3,0,0.5,1,0.7,1,0.5)
star_y = c(0,0.5,1,0.7,1,0.5,0,0.3)

plot(seq(0,1,len=11), seq(0.4,0.6,len=11), xlim=c(0,1), ylim=c(0,1))
polygon(star_x, star_y)
Rwcs_in_poly(seq(0,1,len=11), seq(0.4,0.6,len=11), star_x, star_y)

Rwcs_overlap_poly(star_x, star_y, star_x+0.5, star_y+0.5, plot=TRUE)
Rwcs_overlap_poly(star_x, star_y, star_x+1, star_y+1, plot=TRUE) #touching
Rwcs_overlap_poly(star_x, star_y, star_x+1.01, star_y+1.01, plot=TRUE)
Rwcs_overlap_poly(star_x, star_y, star_x+0.5, star_y+0.7, plot=TRUE) #touching
Rwcs_overlap_poly(star_x, star_y, star_x+0.5, star_y+0.8, plot=TRUE)
}

