# Rwcs (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/Rwcs/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

Core package providing a high level interface to the C **wcslib** library for projecting and deprojecting astronomical image coordinates. It converts between pixel [x,y] and sky [RA,Dec] coordinates using a target World Coordinate System (WCS), and provides functions for WCS-aware image display and annotation. Comes bundled with **wcslib** v7.1.

## Installation

### Build Tools

**Rwcs** requires compilation, so here is what you might need depending on your platform.

#### Linux Users

You know what you are doing. You do you!

#### Mac Users

You should not need to install separate compilers with any **R** after v4.0.0, but in case you are stuck on a museum version you can follow the extra instructions here:

[https://mac.r-project.org/tools/](https://mac.r-project.org/tools/)

#### Windows Users

Windows users might need to install *Rtools*, which are available at [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/). You will know it is working because the following will not be empty:

```R
Sys.which("make")
```

### Getting Rwcs

Source installation from GitHub should be easy:

```R
install.packages('remotes')
remotes::install_github("asgr/Rwcs")
library(Rwcs)
```

If you have trouble with the above you can try:

```R
install.packages('devtools')
devtools::install_github("asgr/Rwcs")
library(Rwcs)
```

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **Rwcs**:

```R
install.packages(c('Rcpp', 'checkmate', 'celestial', 'magicaxis', 'foreach', 'doParallel', 'data.table'))
install.packages('remotes')
remotes::install_github("asgr/Rwcs")
```

Some functionality also uses **Rfits** (for reading FITS headers):

```R
install.packages('remotes')
remotes::install_github("asgr/Rfits")
```

Assuming you have installed all of the packages that you need, you should now be able to load **Rwcs** within **R** with the usual:

```R
library(Rwcs)
```

## Code Example

```R
library(Rwcs)
library(Rfits)

# Read a FITS image with its WCS header
file_image = system.file('extdata', 'image.fits', package = "Rfits")
image = Rfits::Rfits_read_image(file_image)

# Convert sky coordinates (RA, Dec) to pixel coordinates
pixel = Rwcs_s2p(RA = 180.0, Dec = 0.0, keyvalues = image$keyvalues)

# Convert pixel coordinates back to sky coordinates
sky = Rwcs_p2s(x = pixel[,'x'], y = pixel[,'y'], keyvalues = image$keyvalues)

# Display the image with a WCS grid and labels
Rwcs_image(image)
```

To find more long-form examples please check the vignettes provided. You can browse these with:

```R
browseVignettes('Rwcs')
```

## Motivation

**Rwcs** was created to provide a robust and convenient high level R interface to the widely-used **wcslib** library. It handles the full range of sky projections supported by **wcslib** and integrates naturally with **Rfits** for reading and writing FITS files. It is used throughout the **asgr** ecosystem of packages (e.g. **ProPane**, **Rfits**) wherever WCS coordinate transformations are needed.

## Contributors

Code:

Aaron Robotham  
Rodrigo Tobar

## License

LGPL-3
