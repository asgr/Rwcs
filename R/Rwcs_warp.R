Rwcs_warp = function (image_in, keyvalues_out=NULL, keyvalues_in = NULL, dim_out = NULL, 
          direction = "auto", boundary = "dirichlet", interpolation = "cubic", 
          doscale = TRUE, plot = FALSE, header_out = NULL, header_in = NULL, ...) 
{
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", 
         call. = FALSE)
  }
  if(any(names(image_in) == 'imDat') & is.null(keyvalues_in)){
    keyvalues_in = image_in$keyvalues
    header_in = image_in$hdr
    image_in = image_in$imDat
  }else if(any(names(image_in) == 'imDat') & !is.null(keyvalues_in)){
    image_in = image_in$imDat
  }
  if(any(names(image_in) == "image") & is.null(keyvalues_in)){
    keyvalues_in = image_in$keyvalues
    header_in = image_in$header
    image_in = image_in$image
  }else if(any(names(image) == "image") & !is.null(keyvalues_in)){
    image_in = image_in$image
  }
  if(is.character(header_out) & is.null(keyvalues_out)){
    if(requireNamespace("Rfits", quietly = TRUE)){
      keyvalues_out = Rfits::Rfits_hdr_to_keyvalues(header_out)
    }else{
      stop("The Rfits package is need to process the header_out. Install from GitHub asgr/Rfits.")
    }
  }
  if(is.character(header_in) & is.null(keyvalues_in)){
    if(requireNamespace("Rfits", quietly = TRUE)){
      keyvalues_in = Rfits::Rfits_hdr_to_keyvalues(header_in)
    }else{
      stop("The Rfits package is need to process the header_in. Install from GitHub asgr/Rfits.")
    }
  }
  if (!is.null(keyvalues_out) & is.null(dim_out)){
    dim_out = c(keyvalues_out$NAXIS1, keyvalues_out$NAXIS2)
  }else{
    stop('Missing NAXIS1 / NAXIS2 in header keyvalues! Specify dim_out.')
  }
  if (interpolation == "nearest") {
    warpoffset = 0.5
  }
  else {
    warpoffset = 0
  }
  .warpfunc_in2out = function(x, y) {
    radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_in, header = header_in)
    xy_out = Rwcs_s2p(radectemp, keyvalues = keyvalues_out, header = header_out)
    return = list(x = xy_out[, 1] + warpoffset, y = xy_out[,2] + warpoffset)
  }
  .warpfunc_out2in = function(x, y) {
    radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_out, header = header_out)
    xy_out = Rwcs_s2p(radectemp, keyvalues = keyvalues_in, header = header_in)
    return = list(x = xy_out[, 1] + warpoffset, y = xy_out[,2] + warpoffset)
  }
  image_out = matrix(0, max(dim(image_in)[1], dim_out[1]),max(dim(image_in)[2], dim_out[2]))
  image_out[1:dim(image_in)[1], 1:dim(image_in)[2]] = image_in
  pixscale_in = Rwcs_pixscale(keyvalues=keyvalues_in, header=header_in)
  pixscale_out = Rwcs_pixscale(keyvalues=keyvalues_out, header=header_out)
  if (doscale) {
    scale = pixscale_out^2/pixscale_in^2
  }
  else {
    scale = 1L
  }
  if (direction == "auto") {
    if (pixscale_in < pixscale_out) {
      direction = "forward"
    }
    if (pixscale_in >= pixscale_out) {
      direction = "backward"
    }
  }
  if (direction == "forward") {
    out = imager::imwarp(im = imager::as.cimg(image_out), 
                         map = .warpfunc_in2out, direction = direction, coordinates = "absolute", 
                         boundary = boundary, interpolation = interpolation)
  }
  if (direction == "backward") {
    out = imager::imwarp(im = imager::as.cimg(image_out), 
                         map = .warpfunc_out2in, direction = direction, coordinates = "absolute", 
                         boundary = boundary, interpolation = interpolation)
  }
  output = list(image = as.matrix(out)[1:dim_out[1], 1:dim_out[2]] * 
                  scale, keyvalues = keyvalues_out, header = header_out)
  if (plot) {
    Rwcs_image(output, ...)
  }
  return = output
}