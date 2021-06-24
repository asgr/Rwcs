Rwcs_warp = function (image_in, keyvalues_out=NULL, keyvalues_in = NULL, dim_out = NULL, 
          direction = "auto", boundary = "dirichlet", interpolation = "cubic", 
          doscale = TRUE, plot = FALSE, header_out = NULL, header_in = NULL, dotightcrop = TRUE, ...) 
{
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", 
         call. = FALSE)
  }
  
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop("The Rfits package is needed!")
  }
  
  # if(any(names(image_in) == 'imDat') & is.null(keyvalues_in)){
  #   if(is.null(header_in)){
  #     keyvalues_in = image_in$keyvalues
  #     header_in = image_in$hdr
  #   }
  #   image_in = image_in$imDat
  # }else if(any(names(image_in) == 'imDat') & !is.null(keyvalues_in)){
  #   image_in = image_in$imDat
  # }
  # if(any(names(image_in) == "image") & is.null(keyvalues_in)){
  #   if(is.null(header_in)){
  #     keyvalues_in = image_in$keyvalues
  #     header_in = image_in$header
  #   }
  #   image_in = image_in$image
  # }else if(any(names(image) == "image") & !is.null(keyvalues_in)){
  #   image_in = image_in$image
  # }
  
  if(any(names(image_in)=='imDat') | any(names(image_in)=='image')){
    if(is.null(keyvalues_in)){
      keyvalues_in = image_in$keyvalues
    }
    
    if(is.null(header_in)){
      if(is.null(image_in$raw)){
        header_in = image_in$hdr
      }else{
        header_in = image_in$raw
      }
    }
    
    if(any(names(image_in)=='imDat')){
      image_in = image_in$imDat
    }
    
    if(any(names(image_in)=='image')){
      image_in = image_in$image
    }
  }
  
  if(is.character(header_out) & is.null(keyvalues_out)){
    if(length(header_out) > 1){
      if(requireNamespace("Rfits", quietly = TRUE)){
        keyvalues_out = Rfits::Rfits_hdr_to_keyvalues(header_out)
      }else{
        stop("The Rfits package is need to process the header_out. Install from GitHub asgr/Rfits.")
      }
    }
  }
  if(is.character(header_in) & is.null(keyvalues_in)){
    if(length(header_out) > 1){
      if(requireNamespace("Rfits", quietly = TRUE)){
        keyvalues_in = Rfits::Rfits_hdr_to_keyvalues(header_in)
      }else{
        stop("The Rfits package is need to process the header_in. Install from GitHub asgr/Rfits.")
      }
    }
  }
  
  if(is.null(keyvalues_out) & is.null(header_out)){
    keyvalues_out = options()$current_keyvalues
    header_out = options()$current_header
  }
  
  if (!is.null(keyvalues_out) & is.null(dim_out)){
    dim_out = c(keyvalues_out$NAXIS1, keyvalues_out$NAXIS2)
  }
  
  if(is.null(dim_out)){
    stop('Missing NAXIS1 / NAXIS2 in header keyvalues! Specify dim_out.')
  }
  
  if (interpolation == "nearest") {
    warpoffset = 0.5
  }else{
    warpoffset = 0
  }
  
  image_out = list(
    imDat = matrix(NA, max(dim(image_in)[1], dim_out[1]), max(dim(image_in)[2], dim_out[2])),
    keyvalues = keyvalues_out,
    hdr = Rfits::Rfits_keyvalues_to_hdr(keyvalues_out),
    header = Rfits::Rfits_keyvalues_to_header(keyvalues_out),
    raw = Rfits::Rfits_header_to_raw(Rfits::Rfits_keyvalues_to_header(keyvalues_out)),
    keynames = names(keyvalues_out),
    keycomments = as.list(rep('', length(keyvalues_out)))
  )
  names(image_out$keycomments) = image_out$keynames
  class(image_out) = c('Rfits_image', class(image_out))
  
  if(dotightcrop){
    BL = Rwcs_p2s(0, 0, keyvalues = keyvalues_in, header=header_in, pixcen='R')
    TL = Rwcs_p2s(0, dim(image_in)[2], keyvalues = keyvalues_in, header=header_in, pixcen='R')
    TR = Rwcs_p2s(dim(image_in)[1], dim(image_in)[2], keyvalues = keyvalues_in, header=header_in, pixcen='R')
    BR = Rwcs_p2s(dim(image_in)[1], 0, keyvalues = keyvalues_in, header=header_in, pixcen='R')
    corners = rbind(BL, TL, TR, BR)
    tightcrop = ceiling(Rwcs_s2p(corners, keyvalues = keyvalues_out, header=header_out, pixcen='R'))
    min_x = max(1L, min(tightcrop[,1]))
    max_x = max(min_x + dim(image_in)[1], range(tightcrop[,1])[2])
    min_y = max(1L, min(tightcrop[,2]))
    max_y = max(min_y + dim(image_in)[2], range(tightcrop[,2])[2])
    image_out = image_out[c(min_x, max_x), c(min_y, max_y)]
    keyvalues_out = image_out$keyvalues
    header_out = image_out$raw
  }else{
    min_x = 1L
    max_x = dim_out[1]
    min_y = 1L
    max_y = dim_out[2]
  }
  
  .warpfunc_in2out = function(x, y) {
    radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_in, header = header_in)
    xy_out = Rwcs_s2p(radectemp, keyvalues = keyvalues_out, header = header_out)
    return(list(x = xy_out[, 1] + warpoffset, y = xy_out[,2] + warpoffset))
  }
  .warpfunc_out2in = function(x, y) {
    radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_out, header = header_out)
    xy_out = Rwcs_s2p(radectemp, keyvalues = keyvalues_in, header = header_in)
    return(list(x = xy_out[, 1] + warpoffset, y = xy_out[,2] + warpoffset))
  }
  
  image_out$imDat[1:dim(image_in)[1], 1:dim(image_in)[2]] = image_in

  suppressMessages({
    pixscale_in = Rwcs_pixscale(keyvalues=keyvalues_in)
    pixscale_out = Rwcs_pixscale(keyvalues=keyvalues_out)
  })
  
  norm = matrix(1, max(dim(image_in)[1], dim(image_out)[1]),max(dim(image_in)[2], dim(image_out)[2]))
  
  if (direction == "auto") {
    if (pixscale_in < pixscale_out) {
      direction = "forward"
    }
    if (pixscale_in >= pixscale_out) {
      direction = "backward"
    }
  }
  if (direction == "forward") {
    out = imager::imwarp(im = imager::as.cimg(image_out$imDat), 
                         map = .warpfunc_in2out, direction = direction, coordinates = "absolute", 
                         boundary = boundary, interpolation = interpolation)
    if(doscale){
      renorm = imager::imwarp(im = imager::as.cimg(norm), 
                           map = .warpfunc_in2out, direction = direction, coordinates = "absolute", 
                           boundary = boundary, interpolation = interpolation)
      out = out / renorm
    }
  }
  if (direction == "backward") {
    out = imager::imwarp(im = imager::as.cimg(image_out$imDat), 
                         map = .warpfunc_out2in, direction = direction, coordinates = "absolute", 
                         boundary = boundary, interpolation = interpolation)
    if(doscale){
      renorm = imager::imwarp(im = imager::as.cimg(norm), 
                           map = .warpfunc_out2in, direction = direction, coordinates = "absolute", 
                           boundary = boundary, interpolation = interpolation)
      out = out / renorm
    }
  }
  
  # output = list(imDat = as.matrix(out)[1:dim(image_out)[1], 1:dim(image_out)[2]],
  #               keyvalues = keyvalues_out,
  #               hdr = hdr,
  #               header = header,
  #               raw = raw,
  #               keynames = names(keyvalues_out),
  #               keycomments = as.list(rep('', length(keyvalues_out)))
  #               )
  # 
  # names(output$keycomments) = output$keynames
  
  image_out$imDat[] = out
  
  image_out = image_out[c(1L - (min_x - 1L), dim_out[1] - (min_x - 1L)),c(1L - (min_y - 1L), dim_out[2] - (min_y - 1L))]
  
  if (plot) {
    Rwcs_image(image_out, ...)
  }
  return(invisible(image_out))
}
