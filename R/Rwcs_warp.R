.warpfunc_in2out = function(x, y, keyvalues_in=NULL, header_in=NULL, WCSref_in=NULL, keyvalues_out=NULL, raw_out=NULL, WCSref_out=NULL) {
  radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_in, header = header_in, WCSref = WCSref_in)
  xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], keyvalues = keyvalues_out, header = raw_out, WCSref = WCSref_out)
  return(xy_out)
}
.warpfunc_out2in = function(x, y, keyvalues_in=NULL, header_in=NULL, WCSref_in=NULL, keyvalues_out=NULL, raw_out=NULL, WCSref_out=NULL) {
  radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_out, header = raw_out, WCSref = WCSref_in)
  xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], keyvalues = keyvalues_in, header = header_in, WCSref = WCSref_out)
  return(xy_out)
}

Rwcs_warp = function (image_in, keyvalues_out=NULL, keyvalues_in = NULL, dim_out = NULL,
                      pixscale_out = NULL, pixscale_in = NULL,
          direction = "auto", boundary = "dirichlet", interpolation = "cubic", 
          doscale = TRUE, dofinenorm = TRUE, plot = FALSE, header_out = NULL, header_in = NULL, dotightcrop = TRUE,
          keepcrop=FALSE, WCSref_out = NULL, WCSref_in = NULL, magzero_out = NULL, magzero_in = NULL,
          blank=NA, warpfield=NULL, warpfield_return=FALSE, ...) 
{
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", 
         call. = FALSE)
  }
  
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop("The Rfits package is needed!")
  }
  
  if(inherits(image_in, 'Rfits_pointer')){
    image_in = image_in[,]
  }
  
  if(any(names(image_in)=='imDat') | any(names(image_in)=='image')){
    if(is.null(header_in)){
      if(is.null(keyvalues_in)){
        keyvalues_in = image_in$keyvalues
      }
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
  
  if(is.character(header_out)){
    if(requireNamespace("Rfits", quietly = TRUE)){
      if(length(header_out) == 1){
        keyvalues_out = Rfits::Rfits_header_to_keyvalues(Rfits::Rfits_raw_to_header(header_out))
      }else if(length(header_out) > 1){
        keyvalues_out = Rfits::Rfits_hdr_to_keyvalues(header_out)
      }
    }else{
      stop("The Rfits package is need to process the header_out. Install from GitHub asgr/Rfits.")
    }
  }
  if(is.character(header_in)){
    if(requireNamespace("Rfits", quietly = TRUE)){
      if(length(header_in) == 1){
        keyvalues_in = Rfits::Rfits_header_to_keyvalues(Rfits::Rfits_raw_to_header(header_in))
      }else if(length(header_in) > 1){
        keyvalues_in = Rfits::Rfits_hdr_to_keyvalues(header_in)
      }
    }else{
      stop("The Rfits package is need to process the header_in Install from GitHub asgr/Rfits.")
    }
  }
  
  if(is.null(keyvalues_out) & is.null(header_out)){
    keyvalues_out = options()$current_keyvalues
    header_out = options()$current_header
  }
  
  keyvalues_in = keyvalues_in[!is.na(keyvalues_in)]
  keyvalues_out = keyvalues_out[!is.na(keyvalues_out)]
  
  if (!is.null(keyvalues_out) & is.null(dim_out)){
    if(is.null(keyvalues_out$ZNAXIS1)){
      NAXIS1 = keyvalues_out$NAXIS1
    }else{
      NAXIS1 = keyvalues_out$ZNAXIS1
    }
    if(is.null(keyvalues_out$ZNAXIS2)){
      NAXIS2 = keyvalues_out$NAXIS2
    }else{
      NAXIS2 = keyvalues_out$ZNAXIS2
    }
    dim_out = c(NAXIS1, NAXIS2)
  }else{
    keyvalues_out$NAXIS1 = dim_out[1]
    keyvalues_out$NAXIS2 = dim_out[2]
  }
  
  if(is.null(dim_out)){
    stop('Missing NAXIS1 / NAXIS2 in header keyvalues! Specify dim_out.')
  }
  
  if(!is.null(magzero_in) & !is.null(magzero_out)){
    image_in = image_in*10^(-0.4*(magzero_in - magzero_out))
    keyvalues_out$MAGZERO = magzero_out
  }
  
  #seems I don't need this anymore- imager side bug presumably fixed (which is what it seemed to be)
  #leave commented out for now, until 100% happy it doesn't appear elsewhere
  
  #if (interpolation == "nearest") {
  #  warpoffset = 0.5
  #}else{
  #  warpoffset = 0
  #}
  
  # if(is.null(header_out)){
  #   raw = NULL
  # }else{
  #   raw = Rfits::Rfits_header_to_raw(Rfits::Rfits_keyvalues_to_header(keyvalues_out))
  # }
  
  raw_out = header_out
  
  if(is.null(raw_out)){
    header_out = Rfits::Rfits_keyvalues_to_header(keyvalues_out)
  }else{
    header_out = Rfits::Rfits_raw_to_header(raw_out)
  }
  
  if(dotightcrop){
    suppressMessages({
      BL = Rwcs_p2s(0, 0, keyvalues = keyvalues_in, header=header_in, pixcen='R', WCSref=WCSref_in)
      TL = Rwcs_p2s(0, dim(image_in)[2], keyvalues = keyvalues_in, header=header_in, pixcen='R', WCSref=WCSref_in)
      TR = Rwcs_p2s(dim(image_in)[1], dim(image_in)[2], keyvalues = keyvalues_in, header=header_in, pixcen='R', WCSref=WCSref_in)
      BR = Rwcs_p2s(dim(image_in)[1], 0, keyvalues = keyvalues_in, header=header_in, pixcen='R', WCSref=WCSref_in)
    })
    corners = rbind(BL, TL, TR, BR)
    tightcrop = ceiling(Rwcs_s2p(corners, keyvalues = keyvalues_out, header=raw_out, pixcen='R', WCSref=WCSref_out))
    min_x = max(1L, min(tightcrop[,1]))
    max_x = max(min_x + dim(image_in)[1] - 1L, range(tightcrop[,1])[2])
    min_y = max(1L, min(tightcrop[,2]))
    max_y = max(min_y + dim(image_in)[2] - 1L, range(tightcrop[,2])[2])
    #image_out = image_out[c(min_x, max_x), c(min_y, max_y)]
    # new code should be more efficient!
    
    if(!isTRUE(keyvalues_out$ZIMAGE)){
      keyvalues_out$NAXIS1 = max_x - min_x + 1L
      keyvalues_out$NAXIS2 = max_y - min_y + 1L
    }else{
      keyvalues_out$ZNAXIS1 = max_x - min_x + 1L
      keyvalues_out$ZNAXIS2 = max_y - min_y + 1L
    }
    
    keyvalues_out$CRPIX1 = keyvalues_out$CRPIX1 - min_x + 1L
    keyvalues_out$CRPIX2 = keyvalues_out$CRPIX2 - min_y + 1L
      
    image_out = list(
      imDat = matrix(c(blank,image_in[0]), max_x - min_x + 1L, max_y - min_y + 1L),
      keyvalues = keyvalues_out,
      hdr = Rfits::Rfits_keyvalues_to_hdr(keyvalues_out),
      header = header_out,
      raw = raw_out,
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))
    
    keyvalues_out = image_out$keyvalues
    if(!is.null(raw_out)){
      raw_out = image_out$raw
    }
  }else{
    image_out = list(
      imDat = matrix(c(blank,image_in[0]), max(dim(image_in)[1], dim_out[1]), max(dim(image_in)[2], dim_out[2])),
      keyvalues = keyvalues_out,
      hdr = Rfits::Rfits_keyvalues_to_hdr(keyvalues_out),
      header = header_out,
      raw = raw_out,
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))
    
    min_x = 1L
    max_x = dim_out[1]
    min_y = 1L
    max_y = dim_out[2]
  }
  
  # .warpfunc_in2out = function(x, y) {
  #   radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_in, header = header_in, WCSref = WCSref_in)
  #   xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], keyvalues = keyvalues_out, header = raw_out, WCSref = WCSref_out)
  #   return(list(x = xy_out[, 1], y = xy_out[,2]))
  # }
  # .warpfunc_out2in = function(x, y) {
  #   radectemp = Rwcs_p2s(x, y, keyvalues = keyvalues_out, header = raw_out, WCSref = WCSref_in)
  #   xy_out = Rwcs_s2p(radectemp[,1], radectemp[,2], keyvalues = keyvalues_in, header = header_in, , WCSref = WCSref_out)
  #   return(list(x = xy_out[, 1], y = xy_out[,2]))
  # }
  #sort out NULL issues
  #keyvalues_in=NULL, header_in=NULL, WCSref_in=NULL, keyvalues_out=NULL, raw_out=NULL
  # if(!is.null(keyvalues_in)){
  #   formals(.warpfunc_in2out)$keyvalues_in = keyvalues_in
  #   formals(.warpfunc_out2in)$keyvalues_in = keyvalues_in
  # }
  # if(!is.null(header_in)){
  #   formals(.warpfunc_in2out)$header_in = header_in
  #   formals(.warpfunc_out2in)$header_in = header_in
  # }
  # if(!is.null(WCSref_in)){
  #   formals(.warpfunc_in2out)$WCSref_in = WCSref_in
  #   formals(.warpfunc_out2in)$WCSref_in = WCSref_in
  # }
  # if(!is.null(keyvalues_out)){
  #   formals(.warpfunc_in2out)$keyvalues_out = keyvalues_out
  #   formals(.warpfunc_out2in)$keyvalues_out = keyvalues_out
  # }
  # if(!is.null(raw_out)){
  #   formals(.warpfunc_in2out)$raw_out = raw_out
  #   formals(.warpfunc_out2in)$raw_out = raw_out
  # }
  # if(!is.null(WCSref_out)){
  #   formals(.warpfunc_in2out)$WCSref_out = WCSref_out
  #   formals(.warpfunc_out2in)$WCSref_out = WCSref_out
  # }
  
  dim_min_x = min(dim(image_in)[1], dim(image_out$imDat)[1])
  dim_min_y = min(dim(image_in)[2], dim(image_out$imDat)[2])
  
  image_out$imDat[1:dim_min_x, 1:dim_min_y] = image_in[1:dim_min_x, 1:dim_min_y]

  suppressMessages({
    if(is.null(pixscale_in)){
      pixscale_in = Rwcs_pixscale(keyvalues=keyvalues_in)
    }
    if(is.null(pixscale_out)){
      pixscale_out = Rwcs_pixscale(keyvalues=keyvalues_out)
    }
  })
  
  if (direction == "auto") {
    if (pixscale_in < pixscale_out) {
      direction = "forward"
    }
    if (pixscale_in >= pixscale_out) {
      direction = "backward"
    }
  }
  
  if(is.null(warpfield)){
    pix_grid = expand.grid(1:dim(image_out$imDat)[1], 1:dim(image_out$imDat)[2])
    
    #if (direction == "forward") {
    #image_out$imDat = imager::imwarp(im = imager::as.cimg(image_out$imDat),
    #                     map = function(x,y){.warpfunc_in2out(x,y)}, direction = direction, coordinates = "absolute",
    #                     boundary = boundary, interpolation = interpolation)
    
    if (direction == "forward") {
      warp_out = .warpfunc_in2out(
        x = pix_grid[, 1],
        y = pix_grid[, 2],
        keyvalues_in = keyvalues_in,
        header_in = header_in,
        WCSref_in = WCSref_in,
        keyvalues_out = keyvalues_out,
        raw_out = raw_out,
        WCSref_out = WCSref_out
      )
    } else if (direction == 'backward') {
      warp_out = .warpfunc_out2in(
        x = pix_grid[, 1],
        y = pix_grid[, 2],
        keyvalues_in = keyvalues_in,
        header_in = header_in,
        WCSref_in = WCSref_in,
        keyvalues_out = keyvalues_out,
        raw_out = raw_out,
        WCSref_out = WCSref_out
      )
    }
    
    warpfield = imager::imappend(list(imager::as.cimg(matrix(
      warp_out[, 1], dim(image_out$imDat)[1], dim(image_out$imDat)[2]
    )),
    imager::as.cimg(matrix(
      warp_out[, 2], dim(image_out$imDat)[1], dim(image_out$imDat)[2]
    ))), 'c')
    
    rm(pix_grid)
    rm(warp_out)
  }
  
  image_out$imDat = imager::warp(
    im = imager::as.cimg(image_out$imDat),
    warpfield = warpfield,
    mode = switch(direction, backward = 0L, forward =
                    2L),
    interpolation = switch(
      interpolation,
      nearest = 0L,
      linear = 1L,
      cubic = 2L
    ),
    boundary_conditions = switch(
      boundary,
      dirichlet = 0L,
      neumann = 1L,
      periodic = 2L
    )
  )
  
  if (dofinenorm) {
    norm = matrix(1, dim(image_out$imDat)[1], dim(image_out$imDat)[2])
    norm = imager::warp(
      im = imager::as.cimg(norm),
      warpfield = warpfield,
      mode = switch(direction, backward = 0L, forward = 2L),
      interpolation = switch(
        interpolation,
        nearest = 0L,
        linear = 1L,
        cubic = 2L
      ),
      boundary_conditions = switch(
        boundary,
        dirichlet = 0L,
        neumann = 1L,
        periodic = 2L
      )
    )
    # norm = imager::imwarp(im = imager::as.cimg(norm),
    #                       map = function(x,y){.warpfunc_in2out(x,y)}, direction = direction, coordinates = "absolute",
    #                       boundary = boundary, interpolation = interpolation)
    image_out$imDat = image_out$imDat / norm
    rm(norm)
  }
  
  if (doscale) {
    image_out$imDat = image_out$imDat * (pixscale_out / pixscale_in) ^ 2
  }
  #}
  # if (direction == "backward") {
  #   # image_out$imDat = imager::imwarp(im = imager::as.cimg(image_out$imDat),
  #   #                      map = function(x,y){.warpfunc_out2in(x,y)}, direction = direction, coordinates = "absolute",
  #   #                      boundary = boundary, interpolation = interpolation)
  #
  #   warp_out = .warpfunc_out2in(pix_grid[,1],
  #                               pix_grid[,2],
  #                               keyvalues_in=keyvalues_in,
  #                               header_in=header_in,
  #                               WCSref_in=WCSref_in,
  #                               keyvalues_out=keyvalues_out,
  #                               raw_out=raw_out,
  #                               WCSref_out=WCSref_out)
  #
  #   warpfield = imager::imappend(list(imager::as.cimg(matrix(warp_out[,1])), imager::as.cimg(matrix(warp_out[,2]))) ,'c')
  #   rm(pix_grid)
  #   rm(warp_out)
  #
  #   image_out$imDat = imager::warp(im = imager::as.cimg(image_out$imDat),
  #                                  warpfield = warpfield,
  #                                  mode = 2L,
  #                                  interpolation = switch(interpolation, nearest=0L, linear=1L, cubic=2L),
  #                                  boundary_conditions = switch(boundary, dirichlet=0L, neumann=1L, periodic=2L)
  #   )
  #
  #   if(doscale){
  #     image_out$imDat = image_out$imDat*(pixscale_out/pixscale_in)^2
  #
  #     if(dofinenorm){
  #       norm = matrix(1, dim(image_out$imDat)[1], dim(image_out$imDat)[2])
  #       norm = imager::imwarp(im = imager::as.cimg(norm),
  #                             map = function(x,y){.warpfunc_in2out(x,y)}, direction = direction, coordinates = "absolute",
  #                             boundary = boundary, interpolation = interpolation)
  #
  #       image_out$imDat = image_out$imDat/norm
  #       rm(norm)
  #     }
  #   }
  # }
  
  image_out$imDat = as.matrix(image_out$imDat)
  
  if(dotightcrop==FALSE | keepcrop==FALSE){
    image_out = image_out[c(1L - (min_x - 1L), dim_out[1] - (min_x - 1L)),c(1L - (min_y - 1L), dim_out[2] - (min_y - 1L)), box=1] #box=1 just in case we have a single pixel left
    image_out$keyvalues$XCUTLO = 1L
    image_out$keyvalues$XCUTHI = dim_out[1]
    image_out$keyvalues$YCUTLO = 1L
    image_out$keyvalues$YCUTHI = dim_out[2]
  }else{
    
    if(max_x > dim_out[1]){
      trim_x = max_x - dim_out[1]
      image_out = image_out[1:(dim(image_out)[1] - trim_x), , box=1] #box=1 just in case we have a single pixel left
      max_x = dim_out[1]
    }
    
    if(max_y > dim_out[2]){
      trim_y = max_y - dim_out[2]
      image_out = image_out[, 1:(dim(image_out)[2] - trim_y), box=1] #box=1 just in case we have a single pixel left
      max_y = dim_out[2]
    }
    
    image_out$keyvalues$XCUTLO = min_x
    image_out$keyvalues$XCUTHI = max_x
    image_out$keyvalues$YCUTLO = min_y
    image_out$keyvalues$YCUTHI = max_y
    
    image_out$keycomments$XCUTLO = 'Low image x range'
    image_out$keycomments$XCUTHI = 'High image x range'
    image_out$keycomments$YCUTLO = 'Low image y range'
    image_out$keycomments$YCUTHI = 'High image y range'
    
    image_out$keynames['XCUTLO'] = 'XCUTLO'
    image_out$keynames['XCUTHI'] = 'XCUTHI'
    image_out$keynames['YCUTLO'] = 'YCUTLO'
    image_out$keynames['YCUTHI'] = 'YCUTHI'
    
    image_out$crop = c(xlo=min_x, xhi=max_x, ylo=min_y, yhi=max_y) #we want to keep the subset location for potential later writing
  }
  
  if (plot) {
    Rwcs_image(image_out, ...)
  }
  
  if(warpfield_return){
    image_out$warpfield = warpfield
  }
  
  return(invisible(image_out))
}

Rwcs_rebin = function(image, scale = 1,interpolation = 6){
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("The imager package is needed for this function to work. Please install it from CRAN.", 
         call. = FALSE)
  }
  
  if(inherits(image,'Rfits_image')){
    #downsample = 2^floor(log(downsample,2))
    #dim_use = downsample * floor(dim(image) / downsample)
    #image = image[1:dim_use[1], 1:dim_use[2]]
    image_resize = as.matrix(imager::resize(imager::as.cimg(image$imDat), -100*scale, -100*scale))
    norm = matrix(1, dim(image$imDat)[1], dim(image$imDat)[2])
    norm_resize = as.matrix(imager::resize(imager::as.cimg(norm), -100*scale, -100*scale))
    image_resize = (image_resize / norm_resize) / scale^2
    
    # if(all(dim(image$imDat)*scale %% 1 == 0)){
    #   image_resize = image_resize * sum(image$imDat, na.rm = TRUE) / sum(image_resize, na.rm = TRUE)
    # }else{
    #   image_resize = image_resize / (scale^2)
    # }
    
    keyvalues_out = image$keyvalues
    keyvalues_out$NAXIS1 = dim(image_resize)[1]
    keyvalues_out$NAXIS2 = dim(image_resize)[2]
    keyvalues_out$CRPIX1 = keyvalues_out$CRPIX1 * scale
    keyvalues_out$CRPIX2 = keyvalues_out$CRPIX2 * scale
    keyvalues_out$CD1_1 = keyvalues_out$CD1_1 / scale
    keyvalues_out$CD1_2 = keyvalues_out$CD1_2 / scale
    keyvalues_out$CD2_1 = keyvalues_out$CD2_1 / scale
    keyvalues_out$CD2_2 = keyvalues_out$CD2_2 / scale
    
    image_out = list(
      imDat = image_resize,
      keyvalues = keyvalues_out,
      hdr = Rfits::Rfits_keyvalues_to_hdr(keyvalues_out),
      header = Rfits::Rfits_keyvalues_to_header(keyvalues_out),
      raw = Rfits::Rfits_header_to_raw(Rfits::Rfits_keyvalues_to_header(keyvalues_out)),
      keynames = names(keyvalues_out),
      keycomments = as.list(rep('', length(keyvalues_out)))
    )
    names(image_out$keycomments) = image_out$keynames
    class(image_out) = c('Rfits_image', class(image_out))
    
    return(image_out)
  }else{
    image_resize = as.matrix(imager::resize(imager::as.cimg(image), -100*scale, -100*scale))
    norm = matrix(1, dim(image)[1], dim(image)[2])
    norm_resize = as.matrix(imager::resize(imager::as.cimg(norm), -100*scale, -100*scale))
    image_resize = (image_resize / norm_resize) / scale^2
    
    # if(all(dim(image$imDat)*scale %% 1 == 0)){
    #   image_resize = image_resize * sum(image, na.rm = TRUE) / sum(image_resize, na.rm = TRUE)
    # }else{
    #   image_resize = image_resize / (scale^2)
    # }
    return(image_resize)
  }
}

Rwcs_warp_promo = function (promo_in, keyvalues_out=NULL, dim_out = NULL, magzero_out = NULL, ...){
  
  if(!is.null(magzero_out)){
    zero_point_scale = 10^(-0.4*(promo_in$image$keyvalues$MAGZERO - magzero_out))
  }else{
    magzero_out = promo_in$image$keyvalues$MAGZERO
    zero_point_scale = 1
  }
  
  if(!is.null(promo_in$image)){
    message('warping image')
    
    image_warp = Rwcs_warp(
      image_in = promo_in$image*zero_point_scale,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = TRUE,
      ...
    )
    
    image_warp$keyvalues$EXTNAME = 'image'
    image_warp$keyvalues$MAGZERO = magzero_out
  }else{
    image_warp = NULL
  }
  
  if(!is.null(promo_in$weight)){
    message('warping weight')
    
    weight_warp = Rwcs_warp(
      image_in = promo_in$weight,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )
    
    weight_warp$keyvalues$EXTNAME = 'weight'
  }else{
    weight_warp = NULL
  }
  
  if(!is.null(promo_in$inVar)){
    message('warping inVar')
    
    inVar_warp = Rwcs_warp(
      image_in = promo_in$inVar/(zero_point_scale^2),
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )*(Rwcs_pixscale(promo_in$inVar$keyvalues)^4 / Rwcs_pixscale(keyvalues_out)^4)
    
    inVar_warp$keyvalues$EXTNAME = 'inVar'
    inVar_warp$keyvalues$MAGZERO = magzero_out
  }else{
    inVar_warp = NULL
  }
  
  if(!is.null(promo_in$exp)){
    message('warping exp')
    
    exp_warp = Rwcs_warp(
      image_in = promo_in$exp,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )
    
    exp_warp$keyvalues$EXTNAME = 'exp'
  }else{
    exp_warp = NULL
  }
  
  if(!is.null(promo_in$cold)){
    message('warping cold')
    
    cold_warp = Rwcs_warp(
      image_in = promo_in$cold*zero_point_scale,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = TRUE,
      ...
    )
    
    cold_warp$keyvalues$EXTNAME = 'cold'
    cold_warp$keyvalues$MAGZERO = magzero_out
  }else{
    cold_warp = NULL
  }
  
  if(!is.null(promo_in$hot)){
    message('warping hot')
    
    hot_warp = Rwcs_warp(
      image_in = promo_in$hot*zero_point_scale,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = TRUE,
      ...
    )
    
    hot_warp$keyvalues$EXTNAME = 'hot'
    hot_warp$keyvalues$MAGZERO = magzero_out
  }else{
    hot_warp = NULL
  }
  
  if(!is.null(promo_in$clip)){
    message('warping clip')
    
    clip_warp = Rwcs_warp(
      image_in = promo_in$clip,
      keyvalues_out = keyvalues_out,
      dim_out = dim_out,
      doscale = FALSE,
      ...
    )
    
    clip_warp$keyvalues$EXTNAME = 'clip'
  }else{
    clip_warp = NULL
  }
  
  output = list(
    image = image_warp,
    weight = weight_warp,
    inVar = inVar_warp,
    exp = exp_warp,
    cold = cold_warp,
    hot = hot_warp,
    clip=clip_warp
  )
  
  class(output) = "ProMo"
  return(invisible(output))
}