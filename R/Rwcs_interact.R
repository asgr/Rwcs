Rwcs_interact = function(mode='zoom', coord.type="sex", col='green', point.pch=4, point.cex=5, ...){
  if(!requireNamespace("Rfits", quietly=TRUE)){
    stop('The Rfits package is needed to interact.')
  }
  
  if(mode == 'coord'){
    temp_loc = locator(type='p', pch=point.pch, col=col, cex=point.cex)
    temp_coord = Rwcs_p2s(temp_loc$x, temp_loc$y, pixcen = 'R')
    if(coord.type == 'sex'){
      temp_coord = cbind(RA=deg2hms(temp_coord[,1], type='cat'), Dec=deg2dms(temp_coord[,2], type='cat'))
    }
    return(temp_coord)
  }
  
  if(mode == 'zoom'){
    image = get(".current_image", envir = .GlobalEnv)
    temp_loc = locator(n=2, type='p', pch=point.pch, col=col, cex=point.cex)
    
    if(is.null(temp_loc)){
      return(invisible(NULL))
    }
    
    temp_coord_user = Rwcs_p2s(temp_loc$x, temp_loc$y, pixcen = 'R') #so we find the desired centre again
    temp_pix_user = Rwcs_s2p(temp_coord_user, keyvalues=image$keyvalues, header=image$raw, pixcen = 'R')
    
    if(length(temp_loc$x) == 1){
      if(coord.type == 'sex'){
        message(deg2hms(temp_coord_user[1,'RA'], type='cat'), ' ', deg2dms(temp_coord_user[1,'Dec'], type='cat'))
      }else if(coord.type == 'deg'){
        message(signif(temp_coord_user[1,'RA'], digits=8), ' ', signif(temp_coord_user[1,'Dec'], digits=8), type='cat')
      }
      
      temp_coord_corner = Rfits::corners(options()$current_keyvalues)
      temp_pix_corner = Rwcs_s2p(temp_coord_corner, keyvalues=image$keyvalues, header=image$raw, pixcen = 'R')
      limx_corner = ceiling(abs(diff(range(temp_pix_corner[,'x']))))
      limy_corner = ceiling(abs(diff(range(temp_pix_corner[,'y']))))
      box = c(limx_corner, limy_corner)
    }else if(length(temp_loc$x) == 2){
      rect(xleft = temp_loc$x[1], ybottom = temp_loc$y[1], xright = temp_loc$x[2], ytop = temp_loc$y[2], border=col)
    
      limx_user = ceiling(abs(diff(range(temp_pix_user[,'x']))))
      limy_user = ceiling(abs(diff(range(temp_pix_user[,'y']))))
      
      if(temp_loc$x[1] <= temp_loc$x[2]){
        #Zooming in
        box = max(limx_user, limy_user)
      }else{
        #Zooming out
        lines(c(temp_loc$x[1], par()$usr[2]), c(temp_loc$y[1], par()$usr[3]), col=col)
        lines(c(temp_loc$x[1], par()$usr[2]), c(temp_loc$y[2], par()$usr[4]), col=col)
        lines(c(temp_loc$x[2], par()$usr[1]), c(temp_loc$y[1], par()$usr[3]), col=col)
        lines(c(temp_loc$x[2], par()$usr[1]), c(temp_loc$y[2], par()$usr[4]), col=col)
        
        temp_coord_corner = Rfits::corners(options()$current_keyvalues)
        temp_pix_corner = Rwcs_s2p(temp_coord_corner, keyvalues=image$keyvalues, header=image$raw, pixcen = 'R')
        limx_corner = ceiling(abs(diff(range(temp_pix_corner[,'x']))))
        limy_corner = ceiling(abs(diff(range(temp_pix_corner[,'y']))))
        max_diff_corner = max(limx_corner, limy_corner)
        
        max_diff_user = max(limx_user, limy_user)
        
        zoom = max_diff_corner/max_diff_user
        box = min(ceiling(max_diff_corner*zoom), max(dim(image)))
      }
    }
    
    sparse = ceiling(max(box)/1e3)
    
    if(inherits(image, 'Rfits_image')){
      image = image[mean(temp_coord_user[,'RA']), mean(temp_coord_user[,'Dec']), type='coord', box=box]
      plot(image, sparse = sparse, interactive=FALSE, ...)
    }else if(inherits(image, 'Rfits_pointer')){
      image = image[mean(temp_coord_user[,'RA']), mean(temp_coord_user[,'Dec']), type='coord', box=box, sparse=sparse]
      plot(image, interactive=FALSE, ...)
    }
    
    Rwcs_interact(mode=mode, coord.type=coord.type, col=col, point.pch=point.pch, point.cex=point.cex, ...)
  }
}
