Rwcs_tweak = function(image_ref, image_pre_fix, delta_max=3, quan_cut=0.99, Nmeta=3,
                      cores=4, shift_int=TRUE, return_image=FALSE, direction='backward',
                      final_centre=TRUE, verbose=TRUE){
  
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }
  
  timestart = proc.time()[3]
  
  dim_orig = dim(image_ref)
  
  if(verbose){
    message("Trimming comparison region and auto scaling")
  }
  
  trim_lim = which((!is.na(image_ref)) & (!is.na(image_pre_fix)), arr.ind = TRUE)
  x_lo = min(trim_lim[,1])
  x_hi = max(trim_lim[,1])
  y_lo = min(trim_lim[,2])
  y_hi = max(trim_lim[,2])
  
  image_ref = image_ref[x_lo:x_hi, y_lo:y_hi]
  image_ref = image_ref - median(image_ref, na.rm=TRUE)
  image_ref = image_ref / quantile(image_ref, quan_cut, na.rm=TRUE)
  image_ref[image_ref < 1] =NA
  
  image_pre_fix_orig = image_pre_fix[x_lo:x_hi, y_lo:y_hi]
  image_pre_fix = image_pre_fix_orig - median(image_pre_fix_orig, na.rm=TRUE)
  image_pre_fix = image_pre_fix / quantile(image_pre_fix, quan_cut, na.rm=TRUE)
  image_pre_fix[image_pre_fix < 1] = NA
  
  scale = 1 # not really used anymore
  i = NULL #to avoid warnings
  
  #pix_cost_use = which(!is.na(image_ref) & !is.na(image_pre_fix))
  
  if(verbose){
    message('Optimising tweak')
  }
  
  if(shift_int){
    registerDoParallel(cores=cores)
    
    Ndelta = length(-delta_max:delta_max)
    
    for(Nshift in 1:Nmeta){
      
      if(verbose){
        message('  Meta shift ', Nshift)
      }
      
      if(Nshift == 1){
        current_par = c(0L, 0L)
      }
      
      grid_search = expand.grid(-delta_max:delta_max + current_par[1], -delta_max:delta_max + current_par[2])
      
      cost_mat = foreach(i = 1:dim(grid_search)[1], .combine='c')%dopar%{
        cost = .mat_diff_sum(image_ref, image_pre_fix, scale, grid_search[i,1], grid_search[i,2])
        if(verbose){
          message(grid_search[i,1],' ', grid_search[i,2],' ',cost)
        }
        return(cost)
      }
      
      current_par = as.integer(grid_search[which.min(cost_mat),])
      
      if(verbose){
        message('    Current best: ', current_par[1], ' ', current_par[2],': ',round(min(cost_mat)))
      }
      
      at_lim = current_par[1] == min(grid_search[,1]) | current_par[1] == max(grid_search[,1]) | current_par[2] == min(grid_search[,2]) | current_par[2] == max(grid_search[,2])
      if(at_lim == FALSE){
        if(final_centre){
          if(verbose){
            message('  Final centre shift')
          }
          
          grid_search = expand.grid(-delta_max:delta_max + current_par[1], -delta_max:delta_max + current_par[2])
          
          cost_mat = foreach(i = 1:dim(grid_search)[1], .combine='c')%dopar%{
            cost = .mat_diff_sum(image_ref, image_pre_fix, scale, grid_search[i,1], grid_search[i,2])
            if(verbose){
              message(grid_search[i,1],' ', grid_search[i,2],' ',cost)
            }
            return(cost)
          }
          
          current_par = as.integer(grid_search[which.min(cost_mat),])
          
          break
        }else{
          break
        }
      }else{
        if(verbose){
          message('    Limit hit!')
        }
      }
    }
    
    optim_out = list()
    optim_out$par = current_par
    optim_out$value = min(cost_mat)
    optim_out$counts = Ndelta^2
    
    cos_weight = exp(-(cost_mat - min(cost_mat))/2)
    par_weight_x = sum(grid_search[,1] * cos_weight) / sum(cos_weight)
    par_weight_y = sum(grid_search[,2] * cos_weight) / sum(cos_weight)
    optim_out$par_weight = c(par_weight_x, par_weight_y)
    
  }else{
    optim_out = optim(par = c(0,0),
                      fn = .cost_fn,
                      method = "L-BFGS-B",
                      lower = c(-delta_max,-delta_max),
                      upper = c(delta_max,delta_max),
                      image_ref = image_ref,
                      image_pre_fix = image_pre_fix,
                      scale = scale,
                      direction = direction,
                      pix_cost_use = NULL,
                      shift_int = FALSE
    )
  }
  
  if(return_image){
    if(verbose){
      message('Creating output image_post_fix')
    }
    
    image_post_fix_temp = .cost_fn(par = optim_out$par,
                                   image_ref = image_ref,
                                   image_pre_fix = image_pre_fix_orig,
                                   scale = scale,
                                   direction = direction,
                                   return = 'image_post_fix',
                                   shift_int = FALSE)
    
    image_post_fix = matrix(NA, dim_orig[1], dim_orig[2])
    image_post_fix[x_lo:x_hi, y_lo:y_hi] = image_post_fix_temp
  }else{
    image_post_fix = NULL
  }
  
  if(shift_int){
    return(invisible(list(optim_out = optim_out, image_post_fix = image_post_fix, time = proc.time()[3] - timestart, cost_mat=list(x= -delta_max:delta_max + current_par[1], y= -delta_max:delta_max + current_par[2], z=matrix(cost_mat,Ndelta,Ndelta)), at_lim=at_lim)))
  }else{
    return(invisible(list(optim_out = optim_out, image_post_fix = image_post_fix, time = proc.time()[3] - timestart)))
  }
}

Rwcs_tran = function(image, delta_x=0, delta_y=0, direction='backward', padNA=TRUE, shift_int=TRUE){
  #delta refers to the direction we shift the image, not the view point.
  #postive delta_x moves image to the right
  #positive delta_y moves image up
  
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }
  
  if(shift_int){
    if(padNA){
      min_val = min(image, na.rm=TRUE)
      image = image - min_val + 1
      image = as.matrix(imager::imshift(imager::as.cimg(image), delta_x=round(delta_x), delta_y=round(delta_y)))
      image[image == 0] = NA
      image = image + min_val - 1
      return(image)
    }else{
      return(as.matrix(imager::imshift(imager::as.cimg(image), delta_x=round(delta_x), delta_y=round(delta_y))))
    }
  }else{
    if(direction == 'forward'){
      formals(.map.tran)$delta_x = delta_x
      formals(.map.tran)$delta_y = delta_y
    }
    if(direction == 'backward'){
      formals(.map.tran)$delta_x = -delta_x
      formals(.map.tran)$delta_y = -delta_y
    }
    local_fun = function(x, y){
      .map.tran(x = x, y = y)
    }
    norm_mat = matrix(1, dim(image)[1], dim(image)[2])
    norm_mat = as.matrix(imager::imwarp(imager::as.cimg(norm_mat), map=local_fun, direction=direction, coordinates='absolute'))
    if(padNA){
      min_val = min(image, na.rm=TRUE)
      image = image - min_val + 1
      image = as.matrix(imager::imwarp(imager::as.cimg(image), map=local_fun, direction=direction, coordinates='absolute'))
      image[image == 0] = NA
      image = image/norm_mat
      image = image + min_val - 1
      return(image)
    }else{
      return(as.matrix(imager::imwarp(imager::as.cimg(image), map=local_fun, direction=direction, coordinates='absolute'))/norm_mat)
    }
  }
}

.map.tran = function(x=0, y=0, delta_x=0, delta_y=0){
  list(x = x + delta_x, y = y + delta_y)
}

.cost_fn = function(par, image_ref, image_pre_fix, scale=1, direction='backward', pix_cost_use=NULL, shift_int = TRUE, return='cost'){
  if(shift_int){
    cost = .mat_diff_sum(image_ref, image_pre_fix, scale, par[1], par[2])
    message(par[1],' ',par[2],' ',cost)
    return(cost)
  }else{
    image_post_fix = Rwcs_tran(image_pre_fix, delta_x=par[1], delta_y=par[2], padNA=FALSE, shift_int=FALSE)
    if(return=='image_post_fix'){
      return(image_post_fix)
    }
    #frac_good = length(which(!is.na(image_post_fix))) / prod(dim(image_post_fix))
    if(is.null(pix_cost_use)){
      #cost = sum(asinh((image_ref * image_post_fix)/scale/frac_good), na.rm=TRUE)
      cost = sum(((image_ref - image_post_fix)/scale)^2, na.rm=TRUE)
    }else{
      #cost = sum(asinh((image_ref[pix_cost_use] * image_post_fix[pix_cost_use])/scale/frac_good), na.rm=TRUE)
      cost = sum(((image_ref - image_post_fix[pix_cost_use])/scale)^2, na.rm=TRUE)
    }
    #message(par[1],' ',par[2],' ',cost)
    if(return=='cost'){
      return(cost)
    }
    if(return=='all'){
      return(list(cost=cost, image_post_fix=image_post_fix))
    }
  }
}
