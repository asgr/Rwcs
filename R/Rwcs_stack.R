Rwcs_stack = function(image_list=NULL, inVar_list=NULL, exp_list=NULL, weight_list=NULL, mask_list=NULL, magzero_in=0,
                      magzero_out=23.9, keyvalues_out=NULL, dim_out=NULL, cores=4, Nbatch=cores,
                      keep_extreme_pix=FALSE, doclip=FALSE, clip_tol=100, clip_dilate=0, clip_sigma=5,
                      return_all=FALSE, dump_frames=FALSE, dump_dir=tempdir(), ...){
  
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }
  
  timestart = proc.time()[3]
  
  if(dump_frames){
    message('Frames being dumped to ', dump_dir)
  }
  
  registerDoParallel(cores=cores)
  
  Nim = length(image_list)
  
  if(Nbatch > Nim){
    Nbatch = Nim
  }
  
  if(cores > Nbatch){
    cores = Nbatch
  }
  
  if(return_all){
    Nbatch = Nim
  }
  
  message('Stacking ', Nim,' images; Nbatch: ', Nbatch,'; cores: ',cores)
  
  if(length(magzero_in) == 1){
    magzero_in = rep(magzero_in, Nim) 
  }
  
  zero_point_scale = 10^(-0.4*(magzero_in - magzero_out))
  
  if(is.null(keyvalues_out)){
    keyvalues_out = image_list[[1]]$keyvalues
  }
  
  if(is.null(keyvalues_out)){
    stop('Need keyvalues out!') 
  }
  
  if(isTRUE(keyvalues_out$ZIMAGE)){
    dim_im = c(keyvalues_out$ZNAXIS1, keyvalues_out$ZNAXIS2)
  }else{
    dim_im = c(keyvalues_out$NAXIS1, keyvalues_out$NAXIS2)
  }
  
  mask_clip = NULL
  
  post_stack_image = matrix(0, dim_im[1], dim_im[2])
  
  if(!is.null(inVar_list)){
    
    if(length(inVar_list) == 1){
      inVar_list = rep(inVar_list, Nim)
    }
    
    if(length(inVar_list) != Nim){
      stop("Length of inVar_list not equal to length of image_list!")  
    }
    
    post_stack_inVar = matrix(0, dim_im[1], dim_im[2])
  }else{
    post_stack_inVar = NULL
  }
  
  if(!is.null(exp_list)){
    
    if(length(exp_list) == 1){
      exp_list = rep(exp_list, Nim) 
    }
    
    if(length(exp_list) != Nim){
      stop("Length of exp_list not equal to length of image_list!")  
    }
    
    post_stack_exp = matrix(0, dim_im[1], dim_im[2])
  }else{
    post_stack_exp = NULL
  }
  
  if(!is.null(weight_list)){ #new weight_list option. This should allow us to easily continue stacking data
    
    if(length(weight_list) == 1){
      weight_list = rep(weight_list, Nim) 
    }
    
    if(length(weight_list) != Nim){
      stop("Length of weight_list not equal to length of image_list!")  
    }
  }else{
    weight_list = rep(1L, Nim) 
  }
  
  weight_image = sapply(weight_list, is.matrix)
  
  post_stack_weight = matrix(0L, dim_im[1], dim_im[2])
  
  if(keep_extreme_pix | doclip){
    post_stack_cold = matrix(Inf, dim_im[1], dim_im[2])
    post_stack_hot = matrix(-Inf, dim_im[1], dim_im[2])
    mask_clip = matrix(0L, dim_im[1], dim_im[2])
  }else{
    post_stack_cold = NULL
    post_stack_hot = NULL
    mask_clip = NULL
  }
  
  if(doclip){
    
    if(length(clip_tol) == 1){
      clip_tol = rep(clip_tol, 2)
    }
    
    if(!is.null(inVar_list)){
      post_stack_cold_id = matrix(0L, dim_im[1], dim_im[2])
      post_stack_hot_id = matrix(0L, dim_im[1], dim_im[2])
    }else{
      stop('inVar_list required if doclip = TRUE!')
    }
  }
  
  # Check all supplied frames are in WCS:
  
  which_overlap = which(foreach(i = 1:length(image_list), .combine='c')%dopar%{
    Rwcs_overlap(image_list[[i]]$keyvalues, keyvalues_ref = keyvalues_out)
  })
  
  Ncheck = length(which_overlap)
  
  if(Ncheck == 0){
    stop('No frames exist within target WCS!')
  }
  
  if(Nim > Ncheck){
    message('Only ', Ncheck, ' of ', Nim, ' input frames overlap with target WCS!')
    
    Nim = Ncheck
    
    # Trim down all the inputs to just those overlapping with the target WCS:
    
    image_list = image_list[which_overlap]
    
    weight_list = weight_list[which_overlap]
    
    magzero_in = magzero_in[which_overlap]
    
    if(!is.null(exp_list)){
      exp_list = exp_list[which_overlap]
    }
    
    if(!is.null(inVar_list)){
      inVar_list = inVar_list[which_overlap]
    }
    
    if(!is.null(mask_list)){
      mask_list = mask_list[which_overlap]
    }
  }
  
  seq_process = seq(1,Nim,by=Nbatch)
  
  for(seq_start in seq_process){
    seq_end = min(seq_start + Nbatch - 1L, Nim)
    
    Nbatch_sub = length(seq_start:seq_end)
    
    pre_stack_image_list = NULL
    pre_stack_inVar_list = NULL
    pre_stack_exp_list = NULL
    pre_stack_weight_list = NULL
    
    message('Projecting Images ',seq_start,' to ',seq_end,' of ',Nim)
    
    pre_stack_image_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp'))%dopar%{
      if(inherits(image_list[[i]], 'Rfits_pointer')){
        temp_image = image_list[[i]][,]
      }else{
        temp_image = image_list[[i]]
      }
      
      if(zero_point_scale[i] != 1){
        temp_image$imDat = temp_image$imDat*zero_point_scale[i]
      }
      
      if(any(!is.finite(temp_image$imDat))){
        temp_image$imDat[!is.finite(temp_image$imDat)] = NA
      }
      
      if(!is.null(mask_list)){
        if(inherits(mask_list[[i]], 'Rfits_pointer')){
          mask_list[[i]]$header = FALSE
          temp_image$imDat[mask_list[[i]][,] != 0] = NA
        }else if(inherits(mask_list[[i]], 'Rfits_image')){
          temp_image$imDat[mask_list[[i]]$imDat != 0] = NA
        }else{
          temp_image$imDat[mask_list != 0] = NA
        }
      }
      
      suppressMessages({
        temp_warp = Rwcs_warp(
          image_in = temp_image,
          keyvalues_out = keyvalues_out,
          dim_out = dim_out,
          doscale = TRUE,
          dotightcrop = TRUE,
          keepcrop = TRUE,
          warpfield_return = TRUE,
          ...
        )
        
        if(dump_frames){
          Rfits::Rfits_write_image(temp_warp, paste0(dump_dir,'/image_warp_',i,'.fits'))
        }
      })
      return(temp_warp)
    }
    
    if(!is.null(inVar_list)){
      message('Projecting Inverse Variance ',seq_start,' to ',seq_end,' of ',Nim)
      
      pre_stack_inVar_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp'))%dopar%{
        if(inherits(image_list[[i]], 'Rfits_pointer')){
          temp_inVar = image_list[[i]][,]
        }else{
          temp_inVar = image_list[[i]]
        }
        
        if(inherits(inVar_list[[i]], 'Rfits_pointer')){
          inVar_list[[i]]$header = FALSE
          temp_inVar$imDat[] = inVar_list[[i]][,]
        }else if(inherits(inVar_list[[i]], 'Rfits_image')){
          temp_inVar$imDat[] = inVar_list[[i]]$imDat
        }else{
          temp_inVar$imDat[] = inVar_list[[i]]
        }
        
        if(zero_point_scale[i] != 1){
          temp_inVar$imDat = temp_inVar$imDat/(zero_point_scale[i]^2)
        }
        
        suppressMessages({
          temp_warp = Rwcs_warp(
            image_in = temp_inVar,
            keyvalues_out = keyvalues_out,
            dim_out = dim_out,
            doscale = FALSE,
            dotightcrop = TRUE,
            keepcrop = TRUE,
            warpfield = pre_stack_image_list[[i - seq_start + 1L]]$warpfield,
            ...
            )*(Rwcs_pixscale(temp_inVar$keyvalues)^4 / Rwcs_pixscale(keyvalues_out)^4) #this is because RMS scales as linear pixel area. Using Rfits * method here

          if(dump_frames){
            Rfits::Rfits_write_image(temp_warp, paste0(dump_dir,'/inVar_warp_',i,'.fits'))
          }
        })
        return(temp_warp)
      }
    }
    
    if(!is.null(exp_list)){
      message('Projecting Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
      if(length(exp_list) != Nim){
        stop("Length of Exposure Times not equal to length of image_list!")
      }
      pre_stack_exp_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp', 'pre_stack_inVar_list'))%dopar%{
        if(inherits(image_list[[i]], 'Rfits_pointer')){
          temp_exp = image_list[[i]][,]
        }else{
          temp_exp = image_list[[i]]
        }
        
        #need [] because this will assign a single value to all elements of a matrix
        if(inherits(exp_list[[i]], 'Rfits_pointer')){
          exp_list[[i]]$header = FALSE
          temp_exp$imDat[] = exp_list[[i]][,]
        }else if(inherits(exp_list[[i]], 'Rfits_image')){
          temp_exp$imDat[] = exp_list[[i]]$imDat
        }else{
          temp_exp$imDat[] = exp_list[[i]]
        }
        
        suppressMessages({
          temp_warp = Rwcs_warp(
            image_in = temp_exp,
            keyvalues_out = keyvalues_out,
            dim_out = dim_out,
            doscale = FALSE,
            dotightcrop = TRUE,
            keepcrop = TRUE,
            warpfield = pre_stack_image_list[[i - seq_start + 1L]]$warpfield,
            ...
          )
          
          if(dump_frames){
            Rfits::Rfits_write_image(temp_warp, paste0(dump_dir,'/exp_warp_',i,'.fits'))
          }
        })
        return(temp_warp)
      }
    }
    
    #new weight projections (if relevant)
    
    if(any(weight_image)){
      message('Projecting Weights ',seq_start,' to ',seq_end,' of ',Nim)
      pre_stack_weight_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp', 'pre_stack_inVar_list', 'pre_stack_exp_list'))%dopar%{
        if(weight_image[i]){
          if(inherits(image_list[[i]], 'Rfits_pointer')){
            temp_weight = image_list[[i]][,]
          }else{
            temp_weight = image_list[[i]]
          }
          
          #need [] because this will assign a single value to all elements of a matrix
          if(inherits(weight_list[[i]], 'Rfits_pointer')){
            weight_list[[i]]$header = FALSE
            temp_weight$imDat[] = weight_list[[i]][,]
          }else if(inherits(weight_list[[i]], 'Rfits_image')){
            temp_weight$imDat[] = weight_list[[i]]$imDat
          }else{
            temp_weight$imDat[] = weight_list[[i]]
          }
          
          suppressMessages({
            temp_warp = Rwcs_warp(
              image_in = temp_weight,
              keyvalues_out = keyvalues_out,
              dim_out = dim_out,
              doscale = FALSE,
              dotightcrop = TRUE,
              keepcrop = TRUE,
              warpfield = pre_stack_image_list[[i - seq_start + 1L]]$warpfield,
              ...
            )
            
            if(dump_frames){
              Rfits::Rfits_write_image(temp_warp, paste0(dump_dir,'/weight_warp_',i,'.fits'))
            }
          })
          return(temp_warp)
        }else{
          return(weight_list[[i]])
        }
      }
    }else{
      pre_stack_weight_list = weight_list
    }
    
    if(is.null(pre_stack_inVar_list)){
      message('Stacking Images ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:Nbatch_sub){
        if(weight_image[i]){
          pre_weight = pre_stack_weight_list[[i]]$imDat
        }else{
          pre_weight = weight_list[[i]]
        }
        .stack_image_cpp(post_image = post_stack_image,
                         post_weight = post_stack_weight,
                         pre_image = pre_stack_image_list[[i]]$imDat,
                         pre_weight_sexp = pre_weight,
                         offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')])
        )
      }
    }else{
      message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:Nbatch_sub){
        if(weight_image[i]){
          pre_weight = pre_stack_weight_list[[i]]$imDat
        }else{
          pre_weight = weight_list[[i]]
        }
        .stack_image_inVar_cpp(post_image = post_stack_image,
                               post_inVar = post_stack_inVar,
                               post_weight = post_stack_weight,
                               pre_image = pre_stack_image_list[[i]]$imDat,
                               pre_inVar = pre_stack_inVar_list[[i]]$imDat,
                               pre_weight_sexp = pre_weight,
                               offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')])
        )
      }
    }
    
    if(!is.null(pre_stack_exp_list)){
      message('Stacking Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:Nbatch_sub){
        .stack_exp_cpp(post_exp = post_stack_exp,
                       pre_exp = pre_stack_exp_list[[i]]$imDat,
                       offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')])
        )
        # if(anyNA(pre_stack_exp_list[[i]]$imDat)){
        #   addID = which(!is.na(pre_stack_exp_list[[i]]$imDat), arr.ind=TRUE)
        #   #addID_sub = addID
        #   #addID_sub[,1] = addID_sub[,1] + pre_stack_image_list[[i]]$crop['xlo'] - 1L
        #   #addID_sub[,2] = addID_sub[,2] + pre_stack_image_list[[i]]$crop['ylo'] - 1L
        #   
        #   #post_stack_exp[addID_sub] = post_stack_exp[addID_sub] + pre_stack_exp_list[[i]]$imDat[addID]
        #   
        #   .num_mat_add_cpp(post_stack_exp, pre_stack_exp_list[[i]]$imDat, addID, pre_stack_image_list[[i]]$crop[c('xlo','ylo')])
        # }else{
        #   xsub = pre_stack_image_list[[i]]$crop['xlo']:pre_stack_image_list[[i]]$crop['xhi']
        #   ysub = pre_stack_image_list[[i]]$crop['ylo']:pre_stack_image_list[[i]]$crop['yhi']
        #   
        #   post_stack_exp[xsub,ysub] = post_stack_exp[xsub,ysub] + pre_stack_exp_list[[i]]$imDat
        # }
      }
    }
    
    if(keep_extreme_pix | doclip){
      message('Calculating Extreme Pixels ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:Nbatch_sub){
        xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
        ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
        
        new_cold = which(pre_stack_image_list[[i]]$imDat < post_stack_cold[xsub,ysub])
        new_hot = which(pre_stack_image_list[[i]]$imDat > post_stack_hot[xsub,ysub])
        
        post_stack_cold[xsub,ysub][new_cold] = pre_stack_image_list[[i]]$imDat[new_cold]
        post_stack_hot[xsub,ysub][new_hot] = pre_stack_image_list[[i]]$imDat[new_hot]
        
        if(doclip){
          post_stack_cold_id[xsub,ysub][new_cold] = seq_start + i - 1L
          post_stack_hot_id[xsub,ysub][new_hot] = seq_start + i - 1L
        }
      }
    }
  }
  
  if(return_all==FALSE){
    pre_stack_exp_list = NULL
  }
  
  #The below is only done once for the stack, so no need to Rcpp it
  
  weight_sel = (post_stack_weight != 0L)
  if(is.null(pre_stack_inVar_list)){
    post_stack_image[weight_sel] = post_stack_image[weight_sel]/post_stack_weight[weight_sel]
  }else{
    post_stack_image[weight_sel] = post_stack_image[weight_sel]/post_stack_inVar[weight_sel]
    post_stack_inVar[!weight_sel] = NA
  }
  
  if(keep_extreme_pix | doclip){
    post_stack_cold[!weight_sel] = NA
    post_stack_hot[!weight_sel] = NA
  }
  
  post_stack_image[!weight_sel] = NA
  
  #Changed my mind on this- I think I should count exposure as a photon hitting a legal part of a sensor (bad pixel or not). I.e. fully masked regions with NA in the final image might still have a positive exposure time.
  # if(!is.null(post_stack_exp)){
  #   post_stack_exp[!weight_sel] = NA
  # }
  
  if(doclip & !is.null(post_stack_inVar)){
    message('Clipping out extreme cold/hot pixels')
    
    post_stack_inRMS = sqrt(post_stack_inVar)
    
    bad_cold = (post_stack_image - post_stack_cold)*post_stack_inRMS > clip_tol[1] & post_stack_weight > 2
    bad_hot = post_stack_cold*post_stack_inRMS < clip_sigma & (post_stack_hot - post_stack_image)*post_stack_inRMS > clip_tol[2] & post_stack_weight > 2
    
    post_stack_cold_id[!bad_cold] = 0L
    post_stack_hot_id[!bad_hot] = 0L
    
    rm(bad_cold)
    rm(bad_hot)
    rm(post_stack_inRMS)
    
    # post_mask_list = list()
    # 
    # for(i in 1:Nim){
    #   post_mask_list[[i]] = (post_stack_cold_id == i) | (post_stack_hot_id == i)
    #   if(clip_dilate > 0){
    #     post_mask_list[[i]] = .dilate_R(post_mask_list[[i]], size=clip_dilate)
    #   }
    # }
    
    #reset post stack outputs
    
    post_stack_image = matrix(0, dim_im[1], dim_im[2])
    post_stack_weight = matrix(0L, dim_im[1], dim_im[2])
    
    if(!is.null(inVar_list)){
      post_stack_inVar = matrix(0, dim_im[1], dim_im[2])
    }else{
      post_stack_inVar = NULL
    }
    
    post_stack_cold = NULL
    post_stack_hot = NULL
    mask_clip = NULL
    
    if(keep_extreme_pix){
      post_stack_cold = matrix(Inf, dim_im[1], dim_im[2])
      post_stack_hot = matrix(-Inf, dim_im[1], dim_im[2])
      if(doclip){
        mask_clip = matrix(0L, dim_im[1], dim_im[2])
      }
    }
    
    if(Nbatch == Nim){ #if we already have all projections in memory we can just re-stack
      message('Restacking without clipped cold/hot pixels')
      if(is.null(pre_stack_inVar_list)){
        message('Stacking Images ',seq_start,' to ',seq_end,' of ',Nim)
        for(i in 1:Nim){
          temp_mask_clip = (post_stack_cold_id == i) | (post_stack_hot_id == i)
          if(clip_dilate > 0){
            temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
          }
          mode(temp_mask_clip) = 'logical'
          
          if(weight_image[i]){
            pre_weight = pre_stack_weight_list[[i]]$imDat
          }else{
            pre_weight = weight_list[[i]]
          }
          .stack_image_cpp(post_image = post_stack_image,
                           post_weight = post_stack_weight,
                           pre_image = pre_stack_image_list[[i]]$imDat,
                           pre_weight_sexp = pre_weight,
                           offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')]),
                           post_mask = temp_mask_clip
          )
          
          if(keep_extreme_pix){
            mask_clip = mask_clip + temp_mask_clip
            
            xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
            ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
            
            new_cold = which(pre_stack_image_list[[i]]$imDat < post_stack_cold[xsub,ysub] & temp_mask_clip[xsub,ysub]==FALSE)
            new_hot = which(pre_stack_image_list[[i]]$imDat > post_stack_hot[xsub,ysub] & temp_mask_clip[xsub,ysub]==FALSE)
            
            post_stack_cold[xsub,ysub][new_cold] = pre_stack_image_list[[i]]$imDat[new_cold]
            post_stack_hot[xsub,ysub][new_hot] = pre_stack_image_list[[i]]$imDat[new_hot]
          }
        }
      }else{
        message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
        for(i in 1:Nim){
          temp_mask_clip = (post_stack_cold_id == i) | (post_stack_hot_id == i)
          if(clip_dilate > 0){
            temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
          }
          mode(temp_mask_clip) = 'logical'
          
          if(weight_image[i]){
            pre_weight = pre_stack_weight_list[[i]]$imDat
          }else{
            pre_weight = weight_list[[i]]
          }
          .stack_image_inVar_cpp(post_image = post_stack_image,
                                 post_inVar = post_stack_inVar,
                                 post_weight = post_stack_weight,
                                 pre_image = pre_stack_image_list[[i]]$imDat,
                                 pre_inVar = pre_stack_inVar_list[[i]]$imDat,
                                 pre_weight_sexp = pre_weight,
                                 offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')]),
                                 post_mask = temp_mask_clip
          )
          
          if(keep_extreme_pix){
            mask_clip = mask_clip + temp_mask_clip
            
            xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
            ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
            
            new_cold = which(pre_stack_image_list[[i]]$imDat < post_stack_cold[xsub,ysub] & temp_mask_clip[xsub,ysub]==FALSE)
            new_hot = which(pre_stack_image_list[[i]]$imDat > post_stack_hot[xsub,ysub] & temp_mask_clip[xsub,ysub]==FALSE)
            
            post_stack_cold[xsub,ysub][new_cold] = pre_stack_image_list[[i]]$imDat[new_cold]
            post_stack_hot[xsub,ysub][new_hot] = pre_stack_image_list[[i]]$imDat[new_hot]
          }
        }
      }
    }else{
      # If we need to batch process the image_list then we need to re-project everything again
      #would maybe be neater to break this out a distinct function, but that is a job for another day...
      
      #I don't think we need to re-zero this actually?
      # if(!is.null(exp_list)){
      #   
      #   if(length(exp_list) == 1){
      #     exp_list = rep(exp_list, Nim) 
      #   }
      #   
      #   if(length(exp_list) != Nim){
      #     stop("Length of exp_list not equal to length of image_list!")  
      #   }
      #   
      #   post_stack_exp = matrix(0, dim_im[1], dim_im[2])
      # }else{
      #   post_stack_exp = NULL
      # }
      
      message('Reprojecting and restacking without clipped cold/hot pixels')
      for(seq_start in seq_process){
        seq_end = min(seq_start + Nbatch - 1L, Nim)
        
        Nbatch_sub = length(seq_start:seq_end)
        
        pre_stack_image_list = NULL
        pre_stack_inVar_list = NULL
        pre_stack_exp_list = NULL
        pre_stack_weight_list = NULL
        
        if(dump_frames){
          message('Loading Images ',seq_start,' to ',seq_end,' of ',Nim)
          
          pre_stack_image_list = foreach(i = seq_start:seq_end)%dopar%{
            Rfits::Rfits_read_image(paste0(dump_dir,'/image_warp_',i,'.fits'))
          }
        }else{
          message('Projecting Images ',seq_start,' to ',seq_end,' of ',Nim)
          
          pre_stack_image_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp'))%dopar%{
            if(inherits(image_list[[i]], 'Rfits_pointer')){
              temp_image = image_list[[i]][,]
            }else{
              temp_image = image_list[[i]]
            }
            
            if(zero_point_scale[i] != 1){
              temp_image$imDat = temp_image$imDat*zero_point_scale[i]
            }
            
            if(any(!is.finite(temp_image$imDat))){
              temp_image$imDat[!is.finite(temp_image$imDat)] = NA
            }
            
            if(!is.null(mask_list)){
              if(inherits(mask_list[[i]], 'Rfits_pointer')){
                mask_list[[i]]$header = FALSE
                temp_image$imDat[mask_list[[i]][,] != 0] = NA
              }else if(inherits(mask_list[[i]], 'Rfits_image')){
                temp_image$imDat[mask_list[[i]]$imDat != 0] = NA
              }else{
                temp_image$imDat[mask_list != 0] = NA
              }
            }
            
            return(Rwcs_warp(
              image_in = temp_image,
              keyvalues_out = keyvalues_out,
              dim_out = dim_out,
              doscale = TRUE,
              dotightcrop = TRUE,
              keepcrop = TRUE,
              warpfield_return = TRUE,
              ...
            ))
          }
        }
        
        if(!is.null(inVar_list)){
          if(dump_frames){
            message('Loading Inverse Variance ',seq_start,' to ',seq_end,' of ',Nim)
            
            pre_stack_inVar_list = foreach(i = seq_start:seq_end)%dopar%{
              Rfits::Rfits_read_image(paste0(dump_dir,'/inVar_warp_',i,'.fits'))
            }
          }else{
            message('Projecting Inverse Variance ',seq_start,' to ',seq_end,' of ',Nim)
            
            pre_stack_inVar_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp'))%dopar%{
              if(inherits(image_list[[i]], 'Rfits_pointer')){
                temp_inVar = image_list[[i]][,]
              }else{
                temp_inVar = image_list[[i]]
              }
              
              if(inherits(inVar_list[[i]], 'Rfits_pointer')){
                inVar_list[[i]]$header = FALSE
                temp_inVar$imDat[] = inVar_list[[i]][,]
              }else if(inherits(inVar_list[[i]], 'Rfits_image')){
                temp_inVar$imDat[] = inVar_list[[i]]$imDat
              }else{
                temp_inVar$imDat[] = inVar_list[[i]]
              }
              
              if(zero_point_scale[i] != 1){
                temp_inVar$imDat = temp_inVar$imDat/(zero_point_scale[i]^2)
              }
              
              suppressMessages({
                return(Rwcs_warp(
                  image_in = temp_inVar,
                  keyvalues_out = keyvalues_out,
                  dim_out = dim_out,
                  doscale = FALSE,
                  dotightcrop = TRUE,
                  keepcrop = TRUE,
                  warpfield = pre_stack_image_list[[i - seq_start + 1L]]$warpfield,
                  ...
                  )*(Rwcs_pixscale(temp_inVar$keyvalues)^4 / Rwcs_pixscale(keyvalues_out)^4) #this is because RMS scales as linear pixel area
                )
              })
            }
          }
        }
        
        # Shouldn't need to reproject exposure times!
        # if(!is.null(exp_list)){
        #   message('Projecting Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
        #   if(length(exp_list) == 1){
        #     exp_list = rep(exp_list, Nim)
        #   }
        #   if(length(exp_list) != Nim){
        #     stop("Length of Exposure Times not equal to length of image_list!")
        #   }
        #   pre_stack_exp_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp', 'pre_stack_image_list', 'pre_stack_inVar_list'))%dopar%{
        #     if(inherits(image_list[[i]], 'Rfits_pointer')){
        #       temp_exp = image_list[[i]][,]
        #     }else{
        #       temp_exp = image_list[[i]]
        #     }
        #     
        #     #need [] because this will assign a single value to all elements of a matrix
        #     if(inherits(exp_list[[i]], 'Rfits_pointer')){
        #       exp_list[[i]]$header = FALSE
        #       temp_exp$imDat[] = exp_list[[i]][,]
        #     }else if(inherits(exp_list[[i]], 'Rfits_image')){
        #       temp_exp$imDat[] = exp_list[[i]]$imDat
        #     }else{
        #       temp_exp$imDat[] = exp_list[[i]]
        #     }
        #     
        #     return(Rwcs_warp(
        #       image_in = temp_exp,
        #       keyvalues_out = keyvalues_out,
        #       dim_out = dim_out,
        #       doscale = FALSE,
        #       dotightcrop = TRUE,
        #       keepcrop = TRUE,
        #       ...
        #     ))
        #   }
        # }else{
        #   pre_stack_exp_list = NULL
        # }
        
        if(any(weight_image)){
          if(dump_frames){
            message('Loading Weights ',seq_start,' to ',seq_end,' of ',Nim)
          }else{
            message('Projecting Weights ',seq_start,' to ',seq_end,' of ',Nim)
          }
          pre_stack_weight_list = foreach(i = seq_start:seq_end,
            .packages = c('Rwcs', 'Rfits'),
            export = c('image_list', 'dump_dir', 'weight_image')
            #.noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'post_stack_exp', 'pre_stack_inVar_list', 'pre_stack_exp_list')
            )%dopar%{
            if(weight_image[i]){
              if(dump_frames){
                Rfits::Rfits_read_image(paste0(dump_dir,'/weight_warp_',i,'.fits'))
              }else{
                if(inherits(image_list[[i]], 'Rfits_pointer')){
                  temp_weight = image_list[[i]][,]
                }else{
                  temp_weight = image_list[[i]]
                }
                
                #need [] because this will assign a single value to all elements of a matrix
                if(inherits(weight_list[[i]], 'Rfits_pointer')){
                  weight_list[[i]]$header = FALSE
                  temp_weight$imDat[] = weight_list[[i]][,]
                }else if(inherits(weight_list[[i]], 'Rfits_image')){
                  temp_weight$imDat[] = weight_list[[i]]$imDat
                }else{
                  temp_weight$imDat[] = weight_list[[i]]
                }
                
                return(Rwcs_warp(
                  image_in = temp_weight,
                  keyvalues_out = keyvalues_out,
                  dim_out = dim_out,
                  doscale = FALSE,
                  dotightcrop = TRUE,
                  keepcrop = TRUE,
                  warpfield = pre_stack_image_list[[i - seq_start + 1L]]$warpfield,
                  ...
                ))
              }
            }else{
              return(weight_list[[i]])
            }
          }
        }else{
          pre_stack_weight_list = weight_list
        }
        
        if(is.null(pre_stack_inVar_list)){
          message('Stacking Images ',seq_start,' to ',seq_end,' of ',Nim)
          for(i in 1:Nbatch_sub){
            i_stack = seq_start + i - 1L
            
            temp_mask_clip = (post_stack_cold_id == i_stack) | (post_stack_hot_id == i_stack)
            if(clip_dilate > 0){
              temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
            }
            mode(temp_mask_clip) = 'logical'
            
            if(weight_image[i]){
              pre_weight = pre_stack_weight_list[[i]]$imDat
            }else{
              pre_weight = weight_list[[i]]
            }
            .stack_image_cpp(post_image = post_stack_image,
                             post_weight = post_stack_weight,
                             pre_image = pre_stack_image_list[[i]]$imDat,
                             pre_weight_sexp = pre_weight,
                             offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')]),
                             post_mask = temp_mask_clip
            )
          }
        }else{
          message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
          for(i in 1:Nbatch_sub){
            i_stack = seq_start + i - 1L
            
            temp_mask_clip = (post_stack_cold_id == i_stack) | (post_stack_hot_id == i_stack)
            if(clip_dilate > 0){
              temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
            }
            mode(temp_mask_clip) = 'logical'
            
            if(weight_image[i]){
              pre_weight = pre_stack_weight_list[[i]]$imDat
            }else{
              pre_weight = weight_list[[i]]
            }
            .stack_image_inVar_cpp(post_image = post_stack_image,
                                   post_inVar = post_stack_inVar,
                                   post_weight = post_stack_weight,
                                   pre_image = pre_stack_image_list[[i]]$imDat,
                                   pre_inVar = pre_stack_inVar_list[[i]]$imDat,
                                   pre_weight_sexp = pre_weight,
                                   offset = unlist(pre_stack_image_list[[i]]$keyvalues[c('XCUTLO','YCUTLO')]),
                                   post_mask = temp_mask_clip
            )
            
            # xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
            # ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
            # 
            # if(anyNA(pre_stack_image_list[[i]]$imDat) | checkmate::anyInfinite(pre_stack_inVar_list[[i]]$imDat) | any(pre_stack_inVar_list[[i]]$imDat < 0, na.rm=TRUE) | any(temp_mask_clip[xsub,ysub])){
            #   addID = which(!is.na(pre_stack_image_list[[i]]$imDat) & is.finite(pre_stack_inVar_list[[i]]$imDat) & pre_stack_inVar_list[[i]]$imDat > 0 & temp_mask_clip[xsub,ysub]==FALSE, arr.ind=TRUE)
            #   #addID_sub = addID
            #   #addID_sub[,1] = addID_sub[,1] + pre_stack_image_list[[i]]$crop['xlo'] - 1L
            #   #addID_sub[,2] = addID_sub[,2] + pre_stack_image_list[[i]]$crop['ylo'] - 1L
            #   
            #   #post_stack_image[addID_sub] = post_stack_image[addID_sub] + pre_stack_image_list[[i]]$imDat[addID]*pre_stack_inVar_list[[i]]$imDat[addID]
            #   #post_stack_inVar[addID_sub] = post_stack_inVar[addID_sub] + pre_stack_inVar_list[[i]]$imDat[addID]
            #   
            #   .num_mat_add_mult_cpp(post_stack_image, pre_stack_image_list[[i]]$imDat, pre_stack_inVar_list[[i]]$imDat, addID, pre_stack_image_list[[i]]$crop[c('xlo','ylo')])
            #   .num_mat_add_cpp(post_stack_inVar, pre_stack_inVar_list[[i]]$imDat, addID, pre_stack_image_list[[i]]$crop[c('xlo','ylo')])
            #   
            #   if(weight_image[i]){
            #     #post_stack_weight[addID_sub] = post_stack_weight[addID_sub] + pre_stack_weight_list[[i]]$imDat[addID]
            #     .int_mat_add_cpp(post_stack_weight, pre_stack_weight_list[[i]]$imDat, addID, pre_stack_image_list[[i]]$crop[c('xlo','ylo')])
            #   }else{
            #     #post_stack_weight[addID_sub] = post_stack_weight[addID_sub] + weight_list[[i]]
            #     .int_mat_add_sin_cpp(post_stack_weight, weight_list[[i]], addID, pre_stack_image_list[[i]]$crop[c('xlo','ylo')])
            #   }
            # }else{
            #   xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
            #   ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
            #   
            #   post_stack_image[xsub,ysub] = post_stack_image[xsub,ysub] + pre_stack_image_list[[i]]$imDat*pre_stack_inVar_list[[i]]$imDat
            #   post_stack_inVar[xsub,ysub] = post_stack_inVar[xsub,ysub] + pre_stack_inVar_list[[i]]$imDat
            #   if(weight_image[i]){
            #     post_stack_weight[xsub,ysub] = post_stack_weight[xsub,ysub] + pre_stack_weight_list[[i]]$imDat
            #   }else{
            #     post_stack_weight[xsub,ysub] = post_stack_weight[xsub,ysub] + weight_list[[i]]
            #   }
            # }
          }
        }
        
        # this should already exist? Probably need to check this more carefully
        # if(!is.null(pre_stack_exp_list)){
        #   message('Stacking Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
        #   for(i in 1:Nbatch_sub){
        #     if(anyNA(pre_stack_exp_list[[i]]$imDat)){
        #       addID = which(!is.na(pre_stack_exp_list[[i]]$imDat), arr.ind=TRUE)
        #       #addID_sub = addID
        #       #addID_sub[,1] = addID_sub[,1] + pre_stack_image_list[[i]]$crop['xlo'] - 1L
        #       #addID_sub[,2] = addID_sub[,2] + pre_stack_image_list[[i]]$crop['ylo'] - 1L
        #       
        #       #post_stack_exp[addID_sub] = post_stack_exp[addID_sub] + pre_stack_exp_list[[i]]$imDat[addID]
        #       
        #       .num_mat_add_cpp(post_stack_exp, pre_stack_exp_list[[i]]$imDat, addID, pre_stack_image_list[[i]]$crop[c('xlo','ylo')])
        #     }else{
        #       xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
        #       ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
        #       
        #       post_stack_exp[xsub,ysub] = post_stack_exp[xsub,ysub] + pre_stack_exp_list[[i]]$imDat
        #     }
        #   }
        # }
        
        if(keep_extreme_pix){
          message('Calculating Extreme Pixels ',seq_start,' to ',seq_end,' of ',Nim)
          for(i in 1:Nbatch_sub){
            i_stack = seq_start + i - 1L
            
            temp_mask_clip = (post_stack_cold_id == i_stack) | (post_stack_hot_id == i_stack)
            if(clip_dilate > 0){
              temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
            }
            mode(temp_mask_clip) = 'logical'

            mask_clip = mask_clip + temp_mask_clip
            
            xsub = pre_stack_image_list[[i]]$keyvalues$XCUTLO:pre_stack_image_list[[i]]$keyvalues$XCUTHI
            ysub = pre_stack_image_list[[i]]$keyvalues$YCUTLO:pre_stack_image_list[[i]]$keyvalues$YCUTHI
            
            new_cold = which(pre_stack_image_list[[i]]$imDat < post_stack_cold[xsub,ysub] & temp_mask_clip[xsub,ysub]==FALSE)
            new_hot = which(pre_stack_image_list[[i]]$imDat > post_stack_hot[xsub,ysub] & temp_mask_clip[xsub,ysub]==FALSE)
            
            post_stack_cold[xsub,ysub][new_cold] = pre_stack_image_list[[i]]$imDat[new_cold]
            post_stack_hot[xsub,ysub][new_hot] = pre_stack_image_list[[i]]$imDat[new_hot]
          }
        }
      }
    }
    
    if(return_all==FALSE){
      pre_stack_exp_list = NULL
    }
    
    weight_sel = (post_stack_weight != 0L)
    
    if(is.null(pre_stack_inVar_list)){
      post_stack_image[weight_sel] = post_stack_image[weight_sel]/post_stack_weight[weight_sel]
    }else{
      post_stack_image[weight_sel] = post_stack_image[weight_sel]/post_stack_inVar[weight_sel]
      post_stack_inVar[!weight_sel] = NA
    }
    
    if(keep_extreme_pix){
      post_stack_cold[!weight_sel] = NA
      post_stack_hot[!weight_sel] = NA
    }
    
    post_stack_image[!weight_sel] = NA
  }
  
  if(return_all==FALSE){
    pre_stack_image_list = NULL
    pre_stack_inVar_list = NULL
    pre_stack_exp_list = NULL
    pre_stack_weight_list = NULL
  }
  
  keyvalues_out$EXTNAME = 'image'
  keyvalues_out$MAGZERO = magzero_out
  keyvalues_out$R_VER = R.version$version.string
  keyvalues_out$RWCS_VER = as.character(packageVersion('Rwcs'))
  
  image_out = Rfits::Rfits_create_image(image=post_stack_image,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE,
                                 history='Stacked with Rwcs_stack')
  
  keyvalues_out$EXTNAME = 'weight'
  keyvalues_out$MAGZERO = NULL
  weight_out = Rfits::Rfits_create_image(image=post_stack_weight,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE)
  
  if(!is.null(post_stack_inVar)){
    keyvalues_out$EXTNAME = 'inVar'
    keyvalues_out$MAGZERO = magzero_out
    inVar_out = Rfits::Rfits_create_image(image=post_stack_inVar,
                                    keyvalues=keyvalues_out,
                                    keypass=FALSE)
  }else{
    inVar_out = NULL
  }
  
  if(!is.null(post_stack_exp)){
    keyvalues_out$EXTNAME = 'exp'
    keyvalues_out$MAGZERO = NULL
    exp_out = Rfits::Rfits_create_image(image=post_stack_exp,
                                          keyvalues=keyvalues_out,
                                          keypass=FALSE)
  }else{
    exp_out = NULL
  }
  
  if(keep_extreme_pix){
    post_stack_cold[weight_out$imDat == 0L] = NA
    
    keyvalues_out$EXTNAME = 'cold'
    keyvalues_out$MAGZERO = magzero_out
    cold_out = Rfits::Rfits_create_image(image=post_stack_cold,
                                        keyvalues=keyvalues_out,
                                        keypass=FALSE)
    
    post_stack_hot[weight_out$imDat == 0L] = NA
    
    keyvalues_out$EXTNAME = 'hot'
    hot_out = Rfits::Rfits_create_image(image=post_stack_hot,
                                         keyvalues=keyvalues_out,
                                         keypass=FALSE)
    
    if(!is.null(mask_clip)){
      mask_clip[weight_out$imDat == 0L] = NA
      
      keyvalues_out$EXTNAME = 'clip'
      keyvalues_out$MAGZERO = NULL
      mask_clip = Rfits::Rfits_create_image(image=mask_clip,
                                          keyvalues=keyvalues_out,
                                          keypass=FALSE)
    }
  }else{
    cold_out = NULL
    hot_out = NULL
    mask_clip = NULL
  }
  
  time_taken = proc.time()[3] - timestart
  message('Time taken: ',signif(time_taken,4),' seconds')
  
  if(return_all){
    output = list(image = image_out,
                  weight = weight_out,
                  inVar = inVar_out,
                  exp = exp_out,
                  cold = cold_out,
                  hot = hot_out,
                  clip = mask_clip,
                  image_pre_stack = pre_stack_image_list,
                  inVar_pre_stack = pre_stack_inVar_list,
                  exp_pre_stack = pre_stack_exp_list,
                  weight_pre_stack = pre_stack_weight_list,
                  which_overlap = which_overlap,
                  time = time_taken)
  }else{
    output = list(image = image_out,
                  weight = weight_out,
                  inVar = inVar_out,
                  exp = exp_out,
                  cold = cold_out,
                  hot = hot_out,
                  clip = mask_clip,
                  which_overlap = which_overlap,
                  time = time_taken)
  }
  class(output) = "ProMo"
  return(invisible(output))
}

.dilate_R = function(segim=NULL, size=3, shape='disc', expand='all', iters=1){
  kern = .makeBrush(size, shape=shape)
  if(expand[1] == 'all'){
    expand = 0
  }
  
  if(anyNA(segim)){
    NAmask = which(is.na(segim))
    segim[is.na(segim)] = 0L
  }else{
    NAmask = NULL
  }
  
  for(i in 1:iters){
    segim = .dilate_cpp(segim=segim, kern=kern, expand=expand)
  }
  
  if(!is.null(NAmask)){
    segim[NAmask] = NA
  }
  
  return(segim)
}

.makeBrush = function(size, shape=c('box', 'disc', 'diamond', 'Gaussian', 'line'), step=TRUE, sigma=0.3, angle=45) {
  #This is a direct port of EBImage::makeBrush. This reduces code dependencies, and EBImage does not appear to be well maintained.
  if(! (is.numeric(size) && (length(size)==1L) && (size>=1)) ) stop("'size' must be an odd integer.")
  shape = match.arg(arg = tolower(shape), choices = c('box', 'disc', 'diamond', 'gaussian', 'line'))
  
  if(size %% 2 == 0){
    size = size + 1
    warning(paste("'size' was rounded to the next odd number: ", size))
  }
  
  if (shape=='box') z = matrix(1L, size, size)
  else if (shape == 'line') {
    angle = angle %% 180
    angle.radians = angle * pi / 180;
    tg = tan(angle.radians)
    sizeh = (size-1)/2
    if ( angle < 45 || angle > 135) {
      z.x = sizeh
      z.y = round(sizeh*tg)
    }
    else {
      z.y = sizeh
      z.x = round(sizeh/tg)
    }
    z = array(0L, dim=2*c(z.x, z.y)+1);
    for (i in -sizeh:sizeh) {
      if ( angle < 45 || angle > 135) {
        ## scan horizontally
        i.x = i
        i.y = round(i*tg)
      }
      else {
        ## scan vertically
        i.y = i
        i.x = round(i/tg) 
      }
      z[i.x+z.x+1, i.y+z.y+1] = 1L
    }
  }
  else if (shape=='gaussian') {
    x = seq(-(size-1)/2, (size-1)/2, length=size)
    x = matrix(x, size, size)
    z = exp(- (x^2 + t(x)^2) / (2*sigma^2))
    z = z / sum(z)
  } else {
    ## pixel center coordinates
    x = 1:size -((size+1)/2)
    
    ## for each pixel, compute the distance from its center to the origin, using L1 norm ('diamond') or L2 norm ('disc')
    if (shape=='disc') {
      z = outer(x, x, FUN=function(X,Y) (X*X+Y*Y))
      mz = (size/2)^2
      z = (mz - z)/mz
      z = sqrt(ifelse(z>0, z, 0))
    } else {
      z = outer(x, x, FUN=function(X,Y) (abs(X)+abs(Y)))
      mz = (size/2)
      z = (mz - z)/mz
      z = ifelse(z>0, z, 0)
    }
    
    if (step) z = ifelse(z>0, 1L, 0L)
  }
  z
}
