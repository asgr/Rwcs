Rwcs_stack = function(image_list=NULL, inVar_list=NULL, exp_list=NULL, weight_list=NULL, mask_list=NULL, magzero_in=0,
                      magzero_out=23.9, keyvalues_out=NULL, dim_out=NULL, cores=4, Nbatch=cores,
                      keep_extreme_pix=FALSE, doclip=FALSE, clip_tol=100, clip_dilate=0, clip_sigma=5,
                      return_all=FALSE, ...){
  
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }
  
  timestart = proc.time()[3]
    
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
  
  seq_process = seq(1,Nim,by=Nbatch)
  
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
  
  dim_im = c(keyvalues_out$NAXIS1, keyvalues_out$NAXIS2)
  post_stack_image = matrix(0, dim_im[1], dim_im[2])
  mask_clip = NULL
  
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
    
    # if(Nim != Nbatch){
    #   stop('Nbatch must equal number of images')
    # }
    
    if(!is.null(inVar_list)){
      post_stack_cold_id = matrix(0L, dim_im[1], dim_im[2])
      post_stack_hot_id = matrix(0L, dim_im[1], dim_im[2])
    }else{
      stop('inVar_list required if doclip = TRUE!')
    }
  }
  
  for(seq_start in seq_process){
    seq_end = min(seq_start + Nbatch - 1L, Nim)
    message('Projecting Images ',seq_start,' to ',seq_end,' of ',Nim)
    
    Nbatch_sub = length(seq_start:seq_end)
    
    pre_stack_image_list = NULL
    pre_stack_inVar_list = NULL
    pre_stack_exp_list = NULL
    pre_stack_weight_list = NULL
    
    pre_stack_image_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar'))%dopar%{
      if(inherits(image_list[[i]], 'Rfits_pointer')){
        temp_image = image_list[[i]][,]
        temp_image$imDat = temp_image$imDat*zero_point_scale[i]
      }else{
        temp_image = image_list[[i]]
        temp_image$imDat = temp_image$imDat*zero_point_scale[i]
      }
      if(any(!is.finite(temp_image$imDat))){
        temp_image$imDat[!is.finite(temp_image$imDat)] = NA
      }
      if(!is.null(mask_list)){
        temp_image[mask_list[[i]] != 0] = NA
      }
      return(Rwcs_warp(
        image_in = temp_image,
        keyvalues_out = keyvalues_out,
        dim_out = dim_out,
        doscale = TRUE,
        ...
      )$imDat)
    }
    
    if(!is.null(inVar_list)){
      message('Projecting Inverse Variance ',seq_start,' to ',seq_end,' of ',Nim)
      
      pre_stack_inVar_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'pre_stack_image_list'))%dopar%{
        if(inherits(image_list[[i]], 'Rfits_pointer')){
          temp_inVar = image_list[[i]][,]
          temp_inVar$imDat[] = inVar_list[[i]]/(zero_point_scale[i]^2)
        }else{
          temp_inVar = image_list[[i]]
          temp_inVar$imDat[] = inVar_list[[i]]/(zero_point_scale[i]^2)
        }
        suppressMessages({
          return(Rwcs_warp(
            image_in = temp_inVar,
            keyvalues_out = keyvalues_out,
            dim_out = dim_out,
            doscale = FALSE,
            ...
            )$imDat*(Rwcs_pixscale(temp_inVar$keyvalues)^4 / Rwcs_pixscale(keyvalues_out)^4) #this is because RMS scales as linear pixelarea
          )
        })
      }
    }
    
    if(!is.null(exp_list)){
      message('Projecting Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
      if(length(exp_list) == 1){
        exp_list = rep(exp_list, Nim)
      }
      if(length(exp_list) != Nim){
        stop("Length of Exposure Times not equal to length of image_list!")
      }
      pre_stack_exp_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'pre_stack_image_list', 'pre_stack_inVar_list'))%dopar%{
        if(inherits(image_list[[i]], 'Rfits_pointer')){
          temp_exp = image_list[[i]][,]
          temp_exp$imDat[] = exp_list[[i]]
        }else{
          temp_exp = image_list[[i]]
          temp_exp$imDat[] = exp_list[[i]]
        }
        return(Rwcs_warp(
          image_in = temp_exp,
          keyvalues_out = keyvalues_out,
          dim_out = dim_out,
          doscale = FALSE,
          ...
        )$imDat)
      }
    }
    
    #new weight projections (if relevant)
    
    if(any(weight_image)){
      message('Projecting Weights ',seq_start,' to ',seq_end,' of ',Nim)
      if(length(weight_list) == 1){
        weight_list = rep(weight_list, Nim)
      }
      if(length(weight_list) != Nim){
        stop("Length of Weights not equal to length of image_list!")
      }
      pre_stack_weight_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'pre_stack_image_list', 'pre_stack_inVar_list', 'pre_stack_exp_list'))%dopar%{
        if(weight_image[i]){
          if(inherits(image_list[[i]], 'Rfits_pointer')){
            temp_weight = image_list[[i]][,]
            temp_weight$imDat[] = weight_list[[i]]
          }else{
            temp_weight = image_list[[i]]
            temp_weight$imDat[] = weight_list[[i]]
          }
          return(Rwcs_warp(
            image_in = temp_weight,
            keyvalues_out = keyvalues_out,
            dim_out = dim_out,
            doscale = FALSE,
            ...
          )$imDat)
        }else{
          return(weight_list[[i]])
        }
      }
    }else{
      pre_stack_weight_list = weight_list
    }
    
    if(is.null(pre_stack_inVar_list)){
      message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:Nbatch_sub){
        addID = which(!is.na(pre_stack_image_list[[i]]))
        post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]
        if(weight_image[i]){
          post_stack_weight[addID] = post_stack_weight[addID] + pre_stack_weight_list[[i]][addID]
        }else{
          post_stack_weight[addID] = post_stack_weight[addID] + weight_list[[i]]
        }
      }
    }else{
      for(i in 1:Nbatch_sub){
        addID = which(!is.na(pre_stack_image_list[[i]]) & is.finite(pre_stack_inVar_list[[i]]))
        post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]*pre_stack_inVar_list[[i]][addID]
        post_stack_inVar[addID] = post_stack_inVar[addID] + pre_stack_inVar_list[[i]][addID]
        if(weight_image[i]){
          post_stack_weight[addID] = post_stack_weight[addID] + pre_stack_weight_list[[i]][addID]
        }else{
          post_stack_weight[addID] = post_stack_weight[addID] + weight_list[[i]]
        }
      }
    }
    
    if(!is.null(pre_stack_exp_list)){
      message('Stacking Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:Nbatch_sub){
        addID = which(!is.na(pre_stack_exp_list[[i]]))
        post_stack_exp[addID] = post_stack_exp[addID] + pre_stack_exp_list[[i]][addID]
      }
    }
    
    if(keep_extreme_pix | doclip){
      for(i in 1:Nbatch_sub){
        new_cold = which(pre_stack_image_list[[i]] < post_stack_cold)
        new_hot = which(pre_stack_image_list[[i]] > post_stack_hot)
        
        post_stack_cold[new_cold] = pre_stack_image_list[[i]][new_cold]
        post_stack_hot[new_hot] = pre_stack_image_list[[i]][new_hot]
        
        if(doclip){
          post_stack_cold_id[new_cold] = seq_start + i - 1L
          post_stack_hot_id[new_hot] = seq_start + i - 1L
        }
      }
    }
  }
  
  if(return_all==FALSE){
    pre_stack_exp_list = NULL
  }
  
  if(is.null(pre_stack_inVar_list)){
    post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_weight[post_stack_weight > 0]
  }else{
    post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_inVar[post_stack_weight > 0]
    post_stack_inVar[post_stack_weight == 0L] = NA
  }
  
  if(keep_extreme_pix | doclip){
    post_stack_cold[post_stack_weight == 0L] = NA
    post_stack_hot[post_stack_weight == 0L] = NA
  }
  
  post_stack_image[post_stack_weight == 0L] = NA
  
  #Changed my mind on this- I think I should count exposure as a photon hitting a legal part of a sensor (bad pixel or not). I.e. fully masked regions with NA in the final image might still have a positive exposure time.
  # if(!is.null(post_stack_exp)){
  #   post_stack_exp[post_stack_weight == 0L] = NA
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
        message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
        for(i in 1:Nim){
          temp_mask_clip = (post_stack_cold_id == i) | (post_stack_hot_id == i)
          if(clip_dilate > 0){
            temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
          }
          addID = which(!is.na(pre_stack_image_list[[i]]) & temp_mask_clip==FALSE)
          post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]
          if(weight_image[i]){
            post_stack_weight[addID] = post_stack_weight[addID] + pre_stack_weight_list[[i]][addID]
          }else{
            post_stack_weight[addID] = post_stack_weight[addID] + weight_list[[i]]
          }
          
          if(keep_extreme_pix){
            mask_clip = mask_clip + temp_mask_clip
            new_cold = which(pre_stack_image_list[[i]] < post_stack_cold & temp_mask_clip==FALSE)
            new_hot = which(pre_stack_image_list[[i]] > post_stack_hot & temp_mask_clip==FALSE)
            post_stack_cold[new_cold] = pre_stack_image_list[[i]][new_cold]
            post_stack_hot[new_hot] = pre_stack_image_list[[i]][new_hot]
          }
        }
      }else{
        for(i in 1:Nim){
          temp_mask_clip = (post_stack_cold_id == i) | (post_stack_hot_id == i)
          if(clip_dilate > 0){
            temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
          }
          addID = which(!is.na(pre_stack_image_list[[i]]) & is.finite(pre_stack_inVar_list[[i]]) & temp_mask_clip==FALSE)
          post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]*pre_stack_inVar_list[[i]][addID]
          post_stack_inVar[addID] = post_stack_inVar[addID] + pre_stack_inVar_list[[i]][addID]
          if(weight_image[i]){
            post_stack_weight[addID] = post_stack_weight[addID] + pre_stack_weight_list[[i]][addID]
          }else{
            post_stack_weight[addID] = post_stack_weight[addID] + weight_list[[i]]
          }
          
          if(keep_extreme_pix){
            mask_clip = mask_clip + temp_mask_clip
            new_cold = which(pre_stack_image_list[[i]] < post_stack_cold & temp_mask_clip==FALSE)
            new_hot = which(pre_stack_image_list[[i]] > post_stack_hot & temp_mask_clip==FALSE)
            post_stack_cold[new_cold] = pre_stack_image_list[[i]][new_cold]
            post_stack_hot[new_hot] = pre_stack_image_list[[i]][new_hot]
          }
        }
      }
      
      # if(keep_extreme_pix){
      #   for(i in 1:Nim){
      #     temp_mask_clip = (post_stack_cold_id == i) | (post_stack_hot_id == i)
      #     if(clip_dilate > 0){
      #       temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
      #     }
      #     
      #     mask_clip = mask_clip + temp_mask_clip
      #     
      #     new_cold = which(pre_stack_image_list[[i]] < post_stack_cold & temp_mask_clip==FALSE)
      #     new_hot = which(pre_stack_image_list[[i]] > post_stack_hot & temp_mask_clip==FALSE)
      #     
      #     post_stack_cold[new_cold] = pre_stack_image_list[[i]][new_cold]
      #     post_stack_hot[new_hot] = pre_stack_image_list[[i]][new_hot]
      #   }
      # }
    }else{ # If we need to batch process the image_list then we need to re-project everything again
      #would maybe be neater to break this out a distinct function, but that is a job for another day...
      
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
      
      message('Reprojecting and restacking without clipped cold/hot pixels')
      for(seq_start in seq_process){
        seq_end = min(seq_start + Nbatch - 1L, Nim)
        message('Projecting Images ',seq_start,' to ',seq_end,' of ',Nim)
        
        Nbatch_sub = length(seq_start:seq_end)
        
        pre_stack_image_list = NULL
        pre_stack_inVar_list = NULL
        pre_stack_exp_list = NULL
        pre_stack_weight_list = NULL
        
        pre_stack_image_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar'))%dopar%{
          if(inherits(image_list[[i]], 'Rfits_pointer')){
            temp_image = image_list[[i]][,]
            temp_image$imDat = temp_image$imDat*zero_point_scale[i]
          }else{
            temp_image = image_list[[i]]
            temp_image$imDat = temp_image$imDat*zero_point_scale[i]
          }
          if(any(!is.finite(temp_image$imDat))){
            temp_image$imDat[!is.finite(temp_image$imDat)] = NA
          }
          if(!is.null(mask_list)){
            temp_image[mask_list[[i]] != 0] = NA
          }
          return(Rwcs_warp(
            image_in = temp_image,
            keyvalues_out = keyvalues_out,
            dim_out = dim_out,
            doscale = TRUE,
            ...
          )$imDat)
        }
        
        if(!is.null(inVar_list)){
          message('Projecting Inverse Variance ',seq_start,' to ',seq_end,' of ',Nim)
          
          pre_stack_inVar_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'pre_stack_image_list'))%dopar%{
            if(inherits(image_list[[i]], 'Rfits_pointer')){
              temp_inVar = image_list[[i]][,]
              temp_inVar$imDat[] = inVar_list[[i]]/(zero_point_scale[i]^2)
            }else{
              temp_inVar = image_list[[i]]
              temp_inVar$imDat[] = inVar_list[[i]]/(zero_point_scale[i]^2)
            }
            suppressMessages({
              return(Rwcs_warp(
                image_in = temp_inVar,
                keyvalues_out = keyvalues_out,
                dim_out = dim_out,
                doscale = FALSE,
                ...
                )$imDat*(Rwcs_pixscale(temp_inVar$keyvalues)^4 / Rwcs_pixscale(keyvalues_out)^4) #this is because RMS scales as linear pixelarea
              )
            })
          }
        }
        
        if(!is.null(exp_list)){
          message('Projecting Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
          if(length(exp_list) == 1){
            exp_list = rep(exp_list, Nim)
          }
          if(length(exp_list) != Nim){
            stop("Length of Exposure Times not equal to length of image_list!")
          }
          pre_stack_exp_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'pre_stack_image_list', 'pre_stack_inVar_list'))%dopar%{
            if(inherits(image_list[[i]], 'Rfits_pointer')){
              temp_exp = image_list[[i]][,]
              temp_exp$imDat[] = exp_list[[i]]
            }else{
              temp_exp = image_list[[i]]
              temp_exp$imDat[] = exp_list[[i]]
            }
            return(Rwcs_warp(
              image_in = temp_exp,
              keyvalues_out = keyvalues_out,
              dim_out = dim_out,
              doscale = FALSE,
              ...
            )$imDat)
          }
        }else{
          pre_stack_exp_list = NULL
        }
        
        if(any(weight_image)){
          message('Projecting Weights ',seq_start,' to ',seq_end,' of ',Nim)
          if(length(weight_list) == 1){
            weight_list = rep(weight_list, Nim)
          }
          if(length(weight_list) != Nim){
            stop("Length of Weights not equal to length of image_list!")
          }
          pre_stack_weight_list = foreach(i = seq_start:seq_end, .noexport=c('post_stack_image', 'post_stack_weight', 'post_stack_inVar', 'pre_stack_image_list', 'pre_stack_inVar_list', 'pre_stack_exp_list'))%dopar%{
            if(weight_image[i]){
              if(inherits(image_list[[i]], 'Rfits_pointer')){
                temp_weight = image_list[[i]][,]
                temp_weight$imDat[] = weight_list[[i]]
              }else{
                temp_weight = image_list[[i]]
                temp_weight$imDat[] = weight_list[[i]]
              }
              return(Rwcs_warp(
                image_in = temp_weight,
                keyvalues_out = keyvalues_out,
                dim_out = dim_out,
                doscale = FALSE,
                ...
              )$imDat)
            }else{
              return(weight_list[[i]])
            }
          }
        }else{
          pre_stack_weight_list = weight_list
        }
        
        if(is.null(pre_stack_inVar_list)){
          message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
          for(i in 1:Nbatch_sub){
            i_stack = seq_start + i - 1L
            
            temp_mask_clip = (post_stack_cold_id == i_stack) | (post_stack_hot_id == i_stack)
            if(clip_dilate > 0){
              temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
            }
            
            addID = which(!is.na(pre_stack_image_list[[i]]) & temp_mask_clip==FALSE)
            post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]
            if(weight_image[i]){
              post_stack_weight[addID] = post_stack_weight[addID] + pre_stack_weight_list[[i]][addID]
            }else{
              post_stack_weight[addID] = post_stack_weight[addID] + weight_list[[i]]
            }
          }
        }else{
          for(i in 1:Nbatch_sub){
            i_stack = seq_start + i - 1L
            
            temp_mask_clip = (post_stack_cold_id == i_stack) | (post_stack_hot_id == i_stack)
            if(clip_dilate > 0){
              temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
            }
            
            addID = which(!is.na(pre_stack_image_list[[i]]) & is.finite(pre_stack_inVar_list[[i]]) & temp_mask_clip==FALSE)
            post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]*pre_stack_inVar_list[[i]][addID]
            post_stack_inVar[addID] = post_stack_inVar[addID] + pre_stack_inVar_list[[i]][addID]
            if(weight_image[i]){
              post_stack_weight[addID] = post_stack_weight[addID] + pre_stack_weight_list[[i]][addID]
            }else{
              post_stack_weight[addID] = post_stack_weight[addID] + weight_list[[i]]
            }
          }
        }
        
        if(!is.null(pre_stack_exp_list)){
          message('Stacking Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
          for(i in 1:Nbatch_sub){
            addID = which(!is.na(pre_stack_exp_list[[i]]))
            post_stack_exp[addID] = post_stack_exp[addID] + pre_stack_exp_list[[i]][addID]
          }
        }
        
        if(keep_extreme_pix){
          for(i in 1:Nbatch_sub){
            i_stack = seq_start + i - 1L
            
            temp_mask_clip = (post_stack_cold_id == i_stack) | (post_stack_hot_id == i_stack)
            if(clip_dilate > 0){
              temp_mask_clip = .dilate_R(temp_mask_clip, size=clip_dilate)
            }
            
            mask_clip = mask_clip + temp_mask_clip
            
            new_cold = which(pre_stack_image_list[[i]] < post_stack_cold & temp_mask_clip==FALSE)
            new_hot = which(pre_stack_image_list[[i]] > post_stack_hot & temp_mask_clip==FALSE)
            
            post_stack_cold[new_cold] = pre_stack_image_list[[i]][new_cold]
            post_stack_hot[new_hot] = pre_stack_image_list[[i]][new_hot]
          }
        }
      }
    }
    
    if(return_all==FALSE){
      pre_stack_exp_list = NULL
    }
    
    if(is.null(pre_stack_inVar_list)){
      post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_weight[post_stack_weight > 0]
    }else{
      post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_inVar[post_stack_weight > 0]
      post_stack_inVar[post_stack_weight == 0L] = NA
    }
    
    if(keep_extreme_pix){
      post_stack_cold[post_stack_weight == 0L] = NA
      post_stack_hot[post_stack_weight == 0L] = NA
    }
    
    post_stack_image[post_stack_weight == 0L] = NA
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
  weight_out = Rfits::Rfits_create_image(image=post_stack_weight,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE)
  
  if(!is.null(post_stack_inVar)){
    keyvalues_out$EXTNAME = 'inVar'
    inVar_out = Rfits::Rfits_create_image(image=post_stack_inVar,
                                    keyvalues=keyvalues_out,
                                    keypass=FALSE)
  }else{
    inVar_out = NULL
  }
  
  if(!is.null(post_stack_exp)){
    keyvalues_out$EXTNAME = 'exp'
    exp_out = Rfits::Rfits_create_image(image=post_stack_exp,
                                          keyvalues=keyvalues_out,
                                          keypass=FALSE)
  }else{
    exp_out = NULL
  }
  
  if(keep_extreme_pix){
    post_stack_cold[weight_out$imDat == 0L] = NA
    
    keyvalues_out$EXTNAME = 'cold'
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
      mask_clip = Rfits::Rfits_create_image(image=mask_clip,
                                          keyvalues=keyvalues_out,
                                          keypass=FALSE)
    }
  }else{
    cold_out = NULL
    hot_out = NULL
    mask_clip = NULL
  }
  
  message('Time taken: ',signif(proc.time()[3] - timestart,4),' seconds')
  
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
                  weight_pre_stack = pre_stack_weight_list)
  }else{
    output = list(image = image_out,
                  weight = weight_out,
                  inVar = inVar_out,
                  exp = exp_out,
                  cold = cold_out,
                  hot = hot_out,
                  clip = mask_clip)
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
