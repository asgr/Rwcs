Rwcs_stack = function(image_list=NULL, inVar_list=NULL, exp_list=NULL, mask_list=NULL, magzero_in=0,
                      magzero_out=23.9, keyvalues_out=NULL, dim_out=NULL, cores=4, return_all=FALSE,
                      ...){
  
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for smoothing to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }
    
  registerDoParallel(cores=cores)
  
  Nim = length(image_list)
  if(return_all){
    Nbatch = Nim
  }else{
    Nbatch = cores
  }
  seq_process = seq(1,Nim,by=Nbatch)
  
  if(length(magzero_in) == 1){
    magzero_in = rep(magzero_in, Nim) 
  }
  
  zero_point_scale = rep(1, Nim)
  
  for(i in 1:length(magzero_in)){
    zero_point_scale[i] = 10^(-0.4*(magzero_in[i] - magzero_out))
  }
  
  dim_im = c(keyvalues_out$NAXIS1, keyvalues_out$NAXIS2)
  post_stack_image = matrix(0, dim_im[1], dim_im[2])
  post_stack_weight = matrix(0L, dim_im[1], dim_im[2])
  
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
  
  for(seq_start in seq_process){
    seq_end = (min(seq_start + Nbatch - 1L, Nim))
    message('Projecting Images ',seq_start,' to ',seq_end,' of ',Nim)
    
    pre_stack_image_list = NULL
    pre_stack_inVar_list = NULL
    
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
        return(Rwcs_warp(
          image_in = temp_inVar,
          keyvalues_out = keyvalues_out,
          dim_out = dim_out,
          ...
        )$imDat)
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
    
    if(is.null(pre_stack_inVar_list)){
      message('Stacking Images and InVar ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:length(pre_stack_image_list)){
        addID = which(!is.na(pre_stack_image_list[[i]]))
        post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]
        post_stack_weight[addID] = post_stack_weight[addID] + 1L
      }
    }else{
      for(i in 1:length(pre_stack_image_list)){
        addID = which(!is.na(pre_stack_image_list[[i]]) & is.finite(pre_stack_inVar_list[[i]]))
        post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]*pre_stack_inVar_list[[i]][addID]
        post_stack_weight[addID] = post_stack_weight[addID] + 1L
        post_stack_inVar[addID] = post_stack_inVar[addID] + pre_stack_inVar_list[[i]][addID]
      }
    }
    if(!is.null(pre_stack_exp_list)){
      message('Stacking Exposure Times ',seq_start,' to ',seq_end,' of ',Nim)
      for(i in 1:length(pre_stack_image_list)){
        addID = which(!is.na(pre_stack_exp_list[[i]]))
        post_stack_exp[addID] = post_stack_exp[addID] + pre_stack_exp_list[[i]][addID]
      }
    }
  }
  
  if(is.null(pre_stack_inVar_list)){
    post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_weight[post_stack_weight > 0]
  }else{
    post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_inVar[post_stack_weight > 0]
    post_stack_inVar[post_stack_weight == 0L] = NA
  }
  post_stack_image[post_stack_weight == 0L] = NA
  
  # post_stack = .internalStack(image_list=pre_stack_image,
  #                             skyRMS_list=pre_stack_skyRMS,
  #                             magzero_in=magzero_in,
  #                             magzero_out=magzero_out,
  #                             masking=masking
  #                             )
  
  keyvalues_out$EXTNAME = 'image'
  keyvalues_out$MAGZERO = magzero_out
  keyvalues_out$RWCS_VER = packageVersion('Rwcs')
  
  image_out = Rfits::Rfits_create_image(image=post_stack_image,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE,
                                 history='Stacked with Rwcs_stack')
  
  keyvalues_out$EXTNAME = 'weight'
  weight_out = Rfits::Rfits_create_image(image=post_stack_weight,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE)
  
  if(!is.null(pre_stack_inVar_list)){
    keyvalues_out$EXTNAME = 'inVar'
    inVar_out = Rfits::Rfits_create_image(image=post_stack_inVar,
                                    keyvalues=keyvalues_out,
                                    keypass=FALSE)
  }else{
    inVar_out = NULL
  }
  
  if(!is.null(pre_stack_exp_list)){
    keyvalues_out$EXTNAME = 'exp'
    exp_out = Rfits::Rfits_create_image(image=post_stack_exp,
                                          keyvalues=keyvalues_out,
                                          keypass=FALSE)
  }else{
    exp_out = NULL
  }
  
  if(return_all){
    output = list(image = image_out,
                  weight = weight_out,
                  inVar = inVar_out,
                  exp = exp_out,
                  image_pre_stack = pre_stack_image_list,
                  inVar_pre_stack = pre_stack_inVar_list,
                  exp_pre_stack = pre_stack_exp_list)
  }else{
    output = list(image = image_out,
                  weight = weight_out,
                  inVar = inVar_out,
                  exp = exp_out)
  }
  class(output) = "ProMo"
  return(invisible(output))
}
