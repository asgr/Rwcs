Rwcs_stack = function(image_list=NULL, inVar_list=NULL, mask_list=NULL, magzero_in=0, magzero_out=0,
                      keyvalues_out=NULL, dim_out=NULL, cores=4, return_all=FALSE, ...){
  
  registerDoParallel(cores=cores)
  
  # if(!is.null(mask_list)){
  #   for(i in 1:length(image_list)){
  #     image_list[[i]]$imDat[mask_list[[i]] != 0] = NA
  #   }
  # }
  
  if(length(magzero_in) == 1){
    magzero_in = rep(magzero_in, length(image_list)) 
  }
  
  zero_point_scale = rep(1, length(image_list))
  
  for(i in 1:length(magzero_in)){
    zero_point_scale[i] = 10^(-0.4*(magzero_in[i] - magzero_out))
  }
  
  message('Projecting Images')
  pre_stack_image_list = foreach(i = 1:length(image_list))%dopar%{
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
    message('Projecting Inverse Variance')
    if(length(inVar_list) == 1){
      inVar_list = rep(inVar_list, length(image_list))
    }
    if(length(inVar_list) != length(image_list)){
      stop("Length of inVar_list not equal to length of image_list!")  
    }
    pre_stack_inVar_list = foreach(i = 1:length(image_list))%dopar%{
      if(inherits(image_list[[i]], 'Rfits_pointer')){
        temp_inVar = image_list[[i]][,]
        temp_inVar$imDat[] = inVar_list[i]/(zero_point_scale[i]^2)
      }else{
        temp_inVar = image_list[[i]]
        temp_inVar$imDat[] = inVar_list[i]/(zero_point_scale[i]^2)
      }
      return(Rwcs_warp(
        image_in = temp_inVar,
        keyvalues_out = keyvalues_out,
        dim_out = dim_out,
        ...
      )$imDat)
    }
  }else{
    pre_stack_inVar_list = NULL
  }
  
  message('Stacking Images')
  
  post_stack_image = matrix(0, dim(pre_stack_image_list[[1]])[1], dim(pre_stack_image_list[[1]])[2])
  post_stack_weight = matrix(0L, dim(pre_stack_image_list[[1]])[1], dim(pre_stack_image_list[[1]])[2])
  
  if(is.null(pre_stack_inVar_list)){
    for(i in 1:length(pre_stack_image_list)){
      addID = which(!is.na(pre_stack_image_list[[i]]))
      post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]
      post_stack_weight[addID] = post_stack_weight[addID] + 1L
    }
    post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_weight[post_stack_weight > 0]
    post_stack_inVar = NULL
  }else{
    post_stack_inVar = matrix(0, dim(pre_stack_image_list[[1]])[1], dim(pre_stack_image_list[[1]])[2])
    for(i in 1:length(pre_stack_image_list)){
      addID = which(!is.na(pre_stack_image_list[[i]]))
      post_stack_image[addID] = post_stack_image[addID] + pre_stack_image_list[[i]][addID]*pre_stack_inVar_list[[i]][addID]
      post_stack_weight[addID] = post_stack_weight[addID] + 1L
      post_stack_inVar[addID] = post_stack_inVar[addID] + pre_stack_inVar_list[[i]][addID]
    }
    post_stack_image[post_stack_weight > 0] = post_stack_image[post_stack_weight > 0]/post_stack_inVar[post_stack_weight > 0]
    post_stack_inVar[post_stack_weight==0L] = NA
  }
  
  post_stack_image[post_stack_weight==0L] = NA
  
  # post_stack = .internalStack(image_list=pre_stack_image,
  #                             skyRMS_list=pre_stack_skyRMS,
  #                             magzero_in=magzero_in,
  #                             magzero_out=magzero_out,
  #                             masking=masking
  #                             )
  
  keyvalues_out$EXTNAME = 'image'
  keyvalues_out$MAGZERO = magzero_out
  keyvalues_out$RWCS_VER = packageVersion('Rwcs')
  
  image_out = Rfits_create_image(image=post_stack_image,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE,
                                 history='Stacked with Rwcs_stack')
  
  keyvalues_out$EXTNAME = 'weight'
  weight_out = Rfits_create_image(image=post_stack_weight,
                                 keyvalues=keyvalues_out,
                                 keypass=FALSE)
  if(!is.null(pre_stack_inVar_list)){
    keyvalues_out$EXTNAME = 'inVar'
    inVar_out = Rfits_create_image(image=post_stack_inVar,
                                    keyvalues=keyvalues_out,
                                    keypass=FALSE)
  }else{
    inVar_out = NULL
  }
  if(return_all){
    return(invisible(list(image = image_out,
                          weight = weight_out,
                          inVar = inVar_out,
                          image_pre_stack = pre_stack_image_list,
                          inVar_pre_stack = pre_stack_inVar_list)))
  }else{
    return(invisible(list(image = image_out,
                          weight = weight_out,
                          inVar = inVar_out)))
  }
}
