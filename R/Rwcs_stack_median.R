Rwcs_stack_median = function(keyvalues_out = NULL,
                             dump_dir = NULL,
                             cores = 4,
                             chunk = 1e3,
                             doweight = TRUE){
  timestart = proc.time()[3]
  
  j = NULL
  
  assertList(keyvalues_out)
  assertCharacter(dump_dir, len=1)
  assertIntegerish(cores, len=1)
  assertIntegerish(chunk, len=1)
  
  if(!requireNamespace("Rfits", quietly = TRUE)){
    stop('The Rfits package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }
  
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for stacking to work. Please install from GitHub asgr/Rfits.', call. = FALSE)
  }
  
  registerDoParallel(cores=cores)
  
  image_files = list.files(dump_dir, full.names=TRUE, pattern = 'image_warp_')
  
  dump_list = foreach(i = 1:length(image_files))%dopar%{
    Rfits::Rfits_point(filename = image_files[i], ext=1, header=TRUE)
  }
  
  which_overlap = which(foreach(i = 1:length(dump_list), .combine='c')%dopar%{
    Rwcs_overlap(dump_list[[i]]$keyvalues, keyvalues_ref = keyvalues_out)
  })
  
  if(isTRUE(keyvalues_out$ZIMAGE)){
    NAXIS1 = keyvalues_out$ZNAXIS1
    NAXIS2 = keyvalues_out$ZNAXIS2
  }else{
    NAXIS1 = keyvalues_out$NAXIS1
    NAXIS2 = keyvalues_out$NAXIS2
  }
  
  stack_grid = expand.grid(seq(1,NAXIS1,by=chunk), seq(1,NAXIS2,by=chunk))
  stack_grid[,3] = stack_grid[,1] + chunk
  stack_grid[stack_grid[,3] > NAXIS1,3] = NAXIS1
  stack_grid[,4] = stack_grid[,2] + chunk
  stack_grid[stack_grid[,4] > NAXIS2,4] = NAXIS2
  
  stack_med = foreach(i = 1:dim(stack_grid)[1])%dopar%{
  #for(i in 1:dim(stack_grid)[1]){
    message('Stacking sub region ',i,' of ',dim(stack_grid)[1])
    
    xsub = as.integer(stack_grid[i, c(1,3)])
    ysub = as.integer(stack_grid[i, c(2,4)])
    
    keyvalues_sub = Rwcs_keyvalues_sub(keyvalues_out, xsub=xsub, ysub=ysub)
    
    temp_overlap = which(foreach(j = 1:length(dump_list), .combine = 'c')%dopar%{
      Rwcs_overlap(dump_list[[j]]$keyvalues, keyvalues_sub)
    })
    
    if(length(temp_overlap) == 0L){
      return(
        list(
          image = matrix(NA_real_, diff(range(xsub)) + 1L, diff(range(ysub)) + 1L),
          weight = matrix(0L, diff(range(xsub)) + 1L, diff(range(ysub)) + 1L)
        )
      )
    }
    
    image_list = foreach(j = temp_overlap)%do%{ #not sure why this won't work in dopar... Rfits race conditions?
      xrange = c(1,(diff(range(xsub)) + 1L)) + (xsub[1] - dump_list[[j]]$keyvalues$XCUTLO)
      yrange = c(1,(diff(range(ysub)) + 1L)) + (ysub[1] - dump_list[[j]]$keyvalues$YCUTLO)
      return(imager::as.cimg(dump_list[[j]][xrange, yrange]$imDat))
    }
    
    image = as.matrix(imager::parmed(image_list, na.rm=TRUE))
    
    if(doweight){
      weight_list = foreach(j = 1:length(image_list))%do%{
        return(imager::as.cimg(!is.na(image_list[[j]])))
      }
      
      weight = as.matrix(imager::add(weight_list))
    }else{
      weight = NULL
    }
    
    return(
      list(image = image, weight = weight)
    )
  }
  
  big_med = matrix(0, NAXIS1, NAXIS2)
  if(doweight){
    big_weight = matrix(0L, NAXIS1, NAXIS2)
  }else{
    big_weight = NULL
  }
  
  for(i in 1:dim(stack_grid)[1]){
    big_med[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_med[[i]]$image
    if(doweight){
      big_weight[stack_grid[i,1]:stack_grid[i,3], stack_grid[i,2]:stack_grid[i,4]] = stack_med[[i]]$weight
    }
  }
  
  keyvalues_out$EXTNAME = 'image'
  keyvalues_out$MAGZERO = dump_list[[1]]$keyvalues$MAGZERO
  keyvalues_out$R_VER = R.version$version.string
  keyvalues_out$RWCS_VER = as.character(packageVersion('Rwcs'))
  
  image_out = Rfits::Rfits_create_image(image=big_med,
                                        keyvalues=keyvalues_out,
                                        keypass=FALSE,
                                        history='Stacked with Rwcs_stack')
  rm(big_med)
  
  if(doweight){
    keyvalues_out$EXTNAME = 'weight'
    keyvalues_out$MAGZERO = NULL
    weight_out = Rfits::Rfits_create_image(image=big_weight,
                                           keyvalues=keyvalues_out,
                                           keypass=FALSE)
    
    rm(big_weight)
  }else{
    weight_out = NULL
  }
  
  time_taken = proc.time()[3] - timestart
  message('Time taken: ',signif(time_taken,4),' seconds')
  
  output = list(
    image = image_out,
    weight = weight_out,
    which_overlap = which_overlap,
    time = time_taken,
    Nim = length(which_overlap),
    dump_dir = dump_dir
  )
  
  class(output) = "ProMo"
  return(invisible(output))
}
