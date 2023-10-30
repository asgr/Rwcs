Rwcs_imageRGB=function(R, G, B, keyvalues_out=NULL, Rkeyvalues=NULL, Gkeyvalues=NULL, Bkeyvalues=NULL,
                       dowarp='auto', direction = "auto", boundary = "dirichlet", interpolation = "cubic",
                       n, grid.col='grey', grid.lty=2, grid.lwd=0.5, lab.col='green', coord.type='sex',
                       margin=TRUE, loc.diff=c(0,0), xlab='Right Ascension', ylab='Declination',
                       mgp=c(2,0.5,0), mtline=2, position='topright', com.col="green", com.length=0.05,
                       coord.axis='auto', pretty='auto', decorate=TRUE, ...){
  
  if(missing(xlab)){
    if(coord.type=='sex'){
      xlab=paste(xlab,'/ H:M:S')
    }
    if(coord.type=='deg'){
      xlab=paste(xlab,'/ deg')
    }
  }
  if(missing(ylab)){
    if(coord.type=='sex'){
      ylab=paste(ylab,'/ D:M:S')
    }
    if(coord.type=='deg'){
      ylab=paste(ylab,'/ deg')
    }
  }
  
  if(!missing(R)){
    if(is.null(Rkeyvalues)){
      Rkeyvalues = R$keyvalues
    }
  }
  
  if(!missing(G)){
    if(is.null(Gkeyvalues)){
      Gkeyvalues = G$keyvalues
    }
  }
  
  if(!missing(B)){
    if(is.null(Bkeyvalues)){
      Bkeyvalues = B$keyvalues
    }
  }
  
  if(is.null(keyvalues_out)){
    if(exists('Rkeyvalues')){
      keyvalues_out = Rkeyvalues
      dim_out = dim(R)
    }else if(exists('Gkeyvalues')){
      keyvalues_out = Gkeyvalues
      dim_out = dim(G)
    }else if(exists('Bkeyvalues')){
      keyvalues_out = Bkeyvalues
      dim_out = dim(B)
    }
  }else{
    dim_out = NULL
  }
  
  if(dowarp=='auto'){
    dowarp=FALSE
    if(all(dim(R)==dim(G))==FALSE){dowarp=TRUE}
    if(all(dim(R)==dim(B))==FALSE){dowarp=TRUE}
    if(all(as.character(Rkeyvalues)==as.character(Gkeyvalues))==FALSE){dowarp=TRUE}
    if(all(as.character(Rkeyvalues)==as.character(Gkeyvalues))==FALSE){dowarp=TRUE}
  }
  
  if(dowarp){
    if(!requireNamespace("ProPane", quietly = TRUE)){
      stop("The ProPane package is needed for this function to work. Please install it from GitHub asgr/ProPane", call. = FALSE)
    }
    suppressMessages({
      Rim = ProPane::propaneWarp(image_in=R, keyvalues_out=keyvalues_out, keyvalues_in=Rkeyvalues, dim_out=dim_out, direction=direction, boundary=boundary, interpolation=interpolation)$imDat
      Gim = ProPane::propaneWarp(image_in=G, keyvalues_out=keyvalues_out, keyvalues_in=Gkeyvalues, dim_out=dim_out, direction=direction, boundary=boundary, interpolation=interpolation)$imDat
      Bim = ProPane::propaneWarp(image_in=B, keyvalues_out=keyvalues_out, keyvalues_in=Bkeyvalues, dim_out=dim_out, direction=direction, boundary=boundary, interpolation=interpolation)$imDat
    })
  }
  
  output = magimageRGB(R=Rim, G=Gim, B=Bim, axes=FALSE, ...)
  
  if(decorate){
    box()
    suppressMessages({
      Rwcs_grid(keyvalues=keyvalues_out, n=n, grid.col=grid.col, grid.lty=grid.lty, grid.lwd=grid.lwd,
                coord.type=coord.type, loc.diff=loc.diff, pretty=pretty)
      
      Rwcs_labels(keyvalues=keyvalues_out, n=n, lab.col=lab.col, coord.type=coord.type, margin=margin,
                  loc.diff=loc.diff, xlab=xlab, ylab=ylab, mgp=mgp, mtline=mtline,
                  pretty=pretty)
      
      Rwcs_compass(keyvalues=keyvalues_out, position=position, com.col=com.col, com.length=com.length,
                   loc.diff=loc.diff)
    })
  }
  
  return(invisible(output))
}
