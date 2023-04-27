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
    if(any(names(R)=='imDat') | any(names(R)=='image')){
      if(is.null(Rkeyvalues)){
        Rkeyvalues = R$keyvalues
      }
      
      if(any(names(R)=='imDat')){
        Rim = R$imDat
      }
      
      if(any(names(image)=='image')){
        Rim = R$image
      }
    }
  }
  
  if(!missing(G)){
    if(any(names(G)=='imDat') | any(names(G)=='image')){
      if(is.null(Gkeyvalues)){
        Gkeyvalues = G$keyvalues
      }
      
      if(any(names(G)=='imDat')){
        Gim = G$imDat
      }
      
      if(any(names(image)=='image')){
        Gim = G$image
      }
    }
  }
  
  if(!missing(B)){
    if(any(names(B)=='imDat') | any(names(B)=='image')){
      if(is.null(Bkeyvalues)){
        Bkeyvalues = B$keyvalues
      }
      
      if(any(names(B)=='imDat')){
        Bim = B$imDat
      }
      
      if(any(names(image)=='image')){
        Bim = B$image
      }
    }
  }
  
  if(is.null(keyvalues_out)){
    if(exists('Rkeyvalues')){
      keyvalues_out = Rkeyvalues
      dim_out = dim(Rim)
    }else if(exists('Gkeyvalues')){
      keyvalues_out = Gkeyvalues
      dim_out = dim(Gim)
    }else if(exists('Bkeyvalues')){
      keyvalues_out = Bkeyvalues
      dim_out = dim(Bim)
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
    suppressMessages({
      Rim = Rwcs_warp(image_in=Rim, keyvalues_out=keyvalues_out, keyvalues_in=Rkeyvalues, dim_out=dim_out, direction=direction, boundary=boundary, interpolation=interpolation)$imDat
      Gim = Rwcs_warp(image_in=Gim, keyvalues_out=keyvalues_out, keyvalues_in=Gkeyvalues, dim_out=dim_out, direction=direction, boundary=boundary, interpolation=interpolation)$imDat
      Bim = Rwcs_warp(image_in=Bim, keyvalues_out=keyvalues_out, keyvalues_in=Bkeyvalues, dim_out=dim_out, direction=direction, boundary=boundary, interpolation=interpolation)$imDat
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
