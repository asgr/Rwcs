Rwcs_image=function(image, keyvalues=NULL, n, grid.col='grey', grid.lty=2, grid.lwd=0.5,
                    lab.col='green', coord.type='sex', margin=TRUE, loc.diff=c(0,0),
                    xlab='Right Ascension', ylab='Declination', mgp=c(2,0.5,0), mtline=2,
                    position='topright', com.col="green", com.length=0.05,
                    coord.axis='auto', pretty='auto', header=NULL, add=FALSE,  direction='auto', ...){
  
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
  
  if(!missing(image)){
    
    if(add){
      image = Rwcs_warp(image, keyvalues_in=keyvalues, header_in=header)
    }
    
    if(any(names(image)=='imDat') | any(names(image)=='image')){
      if(is.null(keyvalues)){
        keyvalues = image$keyvalues
      }
      
      if(is.null(header)){
        if(is.null(image$raw)){
          header = image$hdr
        }else{
          header = image$raw
        }
      }
      
      if(any(names(image)=='imDat')){
        image = image$imDat
      }
      
      if(any(names(image)=='image')){
        image = image$image
      }
    }
    
    # if(any(names(image)=='imDat') & is.null(keyvalues)){
    #   if(is.null(header)){
    #     if(is.null(image$raw)){
    #       header = image$hdr
    #     }else{
    #       header = image$raw
    #     }
    #   }
    #   keyvalues = image$keyvalues
    #   image = image$imDat
    # }else if(any(names(image)=='imDat') & !is.null(keyvalues)){
    #   image = image$imDat
    # }
    # if(any(names(image) == "image") & is.null(keyvalues)){
    #   if(is.null(header)){
    #     header = image$header
    #   }else{
    #     
    #   }
    #   keyvalues = image$keyvalues
    #   image = image$image
    # }else if(any(names(image) == "image") & !is.null(keyvalues)){
    #   image = image$image
    # }
    
    
    output = magimage(image, axes=FALSE, add=add, ...)
    if(add == FALSE){
      box()
    }
  }
  
  if(add == FALSE){
    suppressMessages({
      Rwcs_grid(keyvalues=keyvalues, n=n, grid.col=grid.col, grid.lty=grid.lty, 
                grid.lwd=grid.lwd, coord.type=coord.type, loc.diff=loc.diff, pretty=pretty, header=header)
        
      Rwcs_labels(keyvalues=keyvalues, n=n, lab.col=lab.col, coord.type=coord.type, 
                  margin=margin, loc.diff=loc.diff, xlab=xlab, ylab=ylab, mgp=mgp, 
                  mtline=mtline, pretty=pretty, header=header)
      
      Rwcs_compass(keyvalues=keyvalues, position=position, com.col=com.col,
                   com.length=com.length, loc.diff=loc.diff, header=header)
    })
    
    if(!is.null(header)){
      options(current_header = header)
    }else{
      options(current_header = NULL)
    }
      
    if(!is.null(keyvalues)){
      options(current_keyvalues = keyvalues)
    }else{
      options(current_keyvalues = NULL)
    }
  }
  return(invisible(output))
}

Rwcs_grid=function(keyvalues=NULL, n, grid.col='grey', grid.lty=2, grid.lwd=0.5, coord.type='sex',
                   loc.diff=c(0,0), pretty='auto', header=NULL, ...){
  
  xlo=min(par()$usr[1:2])
  xhi=max(par()$usr[1:2])
  ylo=min(par()$usr[3:4])
  yhi=max(par()$usr[3:4])
  
  coordlims=rbind(
    Rwcs_p2s(xlo, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xlo, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xhi, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xhi, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  )
  
  rarange=range(coordlims[,1])
  if(abs(diff(rarange))>180){
    rarange[rarange>180]=rarange[rarange>180]-360
    rarange=sort(rarange)
  }else{
    rarange=rarange %% 360
  }
  decrange=range(coordlims[,2])
  decrange=(decrange+90) %% 180 - 90
  
  if(pretty=='auto'){
    if(diff(rarange)>0.5){pretty=1}
    if(diff(rarange)<0.5 & diff(rarange)>0.5/6){pretty=5}
    if(diff(rarange)<0.5/6){pretty=300}
  }
  
  if(coord.type=='sex'){
    ragrid=maglab(rarange, n=n, prettybase = 5/pretty)
    decgrid=maglab(decrange, n=n, prettybase = 5/pretty)
  }
  if(coord.type=='deg'){
    ragrid=maglab(rarange, n=n)
    decgrid=maglab(decrange, n=n)
  }
  
  for(ra in ragrid$tickat){
    tempxy=Rwcs_s2p(cbind(ra, seq(min(decgrid$tickat), max(decgrid$tickat), len=100)),
                    keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
    tempxy[,1]=tempxy[,1]-loc.diff[1]
    tempxy[,2]=tempxy[,2]-loc.diff[2]
    lines(tempxy, col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
  for(dec in decgrid$tickat){
    tempxy=Rwcs_s2p(cbind(seq(min(ragrid$tickat), max(ragrid$tickat),len=100), dec),
                    keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
    tempxy[,1]=tempxy[,1]-loc.diff[1]
    tempxy[,2]=tempxy[,2]-loc.diff[2]
    lines(tempxy, col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
}

Rwcs_labels=function(keyvalues=NULL, n, lab.col='green', coord.type='sex', margin=TRUE,
                     loc.diff=c(0,0), xlab='Right Ascension', ylab='Declination',
                     mgp=c(2,0.5,0), mtline=2, coord.axis='auto', pretty='auto', header=NULL, ...){
  
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
  
  xlo=min(par()$usr[1:2])
  xhi=max(par()$usr[1:2])
  ylo=min(par()$usr[3:4])
  yhi=max(par()$usr[3:4])
  
  lenx = abs(xhi - xlo) + 1L
  leny = abs(yhi - ylo) + 1L
  
  coordlims=rbind(
    Rwcs_p2s(xlo, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xlo, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xhi, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xhi, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  )
  
  if(coord.axis[1]=='auto'){
    if(abs(diff(coordlims[1:2,1]))<abs(diff(coordlims[c(1,3),1]))){
      raaxis=1
      decaxis=2
    }else{
      raaxis=2
      decaxis=1
    }
  }else{
    raaxis=coord.axis[1]
    decaxis=coord.axis[2]
  }
  
  rarange=range(coordlims[,1])
  if(abs(diff(rarange))>180){
    rarange[rarange>180]=rarange[rarange>180]-360
    rarange=sort(rarange)
  }else{
    rarange=rarange %% 360
  }
  decrange=range(coordlims[,2])
  decrange=(decrange+90) %% 180 - 90
  
  if(pretty=='auto'){
    if(diff(rarange)>0.5){pretty=1}
    if(diff(rarange)<0.5 & diff(rarange)>0.5/6){pretty=5}
    if(diff(rarange)<0.5/6){pretty=300}
  }
  
  if(coord.type=='sex'){
    ragrid=maglab(rarange, n=n, prettybase = 5/pretty)
    decgrid=maglab(decrange, n=n, prettybase = 5/pretty)
  }
  if(coord.type=='deg'){
    ragrid=maglab(rarange, n=n)
    decgrid=maglab(decrange, n=n)
  }
  
  rapretty = ragrid$tickat
  rapretty = rapretty[rapretty>min(rarange) & rapretty<max(rarange)]
  rapretty = rapretty %% 360
  decpretty = decgrid$tickat
  decpretty = decpretty[decpretty>min(decrange) & decpretty<max(decrange)]
  
  dec_at_rapretty={}
  ra_at_decpretty={}
  
  coord_on_xaxis = Rwcs_p2s(xlo:xhi, rep(ylo,lenx), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  coord_on_yaxis = Rwcs_p2s(rep(xlo,leny), ylo:yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  
  if(raaxis==1){
    for(ra in rapretty){
      dec_add = coord_on_xaxis[which.min(abs(coord_on_xaxis[,1] - ra)),2]
      dec_at_rapretty = c(dec_at_rapretty, dec_add)
    }
    for(dec in decpretty){
      ra_add = coord_on_yaxis[which.min(abs(coord_on_yaxis[,2] - dec)),1]
      ra_at_decpretty = c(ra_at_decpretty, ra_add)
    }
  }else{
    for(ra in rapretty){
      dec_add = coord_on_yaxis[which.min(abs(coord_on_yaxis[,1] - ra)),2]
      dec_at_rapretty = c(dec_at_rapretty, dec_add)
    }
    for(dec in decpretty){
      ra_add = coord_on_xaxis[which.min(abs(coord_on_xaxis[,2] - dec)),1]
      ra_at_decpretty = c(ra_at_decpretty, ra_add)
    }
  }
  
  if(margin==FALSE){
    if(coord.type=='sex'){
      tempxy = Rwcs_s2p(cbind(rapretty, dec_at_rapretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(raaxis, at=tempxy[,raaxis], labels = deg2hms(rapretty, type='cat', digits=1), mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
      
      tempxy = Rwcs_s2p(cbind(ra_at_decpretty[-1], decpretty[-1]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(decaxis, at=tempxy[,decaxis], labels = deg2dms(decpretty[-1], type='cat', digits=0), mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
    }
    if(coord.type=='deg'){
      tempxy=Rwcs_s2p(cbind(rapretty, dec_at_rapretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(raaxis, at=tempxy[,raaxis], labels = rapretty, mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
      
      tempxy = Rwcs_s2p(cbind(ra_at_decpretty[-1], decpretty[-1]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(decaxis, at=tempxy[,decaxis], labels = decpretty[-1], mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
    }
    mtext(xlab, raaxis, line = -mtline, col=lab.col)
    mtext(ylab, decaxis, line = -mtline, col=lab.col)
  }else{
    if(coord.type=='sex'){
      tempxy = Rwcs_s2p(cbind(rapretty, dec_at_rapretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(raaxis, tempxy[,raaxis], labels=deg2hms(rapretty, type='cat', digits=1), mgp=mgp, tick=FALSE, ...)
      
      tempxy = Rwcs_s2p(cbind(ra_at_decpretty, decpretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(decaxis, tempxy[,decaxis], labels=deg2dms(decpretty, type='cat', digits=0), mgp=mgp, tick=FALSE, ...)
    }
    if(coord.type=='deg'){
      tempxy = Rwcs_s2p(cbind(rapretty, dec_at_rapretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(raaxis, tempxy[,raaxis], labels=rapretty, mgp=mgp, tick=FALSE, ...)
      
      tempxy = Rwcs_s2p(cbind(ra_at_decpretty, decpretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
      tempxy[,1] = tempxy[,1]-loc.diff[1]
      tempxy[,2] = tempxy[,2]-loc.diff[2]
      axis(decaxis, tempxy[,decaxis], labels=decpretty, mgp=mgp, tick=FALSE, ...)
    }
    mtext(xlab, raaxis, line = mtline)
    mtext(ylab, decaxis, line = mtline)
  }
}

Rwcs_compass=function(keyvalues = NULL, position='topright', com.col='green', com.length=0.05,
                      loc.diff=c(0,0), header=NULL, ...){
  xlo=min(par()$usr[1:2])
  xhi=max(par()$usr[1:2])
  ylo=min(par()$usr[3:4])
  yhi=max(par()$usr[3:4])
  
  xdiff=diff(c(xlo, xhi))
  ydiff=diff(c(ylo, yhi))
  
  if(position=='centre'){
    coord=Rwcs_p2s(xlo+xdiff*0.5, ylo+ydiff*0.5, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='bottom'){
    coord=Rwcs_p2s(xlo+xdiff*0.5, ylo+ydiff*0.15, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='bottomleft'){
    coord=Rwcs_p2s(xlo+xdiff*0.15, ylo+ydiff*0.15, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='left'){
    coord=Rwcs_p2s(xlo+xdiff*0.15, ylo+ydiff*0.5, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='topleft'){
    coord=Rwcs_p2s(xlo+xdiff*0.15, ylo+ydiff*0.85, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='top'){
    coord=Rwcs_p2s(xlo+xdiff*0.5, ylo+ydiff*0.85, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='topright'){
    coord=Rwcs_p2s(xlo+xdiff*0.85, ylo+ydiff*0.85, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='right'){
    coord=Rwcs_p2s(xlo+xdiff*0.85, ylo+ydiff*0.5, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  if(position=='bottomright'){
    coord=Rwcs_p2s(xlo+xdiff*0.85, ylo+ydiff*0.15, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)[1,]
  }
  
  startra=coord[1]
  startdec=coord[2]
  
  coordlims=rbind(
    Rwcs_p2s(xlo, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xlo, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xhi, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header),
    Rwcs_p2s(xhi, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  )
  
  rarange=range(coordlims[,1])
  if(abs(diff(rarange))>180){
    rarange[rarange>180]=rarange[rarange>180]-360
    rarange=sort(rarange)
  }else{
    rarange=rarange %% 360
  }
  decrange=range(coordlims[,2])
  decrange=(decrange+90) %% 180 - 90
  
  endra=startra+abs(rarange[2]-rarange[1])*0.05
  enddec=startdec+abs(decrange[2]-decrange[1])*0.05
  
  startxy=Rwcs_s2p(startra, startdec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  endxyN=Rwcs_s2p(startra, enddec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  endxyE=Rwcs_s2p(endra, startdec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  
  arrows(startxy[1,1], startxy[1,2], endxyN[1,1], endxyN[1,2], length=com.length, col=com.col, ...)
  arrows(startxy[1,1], startxy[1,2], endxyE[1,1], endxyE[1,2], length=com.length, col=com.col, ...)
  
  endra=startra+abs(rarange[2]-rarange[1])*0.065
  enddec=startdec+abs(decrange[2]-decrange[1])*0.065
  
  endxyN=Rwcs_s2p(startra, enddec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  endxyE=Rwcs_s2p(endra, startdec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff, header=header)
  
  text(endxyN[1,1], endxyN[1,2], labels='N', col=com.col, adj=c(0.5,0.5))
  text(endxyE[1,1], endxyE[1,2], labels='E', col=com.col, adj=c(0.5,0.5))
}
