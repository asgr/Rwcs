Rwcs_image=function(image, keyvalues, n, grid.col='grey', grid.lty=2, grid.lwd=0.5,
                    lab.col='green', coord.type='sex', margin=TRUE, loc.diff=c(0,0),
                    xlab='Right Ascension', ylab='Declination', mgp=c(2,0.5,0), mtline=2,
                    position='topright', com.col="green", com.length=0.05,
                    coord.axis='auto', pretty='auto', ...){
  
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
    if(any(names(image)=='imDat') & missing(keyvalues)){
      keyvalues=image$keyvalues
      image=image$imDat
    }
    output=magimage(image, axes=FALSE, ...)
    box()
  }
  
  Rwcs_grid(keyvalues=keyvalues, n=n, grid.col=grid.col, grid.lty=grid.lty, 
            grid.lwd=grid.lwd, coord.type=coord.type, loc.diff=loc.diff, pretty=pretty)
  Rwcs_labels(keyvalues=keyvalues, n=n, lab.col=lab.col, coord.type=coord.type, 
              margin=margin, loc.diff=loc.diff, xlab=xlab, ylab=ylab, mgp=mgp, 
              mtline=mtline, pretty=pretty)
  Rwcs_compass(keyvalues=keyvalues, position=position, com.col=com.col,
               com.length=com.length, loc.diff=loc.diff)
  return=output
}

Rwcs_grid=function(keyvalues, n, grid.col='grey', grid.lty=2, grid.lwd=0.5, coord.type='sex',
                   loc.diff=c(0,0), pretty='auto', ...){
  
  xlo=min(par()$usr[1:2])
  xhi=max(par()$usr[1:2])
  ylo=min(par()$usr[3:4])
  yhi=max(par()$usr[3:4])
  
  coordlims=rbind(
    Rwcs_p2s(xlo, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xlo, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xhi, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xhi, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  )
  
  rarange=range(coordlims[,1])
  decrange=range(coordlims[,2])
  rarange=rarange %% 360
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
                    keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
    tempxy[,1]=tempxy[,1]-loc.diff[1]
    tempxy[,2]=tempxy[,2]-loc.diff[2]
    lines(tempxy, col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
  for(dec in decgrid$tickat){
    tempxy=Rwcs_s2p(cbind(seq(min(ragrid$tickat), max(ragrid$tickat),len=100), dec),
                    keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
    tempxy[,1]=tempxy[,1]-loc.diff[1]
    tempxy[,2]=tempxy[,2]-loc.diff[2]
    lines(tempxy, col=grid.col, lty=grid.lty, lwd=grid.lwd, ...)
  }
}

Rwcs_labels=function(keyvalues, n, lab.col='green', coord.type='sex', margin=TRUE,
                     loc.diff=c(0,0), xlab='Right Ascension', ylab='Declination',
                     mgp=c(2,0.5,0), mtline=2, coord.axis='auto', pretty='auto', ...){
  
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
  
  coordlims=rbind(
    Rwcs_p2s(xlo, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xlo, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xhi, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xhi, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
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
  decrange=range(coordlims[,2])
  rarange=rarange %% 360
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
  
  rapretty=ragrid$tickat
  rapretty=rapretty[rapretty>min(rarange) & rapretty<max(rarange)]
  decpretty=decgrid$tickat
  decpretty=decpretty[decpretty>min(decrange) & decpretty<max(decrange)]
  if(margin==FALSE){
    if(coord.type=='sex'){
      tempxy=Rwcs_s2p(cbind(rapretty, coordlims[1,2]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(raaxis, at=tempxy[,raaxis], labels = deg2hms(rapretty, type='cat', digits=1), mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
      
      tempxy=Rwcs_s2p(cbind(coordlims[1,1], decpretty[-1]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(decaxis, at=tempxy[,decaxis], labels = deg2dms(decpretty[-1], type='cat', digits=0), mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
    }
    if(coord.type=='deg'){
      tempxy=Rwcs_s2p(cbind(rapretty, coordlims[1,2]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(raaxis, at=tempxy[,raaxis], labels = rapretty, mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
      
      tempxy=Rwcs_s2p(cbind(coordlims[1,1], decpretty[-1]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(decaxis, at=tempxy[,decaxis], labels = decpretty[-1], mgp=-mgp-3, tick=FALSE, col.axis=lab.col, ...)
    }
    mtext(xlab, raaxis, line = -mtline, col=lab.col)
    mtext(ylab, decaxis, line = -mtline, col=lab.col)
  }else{
    if(coord.type=='sex'){
      tempxy=Rwcs_s2p(cbind(rapretty, coordlims[1,2]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(raaxis, tempxy[,raaxis], labels=deg2hms(rapretty, type='cat', digits=1), mgp=mgp, tick=FALSE, ...)
      
      tempxy=Rwcs_s2p(cbind(coordlims[1,1], decpretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(decaxis, tempxy[,decaxis], labels=deg2dms(decpretty, type='cat', digits=0), mgp=mgp, tick=FALSE, ...)
    }
    if(coord.type=='deg'){
      tempxy=Rwcs_s2p(cbind(rapretty, coordlims[1,2]), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(raaxis, tempxy[,raaxis], labels=rapretty, mgp=mgp, tick=FALSE, ...)
      
      tempxy=Rwcs_s2p(cbind(coordlims[1,1], decpretty), keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
      tempxy[,1]=tempxy[,1]-loc.diff[1]
      tempxy[,2]=tempxy[,2]-loc.diff[2]
      axis(decaxis, tempxy[,decaxis], labels=decpretty, mgp=mgp, tick=FALSE, ...)
    }
    mtext(xlab, raaxis, line = mtline)
    mtext(ylab, decaxis, line = mtline)
  }
}

Rwcs_compass=function(keyvalues, position='topright', com.col='green', com.length=0.05,
                      loc.diff=c(0,0), ...){
  xlo=min(par()$usr[1:2])
  xhi=max(par()$usr[1:2])
  ylo=min(par()$usr[3:4])
  yhi=max(par()$usr[3:4])
  
  xdiff=diff(c(xlo, xhi))
  ydiff=diff(c(ylo, yhi))
  
  if(position=='centre'){
    coord=Rwcs_p2s(xlo+xdiff*0.5, ylo+ydiff*0.5, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='bottom'){
    coord=Rwcs_p2s(xlo+xdiff*0.5, ylo+ydiff*0.15, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='bottomleft'){
    coord=Rwcs_p2s(xlo+xdiff*0.15, ylo+ydiff*0.15, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='left'){
    coord=Rwcs_p2s(xlo+xdiff*0.15, ylo+ydiff*0.5, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='topleft'){
    coord=Rwcs_p2s(xlo+xdiff*0.15, ylo+ydiff*0.85, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='top'){
    coord=Rwcs_p2s(xlo+xdiff*0.5, ylo+ydiff*0.85, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='topright'){
    coord=Rwcs_p2s(xlo+xdiff*0.85, ylo+ydiff*0.85, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='right'){
    coord=Rwcs_p2s(xlo+xdiff*0.85, ylo+ydiff*0.5, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  if(position=='bottomright'){
    coord=Rwcs_p2s(xlo+xdiff*0.85, ylo+ydiff*0.15, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)[1,]
  }
  
  startra=coord[1]
  startdec=coord[2]
  
  coordlims=rbind(
    Rwcs_p2s(xlo, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xlo, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xhi, ylo, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff),
    Rwcs_p2s(xhi, yhi, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  )
  rarange=range(coordlims[,1])
  decrange=range(coordlims[,2])
  
  endra=startra+abs(rarange[2]-rarange[1])*0.05
  enddec=startdec+abs(decrange[2]-decrange[1])*0.05
  
  startxy=Rwcs_s2p(startra, startdec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  endxyN=Rwcs_s2p(startra, enddec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  endxyE=Rwcs_s2p(endra, startdec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  
  arrows(startxy[1,1], startxy[1,2], endxyN[1,1], endxyN[1,2], length=com.length, col=com.col, ...)
  arrows(startxy[1,1], startxy[1,2], endxyE[1,1], endxyE[1,2], length=com.length, col=com.col, ...)
  
  endra=startra+abs(rarange[2]-rarange[1])*0.065
  enddec=startdec+abs(decrange[2]-decrange[1])*0.065
  
  endxyN=Rwcs_s2p(startra, enddec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  endxyE=Rwcs_s2p(endra, startdec, keyvalues=keyvalues, pixcen='R', loc.diff=loc.diff)
  
  text(endxyN[1,1], endxyN[1,2], labels='N', col=com.col, adj=c(0.5,0.5))
  text(endxyE[1,1], endxyE[1,2], labels='E', col=com.col, adj=c(0.5,0.5))
}