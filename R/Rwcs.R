Rwcs_s2p = function(RA, Dec, keyvalues=NULL, pixcen='FITS', loc.diff=c(0,0), coord.type='deg', sep=':', ...){
  
  assertList(keyvalues, null.ok = TRUE)
  assertChoice(pixcen, c('R','FITS'))
  assertNumeric(loc.diff, len=2)
  assertChoice(coord.type, c('deg','sex'))
  assertCharacter(sep, len=1)
  
  if(length(dim(RA))==2){
    Dec=RA[,2]
    RA=RA[,1]
  }
  if(coord.type=='sex'){RA=hms2deg(RA,sep=sep); Dec=dms2deg(Dec,sep=sep)}
  RA=as.numeric(RA)
  Dec=as.numeric(Dec)
  
  assertNumeric(RA)
  assertNumeric(Dec, len = length(RA))
  
  keyvalues = Rwcs_keypass(keyvalues, ...)
  
  output = t(Cwcs_s2p(
    RA = RA,
    Dec = Dec,
    CTYPE1 = keyvalues$CTYPE1,
    CTYPE2 = keyvalues$CTYPE2,
    CRVAL1 = keyvalues$CRVAL1,
    CRVAL2 = keyvalues$CRVAL2,
    CRPIX1 = keyvalues$CRPIX1,
    CRPIX2 = keyvalues$CRPIX2,
    CD1_1 = keyvalues$CD1_1,
    CD1_2 = keyvalues$CD1_2,
    CD2_1 = keyvalues$CD2_1,
    CD2_2 = keyvalues$CD2_2,
    PV1 = keyvalues$PV1,
    PV2 = keyvalues$PV2
  ))

  output[,1]=output[,1]-loc.diff[1]
  output[,2]=output[,2]-loc.diff[2]
  if(pixcen == 'R'){
    output[,1]=output[,1]-0.5
    output[,2]=output[,2]-0.5
  }
  
  colnames(output)=c('x','y')
  
  return(output)
}

Rwcs_p2s = function(x, y, keyvalues=NULL, pixcen='FITS', loc.diff=c(0,0), coord.type='deg', sep=':', ...){
  
  assertList(keyvalues, null.ok = TRUE)
  assertChoice(pixcen, c('R','FITS'))
  assertNumeric(loc.diff, len=2)
  assertChoice(coord.type, c('deg','sex'))
  assertCharacter(sep, len=1)
  
  if(length(dim(x))==2){
    y = x[,2]
    x = x[,1]
  }
  if(pixcen == 'R'){
    x = as.numeric(x) + 0.5
    y = as.numeric(y) + 0.5
  }
  x = x + loc.diff[1]
  y = y + loc.diff[2]
  
  assertNumeric(x)
  assertNumeric(y, len = length(x))
  
  keyvalues = Rwcs_keypass(keyvalues, ...)
  
  output = t(Cwcs_p2s(
    x = x,
    y = y,
    CTYPE1 = keyvalues$CTYPE1,
    CTYPE2 = keyvalues$CTYPE2,
    CRVAL1 = keyvalues$CRVAL1,
    CRVAL2 = keyvalues$CRVAL2,
    CRPIX1 = keyvalues$CRPIX1,
    CRPIX2 = keyvalues$CRPIX2,
    CD1_1 = keyvalues$CD1_1,
    CD1_2 = keyvalues$CD1_2,
    CD2_1 = keyvalues$CD2_1,
    CD2_2 = keyvalues$CD2_2,
    PV1 = keyvalues$PV1,
    PV2 = keyvalues$PV2
  ))
  
  if(coord.type=='sex'){
    RAsex = deg2hms(output[,1], sep=sep)
    Decsex = deg2dms(output[,2], sep=sep)
    output = cbind(RAsex, Decsex)
  }
  
  colnames(output)=c('RA','Dec')
  
  return(output)
}

Rwcs_keypass=function(keyvalues=NULL, CRVAL1=0, CRVAL2=0, CRPIX1=0, CRPIX2=0, CD1_1=1,
                      CD1_2=0, CD2_1=0, CD2_2=1, PV1=0, PV2=0, CTYPE1='RA---TAN',
                      CTYPE2='DEC--TAN'){
  if(!is.null(keyvalues)){
    if(missing(CRVAL1)){if(!is.null(keyvalues$CRVAL1)){CRVAL1 = keyvalues$CRVAL1}else{message('CRVAL1 is not defined!')}}
    if(missing(CRVAL2)){if(!is.null(keyvalues$CRVAL2)){CRVAL2 = keyvalues$CRVAL2}else{message('CRVAL2 is not defined!')}}
    if(missing(CRPIX1)){if(!is.null(keyvalues$CRPIX1)){CRPIX1 = keyvalues$CRPIX1}else{message('CRPIX1 is not defined!')}}
    if(missing(CRPIX2)){if(!is.null(keyvalues$CRPIX2)){CRPIX2 = keyvalues$CRPIX2}else{message('CRPIX2 is not defined!')}}
    if(missing(CD1_1)){if(!is.null(keyvalues$CD1_1)){CD1_1 = keyvalues$CD1_1}else{if(!is.null(keyvalues$CDELT1)){CD1_1 = keyvalues$CDELT1}else{message('CD1_1 is not defined!')}}}
    if(missing(CD1_2)){if(!is.null(keyvalues$CD1_2)){CD1_2 = keyvalues$CD1_2}else{message('CD1_2 is not defined!')}}
    if(missing(CD2_1)){if(!is.null(keyvalues$CD2_1)){CD2_1 = keyvalues$CD2_1}else{message('CD2_1 is not defined!')}}
    if(missing(CD2_2)){if(!is.null(keyvalues$CD2_2)){CD2_2 = keyvalues$CD2_2}else{if(!is.null(keyvalues$CDELT2)){CD2_2 = keyvalues$CDELT2}else{message('CD2_2 is not defined!')}}}
    if(missing(PV1)){if(!is.null(keyvalues$PV1)){PV1 = keyvalues$PV1}}
    if(missing(PV2)){if(!is.null(keyvalues$PV2)){PV2 = keyvalues$PV2}}
    if(missing(CTYPE1)){if(!is.null(keyvalues$CTYPE1)){CTYPE1 = keyvalues$CTYPE1}else{message('CTYPE1 is not defined!')}}
    if(missing(CTYPE2)){if(!is.null(keyvalues$CTYPE2)){CTYPE2 = keyvalues$CTYPE2}else{message('CTYPE2 is not defined!')}}
  }else{
    keyvalues=list()
  }
  
  assertNumeric(CRVAL1, len=1)
  assertNumeric(CRVAL2, len=1)
  assertNumeric(CRPIX1, len=1)
  assertNumeric(CRPIX2, len=1)
  assertNumeric(CD1_1, len=1)
  assertNumeric(CD1_2, len=1)
  assertNumeric(CD2_1, len=1)
  assertNumeric(CD2_2, len=1)
  assertNumeric(PV1, len=1)
  assertNumeric(PV2, len=1)
  assertCharacter(CTYPE1, len=1)
  assertCharacter(CTYPE2, len=1)
  if(nchar(CTYPE1) != 8){stop('CTYPE1 must be 8 characters!')}
  if(nchar(CTYPE2) != 8){stop('CTYPE2 must be 8 characters!')}
  
  keyvalues$CRVAL1 = CRVAL1
  keyvalues$CRVAL2 = CRVAL2
  keyvalues$CRPIX1 = CRPIX1
  keyvalues$CRPIX2 = CRPIX2
  keyvalues$CD1_1 = CD1_1
  keyvalues$CD1_2 = CD1_2
  keyvalues$CD2_1 = CD2_1
  keyvalues$CD2_2 = CD2_2
  keyvalues$PV1 = PV1
  keyvalues$PV2 = PV2
  keyvalues$CTYPE1 = CTYPE1
  keyvalues$CTYPE2 = CTYPE2
  return(keyvalues)
}

