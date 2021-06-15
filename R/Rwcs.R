Rwcs_s2p = function(RA, Dec, keyvalues=NULL, pixcen='FITS', loc.diff=c(0,0), coord.type='deg', sep=':', header=NULL, inherit=TRUE, ...){
  assertList(keyvalues, null.ok = TRUE)
  if(is.character(header) & is.null(keyvalues)){
    if(length(header) > 1){
      if(requireNamespace("Rfits", quietly = TRUE)){
        keyvalues = Rfits::Rfits_hdr_to_keyvalues(header)
      }else{
        stop("The Rfits package is need to process the header. Install from GitHub asgr/Rfits.")
      }
    }
  }
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
  
  if(inherit){
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      header = options()$current_header
    }
    
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      keyvalues = options()$current_keyvalues
    }
  }
  
  if(length(header)==1){
    
    output = Cwcs_head_s2p(
      RA = RA,
      Dec = Dec,
      header = header,
      nkeyrec = nchar(header)/80
    )
    
    if(is.null(dim(output))){
      good = which(output == 0)
      output = matrix(NA, length(RA), 2)
      if(length(good) > 0){
        output[good,] = Cwcs_head_s2p(
          RA = RA[good],
          Dec = Dec[good],
          header = header,
          nkeyrec = nchar(header)/80
        )
      }
    }
  }else{
    
    keyvalues = Rwcs_keypass(keyvalues, ...)
    
    output = Cwcs_s2p(
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
      RADESYS = keyvalues$RADESYS,
      EQUINOX = keyvalues$EQUINOX,
      PV1_1 = keyvalues$PV1_1,
      PV1_2 = keyvalues$PV1_2,
      PV1_3 = keyvalues$PV1_3,
      PV2_1 = keyvalues$PV2_1,
      PV2_2 = keyvalues$PV2_2,
      PV2_3 = keyvalues$PV2_3,
      PV2_4 = keyvalues$PV2_4,
      PV2_5 = keyvalues$PV2_5
    )
    
    if(is.null(dim(output))){
      good = which(output == 0)
      output = matrix(NA, length(RA), 2)
      if(length(good)>0){
        output[good,] = Cwcs_s2p(
          RA = RA[good],
          Dec = Dec[good],
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
          RADESYS = keyvalues$RADESYS,
          EQUINOX = keyvalues$EQUINOX,
          PV1_1 = keyvalues$PV1_1,
          PV1_2 = keyvalues$PV1_2,
          PV1_3 = keyvalues$PV1_3,
          PV2_1 = keyvalues$PV2_1,
          PV2_2 = keyvalues$PV2_2,
          PV2_3 = keyvalues$PV2_3,
          PV2_4 = keyvalues$PV2_4,
          PV2_5 = keyvalues$PV2_5
        )
      }
    }
  }
  
  output[,1]=output[,1]-loc.diff[1]
  output[,2]=output[,2]-loc.diff[2]
  if(pixcen == 'R'){
    output[,1]=output[,1]-0.5
    output[,2]=output[,2]-0.5
  }
  
  colnames(output)=c('x','y')
  
  return(output)
}

Rwcs_p2s = function(x, y, keyvalues=NULL, pixcen='FITS', loc.diff=c(0,0), coord.type='deg', sep=':', header=NULL, inherit=TRUE, ...){
  assertList(keyvalues, null.ok = TRUE)
  if(is.character(header) & is.null(keyvalues)){
    if(length(header) > 1){
      if(requireNamespace("Rfits", quietly = TRUE)){
        keyvalues = Rfits::Rfits_hdr_to_keyvalues(header)
      }else{
        stop("The Rfits package is need to process the header. Install from GitHub asgr/Rfits.")
      }
    }
  }
  assertChoice(pixcen, c('R','FITS'))
  assertNumeric(loc.diff, len=2)
  assertChoice(coord.type, c('deg','sex'))
  assertCharacter(sep, len=1)
  
  if(length(dim(x))==2){
    if(dim(x)[2]==2){
      y = x[,2]
      x = x[,1]
    }else{
      x = expand.grid(1:dim(x)[1], 1:dim(x)[2])
      y = x[,2] - 0.5
      x = x[,1] - 0.5
      pixcen = 'R'
    }
  }
  if(pixcen == 'R'){
    x = as.numeric(x) + 0.5
    y = as.numeric(y) + 0.5
  }
  x = x + loc.diff[1]
  y = y + loc.diff[2]
  
  assertNumeric(x)
  assertNumeric(y, len = length(x))
  
  if(inherit){
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      header = options()$current_header
    }
    
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      keyvalues = options()$current_keyvalues
    }
  }
  
  if(length(header)==1){
    
    output = Cwcs_head_p2s(
      x = x,
      y = y,
      header = header,
      nkeyrec = nchar(header)/80
    )
    
    if(is.null(dim(output))){
      good = which(output == 0)
      output = matrix(NA, length(x), 2)
      if(length(good) > 0){
        output[good,] = Cwcs_head_p2s(
          x = x[good],
          y = y[good],
          header = header,
          nkeyrec = nchar(header)/80
        )
      }
    }
  }else{
    
    keyvalues = Rwcs_keypass(keyvalues, ...)
    
    output = Cwcs_p2s(
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
      RADESYS = keyvalues$RADESYS,
      EQUINOX = keyvalues$EQUINOX,
      PV1_1 = keyvalues$PV1_1,
      PV1_2 = keyvalues$PV1_2,
      PV1_3 = keyvalues$PV1_3,
      PV2_1 = keyvalues$PV2_1,
      PV2_2 = keyvalues$PV2_2,
      PV2_3 = keyvalues$PV2_3,
      PV2_4 = keyvalues$PV2_4,
      PV2_5 = keyvalues$PV2_5
    )
    
    if(is.null(dim(output))){
      good = which(output == 0)
      output = matrix(NA, length(x), 2)
      if(length(good) > 0){
        output[good,] = Cwcs_p2s(
          x = x[good],
          y = y[good],
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
          RADESYS = keyvalues$RADESYS,
          EQUINOX = keyvalues$EQUINOX,
          PV1_1 = keyvalues$PV1_1,
          PV1_2 = keyvalues$PV1_2,
          PV1_3 = keyvalues$PV1_3,
          PV2_1 = keyvalues$PV2_1,
          PV2_2 = keyvalues$PV2_2,
          PV2_3 = keyvalues$PV2_3,
          PV2_4 = keyvalues$PV2_4,
          PV2_5 = keyvalues$PV2_5
        )
      }
    }
  }
  
  if(coord.type=='sex'){
    RAsex = deg2hms(output[,1], type='cat', sep=sep)
    Decsex = deg2dms(output[,2], type='cat', sep=sep)
    output = cbind(RAsex, Decsex)
  }
  
  colnames(output)=c('RA','Dec')
  
  return(output)
}

Rwcs_keypass=function(keyvalues=NULL,
                      CTYPE1='RA---TAN', CTYPE2='DEC--TAN',
                      CRVAL1=0, CRVAL2=0,
                      CRPIX1=0, CRPIX2=0,
                      CD1_1=1, CD1_2=0,
                      CD2_1=0, CD2_2=1,
                      RADESYS='ICRS',
                      EQUINOX='infer',
                      PV1_1=NA, PV1_2=NA, PV1_3=NA,
                      PV2_1=NA, PV2_2=NA, PV2_3=NA, PV2_4=NA, PV2_5=NA,
                      ...){
  if(!is.null(keyvalues)){
    if(missing(CTYPE1)){if(!is.null(keyvalues$CTYPE1)){CTYPE1 = keyvalues$CTYPE1}else{message('CTYPE1 is not defined!')}}
    if(missing(CTYPE2)){if(!is.null(keyvalues$CTYPE2)){CTYPE2 = keyvalues$CTYPE2}else{message('CTYPE2 is not defined!')}}
    if(missing(CRVAL1)){if(!is.null(keyvalues$CRVAL1)){CRVAL1 = keyvalues$CRVAL1}else{message('CRVAL1 is not defined!')}}
    if(missing(CRVAL2)){if(!is.null(keyvalues$CRVAL2)){CRVAL2 = keyvalues$CRVAL2}else{message('CRVAL2 is not defined!')}}
    if(missing(CRPIX1)){if(!is.null(keyvalues$CRPIX1)){CRPIX1 = keyvalues$CRPIX1}else{message('CRPIX1 is not defined!')}}
    if(missing(CRPIX2)){if(!is.null(keyvalues$CRPIX2)){CRPIX2 = keyvalues$CRPIX2}else{message('CRPIX2 is not defined!')}}
    if(missing(CD1_1)){
      if(!is.null(keyvalues$CD1_1)){
        CD1_1 = keyvalues$CD1_1
      }else{
        if((!is.null(keyvalues$CDELT1)) & (is.null(keyvalues$CROTA2))){
          CD1_1 = keyvalues$CDELT1
          CD2_1 = 0 #for clarity
        }else if((!is.null(keyvalues$CDELT1)) & (!is.null(keyvalues$CROTA2))){
          CD1_1 = keyvalues$CDELT1 * cos(keyvalues$CROTA2*pi/180)
          CD2_1 = keyvalues$CDELT1 * sin(keyvalues$CROTA2*pi/180)
        }else{
          stop('CD1_1 is not definable!')
        }
      }
    }
    if(missing(CD2_2)){
      if(!is.null(keyvalues$CD2_2)) {
        CD2_2 = keyvalues$CD2_2
      }else{
        if((!is.null(keyvalues$CDELT2)) & (is.null(keyvalues$CROTA2))){
          CD2_2 = keyvalues$CDELT2
          CD1_2 = 0 #for clarity
        }else if((!is.null(keyvalues$CDELT2)) & (!is.null(keyvalues$CROTA2))){
          CD2_2 = keyvalues$CDELT2 * cos(keyvalues$CROTA2*pi/180)
          CD1_2 = -keyvalues$CDELT2 * sin(keyvalues$CROTA2*pi/180)
        }else{
          stop('CD2_2 is not definable!')
        }
      }
    }
    if(missing(CD1_2)){
      if(!is.null(keyvalues$CD1_2)){
        CD1_2 = keyvalues$CD1_2
      }else{
        CD1_2 = 0
        message('CD1_2 is not definable, setting to 0!')
      }
    }
    if(missing(CD2_1)){
      if(!is.null(keyvalues$CD2_1)){
        CD2_1 = keyvalues$CD2_1
      }else{
        CD2_1 = 0
        message('CD2_1 is not definable, setting to 0!')
      }
    }
    if (missing(RADESYS)) {
      if (!is.null(keyvalues$RADESYS)) {
        RADESYS = keyvalues$RADESYS
      } else{
        if (!is.null(keyvalues$RADECSYS)) {
          RADESYS = keyvalues$RADECSYS
        }else{
          message('RADESYS is not defined (also no RADECSYS)!')
        }
      }
      if (!is.null(keyvalues$EQUINOX)) {
        EQUINOX = keyvalues$EQUINOX
      } else{
        if (!is.null(keyvalues$EPOCH)) {
          EQUINOX = keyvalues$EPOCH
        }else{
          message('EQUINOX is not defined (also no EPOCH)!')
        }
      }
    }
    if(missing(PV1_1)){if(!is.null(keyvalues$PV1_1)){PV1_1 = keyvalues$PV1_1}}
    if(missing(PV1_2)){if(!is.null(keyvalues$PV1_2)){PV1_2 = keyvalues$PV1_2}}
    if(missing(PV1_3)){if(!is.null(keyvalues$PV1_3)){PV1_3 = keyvalues$PV1_3}}
    if(missing(PV2_1)){if(!is.null(keyvalues$PV2_1)){PV2_1 = keyvalues$PV2_1}}
    if(missing(PV2_2)){if(!is.null(keyvalues$PV2_2)){PV2_2 = keyvalues$PV2_2}}
    if(missing(PV2_3)){if(!is.null(keyvalues$PV2_3)){PV2_3 = keyvalues$PV2_3}}
    if(missing(PV2_4)){if(!is.null(keyvalues$PV2_4)){PV2_4 = keyvalues$PV2_4}}
    if(missing(PV2_5)){if(!is.null(keyvalues$PV2_5)){PV2_5 = keyvalues$PV2_5}}
  }else{
    keyvalues=list()
  }
  
  if(EQUINOX == 'infer'){
    if(RADESYS %in% c('ICRS', 'FK5')){EQUINOX = 2000}else{EQUINOX = 1950}
  }
  
  allowed_proj = c(
    "AZP", #zenithal/azimuthal perspective
    "SZP", #slant zenithal perspective
    "TAN", #gnomonic
    "STG", #stereographic
    "SIN", #orthographic/synthesis
    "NCP", #unofficially supported SIN-like projection
    "ARC", #zenithal/azimuthal equidistant
    "ZPN", #zenithal/azimuthal polynomial
    "ZEA", #zenithal/azimuthal equal area
    "AIR", #Airy’s projection
    "CYP", #cylindrical perspective
    "CEA", #cylindrical equal area
    "CAR", #plate carrée
    "MER", #Mercator’s projection
    "COP", #conic perspective
    "COE", #conic equal area
    "COD", #conic equidistant
    "COO", #conic orthomorphic
    "SFL", #Sanson-Flamsteed (“global sinusoid”)
    "PAR", #parabolic
    "MOL", #Mollweide’s projection
    "AIT", #Hammer-Aitoff
    "BON", #Bonne’s projection
    "PCO", #polyconic
    "TSC", #tangential spherical cube
    "CSC", #COBE quadrilateralized spherical cube
    "QSC", #quadrilateralized spherical cube
    "HPX", #HEALPix
    "XPH"  #HEALPix polar, aka “butterfly”
  )
  
  allowed_axes = c(
    "RA", #right ascension
    "DEC", #declination
    "GLON", #galactic longitude
    "GLAT", #galactic latitude
    "ELON", #ecliptic longitude
    "ELAT", #ecliptic latitude
    "HLON", #helioecliptic longitude
    "HLAT", #helioecliptic latitude
    "SLON", #supergalactic longitude
    "SLAT" #supergalactic latitude
  )
  
  allowed_rade = c(
    "ICRS",
    "FK5",
    "FK4",
    "FK4-NO-E",
    "GAPPT"
  )
  
  assertCharacter(CTYPE1, len=1)
  assertCharacter(CTYPE2, len=1)
  if(grepl('-SIP', CTYPE1)){message('SIP not supported for CTYPE1 and ignored!'); CTYPE1=gsub('-SIP', '', CTYPE1)}
  if(grepl('-SIP', CTYPE2)){message('SIP not supported for CTYPE2 and ignored!'); CTYPE2=gsub('-SIP', '', CTYPE2)}
  if(grepl('-TPV', CTYPE1)){message('TPV not supported for CTYPE1 and ignored!'); CTYPE1=gsub('-TPV', '', CTYPE1)}
  if(grepl('-TPV', CTYPE2)){message('TPV not supported for CTYPE2 and ignored!'); CTYPE2=gsub('-TPV', '', CTYPE2)}
  if(grepl('-DSS', CTYPE1)){message('DSS not supported for CTYPE1 and ignored!'); CTYPE1=gsub('-DSS', '', CTYPE1)}
  if(grepl('-DSS', CTYPE2)){message('DSS not supported for CTYPE2 and ignored!'); CTYPE2=gsub('-DSS', '', CTYPE2)}
  if(grepl('-WAT', CTYPE1)){message('WAT not supported for CTYPE1 and ignored!'); CTYPE1=gsub('-WAT', '', CTYPE1)}
  if(grepl('-WAT', CTYPE2)){message('WAT not supported for CTYPE2 and ignored!'); CTYPE2=gsub('-WAT', '', CTYPE2)}
  if(grepl('-TPD', CTYPE1)){message('TPD not supported for CTYPE1 and ignored!'); CTYPE1=gsub('-TPD', '', CTYPE1)}
  if(grepl('-TPD', CTYPE2)){message('TPD not supported for CTYPE2 and ignored!'); CTYPE2=gsub('-TPD', '', CTYPE2)}
  if(nchar(CTYPE1) != 8){stop('CTYPE1 must be 8 characters!')}
  if(nchar(CTYPE2) != 8){stop('CTYPE2 must be 8 characters!')}
  split1=strsplit(CTYPE1, '-+')[[1]]
  split2=strsplit(CTYPE2, '-+')[[1]]
  assertCharacter(split1, len = 2)
  assertCharacter(split2, len = 2)
  assertChoice(split1[1], allowed_axes)
  assertChoice(split2[1], allowed_axes)
  assertChoice(split1[2], allowed_proj)
  assertChoice(split2[2], allowed_proj)
  assertNumeric(CRVAL1, len=1)
  assertNumeric(CRVAL2, len=1)
  assertNumeric(CRPIX1, len=1)
  assertNumeric(CRPIX2, len=1)
  assertNumeric(CD1_1, len=1)
  assertNumeric(CD1_2, len=1)
  assertNumeric(CD2_1, len=1)
  assertNumeric(CD2_2, len=1)
  assertCharacter(RADESYS, len=1)
  assertChoice(RADESYS, choices = allowed_rade)
  assertChoice(EQUINOX, choices = c(1950, 2000))
  assertNumeric(PV1_1, len=1)
  assertNumeric(PV1_2, len=1)
  assertNumeric(PV1_3, len=1)
  assertNumeric(PV2_2, len=1)
  assertNumeric(PV2_2, len=1)
  assertNumeric(PV2_3, len=1)
  assertNumeric(PV2_4, len=1)
  assertNumeric(PV2_5, len=1)
  
  keyvalues$CTYPE1 = CTYPE1
  keyvalues$CTYPE2 = CTYPE2
  keyvalues$CRVAL1 = CRVAL1
  keyvalues$CRVAL2 = CRVAL2
  keyvalues$CRPIX1 = CRPIX1
  keyvalues$CRPIX2 = CRPIX2
  keyvalues$CD1_1 = CD1_1
  keyvalues$CD1_2 = CD1_2
  keyvalues$CD2_1 = CD2_1
  keyvalues$CD2_2 = CD2_2
  keyvalues$RADESYS = RADESYS
  keyvalues$EQUINOX = EQUINOX
  keyvalues$PV1_1 = PV1_1
  keyvalues$PV1_2 = PV1_2
  keyvalues$PV1_3 = PV1_3
  keyvalues$PV2_1 = PV2_1
  keyvalues$PV2_2 = PV2_2
  keyvalues$PV2_3 = PV2_3
  keyvalues$PV2_4 = PV2_4
  keyvalues$PV2_5 = PV2_5
  return(keyvalues)
}

Rwcs_pixscale = function(keyvalues=NULL, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1){
  assertList(keyvalues, null.ok = TRUE)
  
  if(is.null(keyvalues)){
    keyvalues = options()$current_keyvalues
  }
  
  if(!is.null(keyvalues)){
    keyvalues = Rwcs_keypass(keyvalues)
    CD1_1 = keyvalues$CD1_1
    CD1_2 = keyvalues$CD1_2
    CD2_1 = keyvalues$CD2_1
    CD2_2 = keyvalues$CD2_2
  }
  assertNumeric(CD1_1, len=1)
  assertNumeric(CD1_2, len=1)
  assertNumeric(CD2_1, len=1)
  assertNumeric(CD2_2, len=1)
  return(3600*(sqrt(CD1_1^2+CD1_2^2)+sqrt(CD2_1^2+CD2_2^2))/2)
}