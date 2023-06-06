.sph2car = function(long, lat, radius=1, deg=TRUE){
  if (is.matrix(long) || is.data.frame(long)) {
    if (ncol(long) == 1) {
      long = long[, 1]
    }
    else if (ncol(long) == 2) {
      lat = long[, 2]
      long = long[, 1]
    }
    else if (ncol(long) == 3) {
      radius = long[, 3]
      lat = long[, 2]
      long = long[, 1]
    }
  }
  if (missing(long) | missing(lat)) {
    stop("Missing full spherical 3D input data.")
  }
  if (deg) {
    long = long * pi/180
    lat = lat * pi/180
  }
  return = cbind(x = radius * cos(long) * cos(lat), y = radius * 
                   sin(long) * cos(lat), z = radius * sin(lat))
}

Rwcs_s2p = function(RA, Dec, keyvalues=NULL, pixcen='FITS', loc.diff=c(0,0), coord.type='deg',
                    sep=':', header=NULL, inherit=TRUE, WCSref=NULL, ctrl=2L, cores=1, ...){
  assertList(keyvalues, null.ok = TRUE)
  if(is.character(header) & is.null(keyvalues)){
    if(length(header) > 1){
      if(requireNamespace("Rfits", quietly = TRUE)){
        keyvalues = Rfits::Rfits_hdr_to_keyvalues(header)
        header = Rfits::Rfits_header_to_raw(Rfits::Rfits_keyvalues_to_header(keyvalues))
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
    Dec = RA[,2]
    RA = RA[,1]
  }
  
  if(coord.type=='sex'){
    RA = hms2deg(RA,sep=sep)
    Dec = dms2deg(Dec,sep=sep)
  }
  
  assertNumeric(RA)
  assertNumeric(Dec, len = length(RA))
  
  if(inherit){
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      header = options()$current_header
    }
    
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      keyvalues = options()$current_keyvalues
    }
    
    if(is.null(WCSref)){
      WCSref = options()$current_WCSref
    }
  }
  
  if(length(header)==1){
    
    nkey = nchar(header)/80
    
    if(!is.null(WCSref)){
      WCSref = tolower(WCSref)
      reflet = 1:26
      names(reflet) = letters
      if(! WCSref %in% letters){
        stop('WCS ref must be 0 (base WCS) or a letter [a-z]!')
      }
      WCSref = reflet[WCSref]
    }else{
      WCSref = 0
    }
    
    if(cores == 1L){
      output = Cwcs_head_s2p(
        RA = RA,
        Dec = Dec,
        header = header,
        nkey = nkey,
        WCSref = WCSref,
        ctrl=ctrl
      )
      
      if(is.null(dim(output))){
        good = which(output == 0)
        output = matrix(NA, length(RA), 2)
        if(length(good) > 0){
          output[good,] = Cwcs_head_s2p(
            RA = RA[good],
            Dec = Dec[good],
            header = header,
            nkey = nkey,
            WCSref = WCSref
          )
        }
      }
      
      if(anyInfinite(output)){ #catch for weird inversion problems
        bad = unique(which(is.infinite(output), arr.ind = TRUE)[,1])
        output[bad,] = Cwcs_head_s2p(
          RA = RA[bad] + 1e-12,
          Dec = Dec[bad] + 1e-12,
          header = header,
          nkey = nkey,
          WCSref = WCSref
        )
      }
    }else{
      registerDoParallel(cores=cores)
      
      maxlen = length(RA)
      chunk = ceiling(maxlen/cores)
      
      i = RAsub = Decsub = NULL
      
      RA = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        RA[lo:hi]
      }
      
      Dec = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        Dec[lo:hi]
      }
      
      output = foreach(RAsub = RA, Decsub = Dec, .combine='rbind')%dopar%{
        temp = Cwcs_head_s2p(
          RA = RAsub,
          Dec = Decsub,
          header = header,
          nkey = nkey,
          WCSref = WCSref,
          ctrl=ctrl
        )
        
        if(is.null(dim(temp))){
          good = which(temp == 0)
          temp = matrix(NA, length(RAsub), 2)
          if(length(good)>0){
            temp[good,] = Cwcs_head_s2p(
              RA = RAsub[good],
              Dec = Decsub[good],
              header = header,
              nkey = nkey,
              WCSref = WCSref
            )
          }
        }
        
        if(anyInfinite(temp)){ #catch for weird inversion problems
          bad = unique(which(is.infinite(temp), arr.ind = TRUE)[,1])
          temp[bad,] = Cwcs_head_s2p(
            RA = RAsub[bad] + 1e-12,
            Dec = Decsub[bad] + 1e-12,
            header = header,
            nkey = nkey,
            WCSref = WCSref
          )
        }
        
        return(temp)
      }
    }
  }else{
    
    keyvalues = Rwcs_keypass(keyvalues, ...)
    
    if(cores == 1L){
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
        PV1_0 = keyvalues$PV1_0,
        PV1_1 = keyvalues$PV1_1,
        PV1_2 = keyvalues$PV1_2,
        PV1_3 = keyvalues$PV1_3,
        PV1_4 = keyvalues$PV1_4,
        PV2_0 = keyvalues$PV2_0,
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
            PV1_0 = keyvalues$PV1_0,
            PV1_1 = keyvalues$PV1_1,
            PV1_2 = keyvalues$PV1_2,
            PV1_3 = keyvalues$PV1_3,
            PV1_4 = keyvalues$PV1_4,
            PV2_0 = keyvalues$PV2_0,
            PV2_1 = keyvalues$PV2_1,
            PV2_2 = keyvalues$PV2_2,
            PV2_3 = keyvalues$PV2_3,
            PV2_4 = keyvalues$PV2_4,
            PV2_5 = keyvalues$PV2_5
          )
        }
      }
    }else{
      registerDoParallel(cores=cores)
  
      maxlen = length(RA)
      chunk = ceiling(maxlen/cores)
      
      i = RAsub = Decsub = NULL
      
      RA = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        RA[lo:hi]
      }
      
      Dec = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        Dec[lo:hi]
      }
      
      output = foreach(RAsub = RA, Decsub = Dec, .combine='rbind')%dopar%{
        temp = Cwcs_s2p(
          RA = RAsub,
          Dec = Decsub,
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
          PV1_0 = keyvalues$PV1_0,
          PV1_1 = keyvalues$PV1_1,
          PV1_2 = keyvalues$PV1_2,
          PV1_3 = keyvalues$PV1_3,
          PV1_4 = keyvalues$PV1_4,
          PV2_0 = keyvalues$PV2_0,
          PV2_1 = keyvalues$PV2_1,
          PV2_2 = keyvalues$PV2_2,
          PV2_3 = keyvalues$PV2_3,
          PV2_4 = keyvalues$PV2_4,
          PV2_5 = keyvalues$PV2_5
        )
        
        if(is.null(dim(temp))){
          good = which(temp == 0)
          temp = matrix(NA, length(RAsub), 2)
          if(length(good)>0){
            temp[good,] = Cwcs_s2p(
              RA = RAsub[good],
              Dec = Decsub[good],
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
              PV1_0 = keyvalues$PV1_0,
              PV1_1 = keyvalues$PV1_1,
              PV1_2 = keyvalues$PV1_2,
              PV1_3 = keyvalues$PV1_3,
              PV1_4 = keyvalues$PV1_4,
              PV2_0 = keyvalues$PV2_0,
              PV2_1 = keyvalues$PV2_1,
              PV2_2 = keyvalues$PV2_2,
              PV2_3 = keyvalues$PV2_3,
              PV2_4 = keyvalues$PV2_4,
              PV2_5 = keyvalues$PV2_5
            )
          }
        }
        return(temp)
      }
    }
  }
  
  if(loc.diff[1] != 0){
    output[,1] = output[,1] - loc.diff[1]
  }
  if(loc.diff[2] != 0){
    output[,2] = output[,2] - loc.diff[2]
  }
  
  if(pixcen == 'R'){
    output[,1] = output[,1] - 0.5
    output[,2] = output[,2] - 0.5
  }
  
  colnames(output) = c('x','y')
  
  return(output)
}

Rwcs_p2s = function(x, y, keyvalues=NULL, pixcen='FITS', loc.diff=c(0,0), coord.type='deg',
                    sep=':', header=NULL, inherit=TRUE, WCSref=NULL, ctrl=2L, cores=1, ...){
  assertList(keyvalues, null.ok = TRUE)
  if(is.character(header) & is.null(keyvalues)){
    if(length(header) > 1){
      if(requireNamespace("Rfits", quietly = TRUE)){
        keyvalues = Rfits::Rfits_hdr_to_keyvalues(header)
        header = Rfits::Rfits_header_to_raw(Rfits::Rfits_keyvalues_to_header(keyvalues))
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
  
  if(loc.diff[1] != 0){
    x = x + loc.diff[1]
  }
  
  if(loc.diff[2] != 0){
    y = y + loc.diff[2]
  }
  
  assertNumeric(x)
  assertNumeric(y, len = length(x))
  
  if(inherit){
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      header = options()$current_header
    }
    
    if(is.null(keyvalues) & is.null(header) & length(list(...))==0){
      keyvalues = options()$current_keyvalues
    }
    
    if(is.null(WCSref)){
      WCSref = options()$current_WCSref
    }
  }
  
  if(length(header)==1){
    
    nkey = nchar(header)/80
    
    if(!is.null(WCSref)){
      WCSref = tolower(WCSref)
      reflet = 1:26
      names(reflet) = letters
      if(! WCSref %in% letters){
        stop('WCS ref must be 0 (base WCS) or a letter [a-z]!')
      }
      WCSref = reflet[WCSref]
    }else{
      WCSref = 0
    }
    
    if(cores == 1L){
      output = Cwcs_head_p2s(
        x = x,
        y = y,
        header = header,
        nkey = nkey,
        WCSref = WCSref,
        ctrl=ctrl
      )
      
      if(is.null(dim(output))){
        good = which(output == 0)
        output = matrix(NA, length(x), 2)
        if(length(good) > 0){
          output[good,] = Cwcs_head_p2s(
            x = x[good],
            y = y[good],
            header = header,
            nkey = nkey,
            WCSref = WCSref
          )
        }
      }
      
      if(anyInfinite(output)){ #catch for weird inversion problems
        bad = unique(which(is.infinite(output), arr.ind = TRUE)[,1])
        output[bad,] = Cwcs_head_p2s(
          x = x[bad] + 1e-6,
          y = y[bad] + 1e-6,
          header = header,
          nkey = nkey,
          WCSref = WCSref
        )
      }
    }else{
      registerDoParallel(cores=cores)
      
      maxlen = length(x)
      chunk = ceiling(maxlen/cores)
      
      i = xsub = ysub = NULL
      
      x = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        x[lo:hi]
      }
      
      y = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        y[lo:hi]
      }
      
      output = foreach(xsub = x, ysub = y, .combine='rbind')%dopar%{
        temp = Cwcs_head_p2s(
          x = xsub,
          y = ysub,
          header = header,
          nkey = nkey,
          WCSref = WCSref,
          ctrl=ctrl
        )
        
        if(is.null(dim(temp))){
          good = which(temp == 0)
          temp = matrix(NA, length(xsub), 2)
          if(length(good)>0){
            temp[good,] = Cwcs_head_p2s(
              x = xsub[good],
              y = ysub[good],
              header = header,
              nkey = nkey,
              WCSref = WCSref
            )
          }
        }
        
        if(anyInfinite(temp)){ #catch for weird inversion problems
          bad = unique(which(is.infinite(temp), arr.ind = TRUE)[,1])
          temp[bad,] = Cwcs_head_p2s(
            x = xsub[bad] + 1e-6,
            y = ysub[bad] + 1e-6,
            header = header,
            nkey = nkey,
            WCSref = WCSref
          )
        }
        
        return(temp)
      }
    }
  }else{
    
    keyvalues = Rwcs_keypass(keyvalues, ...)
    
    if(cores == 1L){
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
        PV1_0 = keyvalues$PV1_0,
        PV1_1 = keyvalues$PV1_1,
        PV1_2 = keyvalues$PV1_2,
        PV1_3 = keyvalues$PV1_3,
        PV1_4 = keyvalues$PV1_4,
        PV2_0 = keyvalues$PV2_0,
        PV2_1 = keyvalues$PV2_1,
        PV2_2 = keyvalues$PV2_2,
        PV2_3 = keyvalues$PV2_3,
        PV2_4 = keyvalues$PV2_4,
        PV2_5 = keyvalues$PV2_5
      )
      
      if(is.null(dim(output))){
        good = which(output == 0)
        output = matrix(NA, length(x), 2)
        if(length(good)>0){
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
            PV1_0 = keyvalues$PV1_0,
            PV1_1 = keyvalues$PV1_1,
            PV1_2 = keyvalues$PV1_2,
            PV1_3 = keyvalues$PV1_3,
            PV1_4 = keyvalues$PV1_4,
            PV2_0 = keyvalues$PV2_0,
            PV2_1 = keyvalues$PV2_1,
            PV2_2 = keyvalues$PV2_2,
            PV2_3 = keyvalues$PV2_3,
            PV2_4 = keyvalues$PV2_4,
            PV2_5 = keyvalues$PV2_5
          )
        }
      }
    }else{
      registerDoParallel(cores=cores)
      
      maxlen = length(x)
      chunk = ceiling(maxlen/cores)
      
      i = xsub = ysub = NULL
      
      x = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        x[lo:hi]
      }
      
      y = foreach(i = 1:cores)%do%{
        lo = (i - 1L)*chunk + 1L
        hi = min(lo + chunk - 1L, maxlen)
        y[lo:hi]
      }

      output = foreach(xsub = x, ysub = y, .combine='rbind')%dopar%{
        temp = Cwcs_p2s(
          x = xsub,
          y = ysub,
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
          PV1_0 = keyvalues$PV1_0,
          PV1_1 = keyvalues$PV1_1,
          PV1_2 = keyvalues$PV1_2,
          PV1_3 = keyvalues$PV1_3,
          PV1_4 = keyvalues$PV1_4,
          PV2_0 = keyvalues$PV2_0,
          PV2_1 = keyvalues$PV2_1,
          PV2_2 = keyvalues$PV2_2,
          PV2_3 = keyvalues$PV2_3,
          PV2_4 = keyvalues$PV2_4,
          PV2_5 = keyvalues$PV2_5
        )
        
        if(is.null(dim(temp))){
          good = which(temp == 0)
          temp = matrix(NA, length(xsub), 2)
          if(length(good)>0){
            temp[good,] = Cwcs_p2s(
              x = xsub[good],
              y = ysub[good],
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
              PV1_0 = keyvalues$PV1_0,
              PV1_1 = keyvalues$PV1_1,
              PV1_2 = keyvalues$PV1_2,
              PV1_3 = keyvalues$PV1_3,
              PV1_4 = keyvalues$PV1_4,
              PV2_0 = keyvalues$PV2_0,
              PV2_1 = keyvalues$PV2_1,
              PV2_2 = keyvalues$PV2_2,
              PV2_3 = keyvalues$PV2_3,
              PV2_4 = keyvalues$PV2_4,
              PV2_5 = keyvalues$PV2_5
            )
          }
        }
        return(temp)
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
                      CUNIT1='deg', CUNIT2='deg',
                      PV1_0=NULL, PV1_1=NULL, PV1_2=NULL, PV1_3=NULL, PV1_4=NULL,
                      #PV1_5=NULL, PV1_6=NULL, PV1_7=NULL, PV1_8=NULL, PV1_9=NULL, PV1_10=NULL,
                      PV2_0=NULL, PV2_1=NULL, PV2_2=NULL, PV2_3=NULL, PV2_4=NULL, PV2_5=NULL,
                      #PV2_6=NULL, PV2_7=NULL, PV2_8=NULL, PV2_9=NULL, PV2_10=NULL,
                      ...){
  if(!is.null(keyvalues)){
    if(missing(CTYPE1)){if(!is.null(keyvalues$CTYPE1)){CTYPE1 = keyvalues$CTYPE1}else{message('CTYPE1 is not defined!')}}
    if(missing(CTYPE2)){if(!is.null(keyvalues$CTYPE2)){CTYPE2 = keyvalues$CTYPE2}else{message('CTYPE2 is not defined!')}}
    if(missing(CRVAL1)){if(!is.null(keyvalues$CRVAL1)){CRVAL1 = keyvalues$CRVAL1}else{message('CRVAL1 is not defined!')}}
    if(missing(CRVAL2)){if(!is.null(keyvalues$CRVAL2)){CRVAL2 = keyvalues$CRVAL2}else{message('CRVAL2 is not defined!')}}
    if(missing(CRPIX1)){if(!is.null(keyvalues$CRPIX1)){CRPIX1 = keyvalues$CRPIX1}else{message('CRPIX1 is not defined!')}}
    if(missing(CRPIX2)){if(!is.null(keyvalues$CRPIX2)){CRPIX2 = keyvalues$CRPIX2}else{message('CRPIX2 is not defined!')}}
    if(missing(CUNIT1)){if(!is.null(keyvalues$CUNIT1)){CUNIT1 = keyvalues$CUNIT1}else{message('CUNIT1 is not defined!')}}
    if(missing(CUNIT2)){if(!is.null(keyvalues$CUNIT2)){CUNIT2 = keyvalues$CUNIT2}else{message('CUNIT2 is not defined!')}}
    if(missing(CD1_1)){
      if(!is.null(keyvalues$CD1_1)){
        CD1_1 = keyvalues$CD1_1
      }else{
        if((!is.null(keyvalues$CDELT1)) & (!is.null(keyvalues$PC1_1)) &  (!is.null(keyvalues$PC2_1))){
          CD1_1 = keyvalues$CDELT1 * keyvalues$PC1_1
          CD2_1 = keyvalues$CDELT1 * keyvalues$PC2_1
        }else if((!is.null(keyvalues$CDELT1)) & (!is.null(keyvalues$PC1_1)) &  (is.null(keyvalues$PC2_1))){
          CD1_1 = keyvalues$CDELT1 * keyvalues$PC1_1
          CD2_1 = 0
        }else if((!is.null(keyvalues$CDELT1)) & (!is.null(keyvalues$CROTA2))){
          CD1_1 = keyvalues$CDELT1 * cos(keyvalues$CROTA2*pi/180)
          CD2_1 = keyvalues$CDELT1 * sin(keyvalues$CROTA2*pi/180)
        }else if((!is.null(keyvalues$CDELT1)) & (is.null(keyvalues$CROTA2))){
          CD1_1 = keyvalues$CDELT1
          CD2_1 = 0 #for clarity
        }else{
          stop('CD1_1 and/or CD2_1 is not definable!')
        }
      }
    }
    if(missing(CD2_2)){
      if(!is.null(keyvalues$CD2_2)) {
        CD2_2 = keyvalues$CD2_2
      }else{
        if((!is.null(keyvalues$CDELT2)) & (!is.null(keyvalues$PC2_2)) & (!is.null(keyvalues$PC1_2))){
          CD2_2 = keyvalues$CDELT2 * keyvalues$PC2_2
          CD1_2 = keyvalues$CDELT2 * keyvalues$PC1_2
        }else if((!is.null(keyvalues$CDELT2)) & (!is.null(keyvalues$PC2_2)) & (is.null(keyvalues$PC1_2))){
          CD2_2 = keyvalues$CDELT2 * keyvalues$PC2_2
          CD1_2 = 0
        }else if((!is.null(keyvalues$CDELT2)) & (!is.null(keyvalues$CROTA2))){
          CD2_2 = keyvalues$CDELT2 * cos(keyvalues$CROTA2*pi/180)
          CD1_2 = -keyvalues$CDELT2 * sin(keyvalues$CROTA2*pi/180)
        }else  if((!is.null(keyvalues$CDELT2)) & (is.null(keyvalues$CROTA2))){
          CD2_2 = keyvalues$CDELT2
          CD1_2 = 0 #for clarity
        }else{
          stop('CD2_2 and/or CD1_2 is not definable!')
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
    if(CUNIT1=='        '){
      CUNIT1 = 'deg'
      message('CUNIT1 is blank, setting to \'deg\'!')
    }
    if(CUNIT2=='        '){
      CUNIT2 = 'deg'
      message('CUNIT1 is blank, setting to \'deg\'!')
    }
    
    if(is.null(PV1_0)){if(!is.null(keyvalues$PV1_0)){PV1_0 = keyvalues$PV1_0}else{PV1_0 = NA}}
    if(is.null(PV1_1)){if(!is.null(keyvalues$PV1_1)){PV1_1 = keyvalues$PV1_1}else{PV1_1 = NA}}
    if(is.null(PV1_2)){if(!is.null(keyvalues$PV1_2)){PV1_2 = keyvalues$PV1_2}else{PV1_2 = NA}}
    if(is.null(PV1_3)){if(!is.null(keyvalues$PV1_3)){PV1_3 = keyvalues$PV1_3}else{PV1_3 = NA}}
    if(is.null(PV1_4)){if(!is.null(keyvalues$PV1_4)){PV1_4 = keyvalues$PV1_4}else{PV1_4 = NA}}
    # Beyond this appears to be non-standard
    # if(is.null(PV1_5)){if(!is.null(keyvalues$PV1_5)){PV1_5 = keyvalues$PV1_5}else{PV1_5 = NA}}
    # if(is.null(PV1_6)){if(!is.null(keyvalues$PV1_6)){PV1_6 = keyvalues$PV1_6}else{PV1_6 = NA}}
    # if(is.null(PV1_7)){if(!is.null(keyvalues$PV1_7)){PV1_7 = keyvalues$PV1_7}else{PV1_7 = NA}}
    # if(is.null(PV1_8)){if(!is.null(keyvalues$PV1_8)){PV1_8 = keyvalues$PV1_8}else{PV1_8 = NA}}
    # if(is.null(PV1_9)){if(!is.null(keyvalues$PV1_9)){PV1_9 = keyvalues$PV1_9}else{PV1_9 = NA}}
    # if(is.null(PV1_10)){if(!is.null(keyvalues$PV1_10)){PV1_10 = keyvalues$PV1_10}else{PV1_10 = NA}}
    
    if(is.null(PV2_0)){if(!is.null(keyvalues$PV2_0)){PV2_0 = keyvalues$PV2_0}else{PV2_0 = NA}}
    if(is.null(PV2_1)){if(!is.null(keyvalues$PV2_1)){PV2_1 = keyvalues$PV2_1}else{PV2_1 = NA}}
    if(is.null(PV2_2)){if(!is.null(keyvalues$PV2_2)){PV2_2 = keyvalues$PV2_2}else{PV2_2 = NA}}
    if(is.null(PV2_3)){if(!is.null(keyvalues$PV2_3)){PV2_3 = keyvalues$PV2_3}else{PV2_3 = NA}}
    if(is.null(PV2_4)){if(!is.null(keyvalues$PV2_4)){PV2_4 = keyvalues$PV2_4}else{PV2_4 = NA}}
    if(is.null(PV2_5)){if(!is.null(keyvalues$PV2_5)){PV2_5 = keyvalues$PV2_5}else{PV2_5 = NA}}
    # Beyond this appears to be non-standard
    # if(is.null(PV2_6)){if(!is.null(keyvalues$PV2_6)){PV2_6 = keyvalues$PV2_6}else{PV2_6 = NA}}
    # if(is.null(PV2_7)){if(!is.null(keyvalues$PV2_7)){PV2_7 = keyvalues$PV2_7}else{PV2_7 = NA}}
    # if(is.null(PV2_8)){if(!is.null(keyvalues$PV2_8)){PV2_8 = keyvalues$PV2_8}else{PV2_8 = NA}}
    # if(is.null(PV2_9)){if(!is.null(keyvalues$PV2_9)){PV2_9 = keyvalues$PV2_9}else{PV2_9 = NA}}
    # if(is.null(PV2_10)){if(!is.null(keyvalues$PV2_10)){PV2_10 = keyvalues$PV2_10}else{PV2_10 = NA}}
  }else{
    keyvalues=list()
    
    if(is.null(PV1_0)){PV1_0 = NA}
    if(is.null(PV1_1)){PV1_1 = NA}
    if(is.null(PV1_2)){PV1_2 = NA}
    if(is.null(PV1_3)){PV1_3 = NA}
    if(is.null(PV1_4)){PV1_4 = NA}
    
    if(is.null(PV2_0)){PV2_0 = NA}
    if(is.null(PV2_1)){PV2_1 = NA}
    if(is.null(PV2_2)){PV2_2 = NA}
    if(is.null(PV2_3)){PV2_3 = NA}
    if(is.null(PV2_4)){PV2_4 = NA}
    if(is.null(PV2_5)){PV2_5 = NA}
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
  assertNumeric(PV1_0, len=1, null.ok = FALSE)
  assertNumeric(PV1_1, len=1, null.ok = FALSE)
  assertNumeric(PV1_2, len=1, null.ok = FALSE)
  assertNumeric(PV1_3, len=1, null.ok = FALSE)
  assertNumeric(PV1_4, len=1, null.ok = FALSE)
  # assertNumeric(PV1_5, len=1, null.ok = FALSE)
  # assertNumeric(PV1_6, len=1, null.ok = FALSE)
  # assertNumeric(PV1_7, len=1, null.ok = FALSE)
  # assertNumeric(PV1_8, len=1, null.ok = FALSE)
  # assertNumeric(PV1_9, len=1, null.ok = FALSE)
  # assertNumeric(PV1_10, len=1, null.ok = FALSE)
  assertNumeric(PV2_0, len=1, null.ok = FALSE)
  assertNumeric(PV2_1, len=1, null.ok = FALSE)
  assertNumeric(PV2_2, len=1, null.ok = FALSE)
  assertNumeric(PV2_3, len=1, null.ok = FALSE)
  assertNumeric(PV2_4, len=1, null.ok = FALSE)
  assertNumeric(PV2_5, len=1, null.ok = FALSE)
  # assertNumeric(PV2_6, len=1, null.ok = FALSE)
  # assertNumeric(PV2_7, len=1, null.ok = FALSE)
  # assertNumeric(PV2_8, len=1, null.ok = FALSE)
  # assertNumeric(PV2_9, len=1, null.ok = FALSE)
  # assertNumeric(PV2_10, len=1, null.ok = FALSE)
  
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
  keyvalues$CUNIT1 = CUNIT1
  keyvalues$CUNIT2 = CUNIT2
  keyvalues$PV1_0 = PV1_0
  keyvalues$PV1_1 = PV1_1
  keyvalues$PV1_2 = PV1_2
  keyvalues$PV1_3 = PV1_3
  keyvalues$PV1_4 = PV1_4
  # keyvalues$PV1_5 = PV1_5
  # keyvalues$PV1_6 = PV1_6
  # keyvalues$PV1_7 = PV1_7
  # keyvalues$PV1_8 = PV1_8
  # keyvalues$PV1_9 = PV1_9
  # keyvalues$PV1_10 = PV1_10
  keyvalues$PV2_0 = PV2_0
  keyvalues$PV2_1 = PV2_1
  keyvalues$PV2_2 = PV2_2
  keyvalues$PV2_3 = PV2_3
  keyvalues$PV2_4 = PV2_4
  keyvalues$PV2_5 = PV2_5
  # keyvalues$PV2_6 = PV2_6
  # keyvalues$PV2_7 = PV2_7
  # keyvalues$PV2_8 = PV2_8
  # keyvalues$PV2_9 = PV2_9
  # keyvalues$PV2_10 = PV2_10
  
  class(keyvalues) = 'Rfits_keylist'
  
  return(keyvalues)
}

Rwcs_pixscale = function(keyvalues=NULL, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1, type='old', dim=NULL){
  assertList(keyvalues, null.ok = TRUE)
  if(is.null(keyvalues)){
    keyvalues = options()$current_keyvalues
  }
  
  if(type=='old'){
    
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
  if(type=='new'){
    if(is.null(dim)){
      dim = c(keyvalues$NAXIS1, keyvalues$NAXIS2)
    }
    output = Rwcs_p2s(dim[1]/2 + c(-0.5,0.5), dim[2]/2 + c(-0.5,0.5), keyvalues = keyvalues)
    output[,1] = output[,1] * cos(mean(output[,2])*pi/180)
    return(2545.584412*sqrt(diff(output[,1])^2 + diff(output[,2])^2)) # 2545.584412 = 3600/sqrt(2)
  }
}

Rwcs_in_image = function(RA, Dec, xlim, ylim, buffer=0, plot=FALSE, style='points', pad=0, add=FALSE, ...){
  dots = list(...)
  
  if(missing(xlim) & !is.null(dots$keyvalues)){
    if(isTRUE(dots$keyvalues$ZIMAGE)){
      xlim = c(0, dots$keyvalues$ZNAXIS1)
    }else if(!is.null(dots$keyvalues$NAXIS1)){
      xlim = c(0, dots$keyvalues$NAXIS1)
    }else{
      stop('Missing NAXIS1 in keyvalues, please specify xlim manually!')
    }
  }
  
  if(missing(ylim) & !is.null(dots$keyvalues)){
    if(isTRUE(dots$keyvalues$ZIMAGE)){
      ylim = c(0, dots$keyvalues$ZNAXIS2)
    }else if(!is.null(dots$keyvalues$NAXIS2)){
      ylim = c(0, dots$keyvalues$NAXIS2)
    }else{
      stop('Missing NAXIS2 in keyvalues, please specify ylim manually!')
    }
  }
  
  suppressMessages({
    test_xy = Rwcs_s2p(RA=RA, Dec=Dec, pixcen='R', ...)
  })
  if(plot){
    if(add==FALSE){
      magplot(NA, NA, xlim=xlim + c(-pad,pad), ylim=ylim + c(-pad,pad), pch='.', asp=1, side=FALSE)
    }
    if(style=='points'){
      points(test_xy, col='red')
    }else if(style=='polygon'){
      polygon(test_xy, col=hsv(alpha=0.2), border='red')
    }else{
      stop('style must be points or polygon!')
    }
    if(add==FALSE){
      rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2])
      suppressMessages({
        Rwcs_grid(...)
        Rwcs_labels(...)
        Rwcs_compass(...)
      })
    }
    
  }
  return(as.logical(test_xy[,'x'] >= xlim[1] - buffer & test_xy[,'x'] <= xlim[2] + buffer & test_xy[,'y'] >= ylim[1] - buffer & test_xy[,'y'] <= ylim[2] + buffer))
}

Rwcs_overlap = function(keyvalues_test, keyvalues_ref=NULL, buffer=0, plot=FALSE, pad=0, add=FALSE){
  
  if(is.null(keyvalues_ref)){
    message('Using last provided keyvalues_ref!')
    keyvalues_ref = options()$current_keyvalues_ref
    if(is.null(keyvalues_ref)){
      stop('User must provide keyvalues_ref!')
    }
  }else{
    options(current_keyvalues_ref = keyvalues_ref)
  }
  
  if(isTRUE(keyvalues_test$ZIMAGE)){
    NAXIS1_test = keyvalues_test$ZNAXIS1
    NAXIS2_test = keyvalues_test$ZNAXIS2
  }else{
    NAXIS1_test = keyvalues_test$NAXIS1
    NAXIS2_test = keyvalues_test$NAXIS2
  }
  
  if(isTRUE(keyvalues_ref$ZIMAGE)){
    NAXIS1_ref = keyvalues_ref$ZNAXIS1
    NAXIS2_ref = keyvalues_ref$ZNAXIS2
  }else{
    NAXIS1_ref = keyvalues_ref$NAXIS1
    NAXIS2_ref = keyvalues_ref$NAXIS2
  }
  
  suppressMessages({
    pixscale_test = Rwcs_pixscale(keyvalues_test) # in asec
    pixscale_ref = Rwcs_pixscale(keyvalues_ref) # in asec
    
    centre_test = Rwcs_p2s(NAXIS1_test/2, NAXIS2_test/2, keyvalues=keyvalues_test, pixcen='R')
    centre_ref = Rwcs_p2s(NAXIS1_ref/2, NAXIS2_ref/2, keyvalues=keyvalues_ref, pixcen='R')
  })
  
  sph_test = .sph2car(centre_test)[1,]
  sph_ref = .sph2car(centre_ref)[1,]
  
  dot_prod = sph_test[1]*sph_ref[1] + sph_test[2]*sph_ref[2] + sph_test[3]*sph_ref[3]
  dot_prod[dot_prod < -1] = -1
  dot_prod[dot_prod > 1] = 1
  ang_sep = acos(dot_prod)/((pi/180)/3600) # in asec
  
  #using a 10% buffer to be safe
  max_sep = 1.1*(sqrt(NAXIS1_test^2 + NAXIS2_test^2)*pixscale_test + sqrt(NAXIS1_ref^2 + NAXIS2_ref^2)*pixscale_ref)/2
  
  if(ang_sep > max_sep){
    return(FALSE)
  }
  
  suppressMessages({
    left = Rwcs_p2s(rep(0,NAXIS2_test + 1L), 0:NAXIS2_test, keyvalues=keyvalues_test, pixcen='R')
    top = Rwcs_p2s(0:NAXIS1_test, rep(NAXIS2_test,NAXIS1_test + 1L), keyvalues=keyvalues_test, pixcen='R')
    right = Rwcs_p2s(rep(NAXIS1_test,NAXIS2_test + 1L), NAXIS2_test:0 , keyvalues=keyvalues_test, pixcen='R')
    bottom = Rwcs_p2s(NAXIS1_test:0, rep(0,NAXIS1_test + 1L), keyvalues=keyvalues_test, pixcen='R')
  })
  
  test_in_ref = any(Rwcs_in_image(RA=c(left[,'RA'], top[,'RA'], right[,'RA'], bottom[,'RA']), Dec=c(left[,'Dec'], top[,'Dec'], right[,'Dec'], bottom[,'Dec']), buffer=buffer, plot=plot, style='polygon', pad=pad, add=add, keyvalues=keyvalues_ref))
  
  if(test_in_ref){
    return(test_in_ref)
  }else{
    suppressMessages({
      left = Rwcs_p2s(rep(0,NAXIS2_ref + 1L), 0:NAXIS2_ref, keyvalues=keyvalues_ref, pixcen='R')
      top = Rwcs_p2s(0:NAXIS1_ref, rep(NAXIS2_ref,NAXIS1_ref + 1L), keyvalues=keyvalues_ref, pixcen='R')
      right = Rwcs_p2s(rep(NAXIS1_ref,NAXIS2_ref + 1L), NAXIS2_ref:0 , keyvalues=keyvalues_ref, pixcen='R')
      bottom = Rwcs_p2s(NAXIS1_ref:0, rep(0,NAXIS1_ref + 1L), keyvalues=keyvalues_ref, pixcen='R')
    })
    
    ref_in_test = any(Rwcs_in_image(RA=c(left[,'RA'], top[,'RA'], right[,'RA'], bottom[,'RA']), Dec=c(left[,'Dec'], top[,'Dec'], right[,'Dec'], bottom[,'Dec']), buffer=buffer, keyvalues=keyvalues_test))
    
    return(ref_in_test)
  }
}