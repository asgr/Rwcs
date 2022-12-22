Rwcs_setscreen = function(CRVAL1=0, CRVAL2=0,
                          pixscale=1,
                          NAXIS1=1000, NAXIS2=1000,
                          CRPIX1 = NAXIS1/2 + 0.5, CRPIX2 = NAXIS2/2 + 0.5, # +0.5 for FITS standard
                          CTYPE1='RA---TAN', CTYPE2='DEC--TAN',
                          CD1_1=-pixscale/3600, CD1_2=0, CD2_1=0, CD2_2=pixscale/3600,
                          CUNIT1 = "deg", CUNIT2 = "deg",
                          ...){
  keyvalues = list(
    NAXIS1 = NAXIS1,
    NAXIS2 = NAXIS2,
    CTYPE1 = CTYPE1,
    CTYPE2 = CTYPE2,
    CRVAL1 = CRVAL1,
    CRVAL2 = CRVAL2,
    CRPIX1 = CRPIX1,
    CRPIX2 = CRPIX2,
    CD1_1 = CD1_1,
    CD1_2 = CD1_2,
    CD2_1 = CD2_1,
    CD2_2 = CD2_2,
    RADESYS = 'ICRS',
    EQUINOX = 2000,
    CUNIT1 = CUNIT1,
    CUNIT2 = CUNIT1
  )
  keyvalues = Rwcs_keypass(keyvalues)
  
  #options()$current_keyvalues = keyvalues
  #options()$current_header = NULL
  
  magplot(NA, NA, side=FALSE, xlim=c(0,NAXIS1), ylim=c(0,NAXIS2), asp=1)
  box()
  Rwcs_image(keyvalues = keyvalues, ...)
}

Rwcs_setkeyvalues = function(CRVAL1=0, CRVAL2=0,
                             pixscale=1,
                             NAXIS1=1000, NAXIS2=1000,
                             CRPIX1 = NAXIS1/2 + 0.5, CRPIX2 = NAXIS2/2 + 0.5, # +0.5 for FITS standard
                             CTYPE1='RA---TAN', CTYPE2='DEC--TAN',
                             CD1_1=-pixscale/3600, CD1_2=0, CD2_1=0, CD2_2=pixscale/3600,
                             CUNIT1 = "deg", CUNIT2 = "deg"){
  keyvalues = list(
    NAXIS1 = NAXIS1,
    NAXIS2 = NAXIS2,
    CTYPE1 = CTYPE1,
    CTYPE2 = CTYPE2,
    CRVAL1 = CRVAL1,
    CRVAL2 = CRVAL2,
    CRPIX1 = CRPIX1,
    CRPIX2 = CRPIX2,
    CD1_1 = CD1_1,
    CD1_2 = CD1_2,
    CD2_1 = CD2_1,
    CD2_2 = CD2_2,
    RADESYS = 'ICRS',
    EQUINOX = 2000,
    CUNIT1 = CUNIT1,
    CUNIT2 = CUNIT1
  )
  
  return(keyvalues)
}
