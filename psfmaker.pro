xPRO CALEVAL, X, P, YMOD

  YMOD = P[0]*X
  YMOD = image_shift(YMOD,P[1],P[2])

END

;---------------------

; Contrary to StarFinder PSF_REPEAT_EXTRACT
; this routine takes into account that the references
; stars can also be secondary sources that must be subtracted.
; This is helpful in very dense fields when the reference stars are
; closely spaced.
;
; KEYWORDS
;
; MINDIST    Secondary stars at less than MINDIST separation of a PSF
;            reference star will not be subtracted. This is an
;            untested feature, but I suspect that it avoids rarely
;            occurring problems with PSF extraction in very crowded
;            fields.
;            The default value of MINDIST is 2 pixels.
; LAST MODIFICATION
;
; Introduced MINDIST keyword to avoid that stars very close to the
; reference star (or the reference star itself) can be subtracted as secondaries.
; 
; Activated the CALEVAL procedure again that scales the photometry of
; secondaries to the current best-estimate of the PSF.
; Rainer Schoedel 4 March 2015
;
; Changed nme from MYPSF_REPEAT_EXTRACT to MAKEPSF
; Changed normalisation and weighting of sources.
; Can also be used for first estimation of PSF
; Removed determination of constant background from CALEVAL 
; because it can create problems with nearby unresolved 
; sources or nearby sources that are not included in the list of
; secondaries
; Included potential use of mask
; PSF must be fully covered (important in case of unregularly masked
; image
; Note that iterations appear to make the result slightly worse
;
; 24 July 2024, Rainer Schoedel
; USE_CENTROID can be set to sub-pixel shift the reference stars with the 
; StarFinder CENTROID function. This is necessary if the positions
; of the reference stars have not yet been determined with high
; accuracy.
; 
; 26 MAY 2025, Rainer Schoedel
; 
; Made sure no crash occurs when the reference stars is not found
; during cleaning of secondaries.
; IN that (rare) case, the frame corresponding to this reference star
; is deleted (weight = 0)

PRO PSFMAKER, x_ref, y_ref, x_stars, y_stars, f_stars, image, noise, nrad, FOVMASK = thisfov, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, THRESHOLD = Threshold, MINDIST = mindist, NOISE_PSF = psf_sigma, MASKRAD = mrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=local_sky, USE_CENTROID=use_centroid, oversamp = oversamp


;  print, 'HERE'
 if not(keyword_set(Threshold)) then Threshold = [5.]
 if not(keyword_set(debug)) then debug = 0
 if not(keyword_set(oversamp)) then oversamp = 1
 if not(keyword_set(mindist)) then mindist = 2.
 if keyword_set(mask) then withmask = 1 else withmask = 0
 if (not KEYWORD_SET(satlevel)) then satlevel = 1.e9
 if (not KEYWORD_SET(local_sky)) then local_sky = 0
 if (not KEYWORD_SET(thisfov)) then begin
    thisfov = image
    thisfov[*,*] = 1
 endif
 if not(keyword_set(use_centroid)) then use_centroid = 0

 if (oversamp gt 1) then begin
    sz = size(psf)
    xax_psf = sz[1]
    yax_psf = sz[2]
    psf = CREBIN(psf,sz[1]*oversamp,sz[2]*oversamp,/TOTAL)
    sz = size(image)
    image = CREBIN(image,sz[1]*oversamp,sz[1]*oversamp)
    x_ref = x_ref * oversamp
    y_ref = y_ref * oversamp
    x_stars = x_stars * oversamp
    y_stars = y_stars * oversamp
    thisfov = rebin(thisfov,sz[1]*oversamp,sz[2]*oversamp,/SAMPLE)
    mindist = mindist*oversamp
    nrad = nrad * oversamp
    mrad = mrad * oversamp
 endif
 back_box = mrad
 
 sz = size(psf)
 boxsize = sz[1]
 boxhw = boxsize/2 ; careful: MUST NOT BE FLOAT!
 
 sz = size(image)
 nax1 = sz[1]
 nax2= sz[2]

 nref = n_elements(x_ref)
 xint = round(x_ref)
 yint = round(y_ref)

    
bgring = replicate(1,boxsize,boxsize)
bgring = CIRC_MASK(bgring,boxhw,boxhw,mrad,/INNER)
bgind = where(bgring gt 0, n_bgind)

dummy = replicate(1,boxsize,boxsize)
dummy = CIRC_MASK(bgring,boxhw,boxhw,mrad)
refind = where(dummy gt 0, npix_ref)

if KEYWORD_SET(background) then begin
  sz = size(psfim)
  if (oversamp gt 1) then background = rebin(background,sz[1],sz[2])
  sub_arrays, background, xint, yint, boxsize, bg_stack, masks
endif
 
; 1) Select stars in FoV
;use psfmasks if a reference star may not be contained in all
; cubes because of dithering
; ---------------------------
  sub_arrays, noise, xint, yint, boxsize, noise_stack, masks
  sub_arrays, image, xint, yint, boxsize, stack, masks
  if withmask then begin
    sub_arrays, thisfov, xint, yint, boxsize, psfmasks, masks
  endif
  if debug then begin
    writefits, tmpdir + 'rawstack.fits', stack
;    writefits, 'masks.fits', masks
  endif
  valid_ref = replicate(1,nref)
  for iref = 0, nref-1 do begin
   if withmask then masks[*,*,iref] = masks[*,*,iref]*psfmasks[*,*,iref]
   dummy = circ_mask(masks[*,*,iref],boxhw,boxhw,mrad)
   if total(dummy) lt (0.99*npix_ref) then begin
     valid_ref[iref] = 0
   endif
  endfor
  acceptref = where(valid_ref gt 0, nref_accept)
  xint_accept = xint[acceptref]
  yint_accept = yint[acceptref]
  x_psf_accept = x_ref[acceptref]
  y_psf_accept = y_ref[acceptref]
  stack = stack[*,*,acceptref]
  masks = masks[*,*,acceptref]
 
  print, 'Using: ' + strn(nref_accept) + ' reference stars. '
  psf_weights = fltarr(nref_accept)
  psfnorms = fltarr(nref_accept)

   ; save current stack in rawstack
  rawstack = stack
  subtracted_stack = stack
  
  ; set StarFinder parameters
  estim_bg = 0
  rel_thresh = 1
  min_correlation = 0.7
  deblend = 0
  deblost = 0
  backbox = 0
  niter = 2
  guide_x = ""     ; completely irrelevant, but necessary for our StarFinder version
  guide_y = ""     ; completely irrelevant, but necessary for our StarFinder version
  correl_mag = 4.0 ; as far as I remember this parameter is not very important

 
; (1) Use current PSF and information on image
; to improve PSF estimate
; -----------------------------------------
if (iter gt 0) then begin

  psf = psf/total(psf)

 ; improval of PSF with known sources
  ; ####################################
  for it = 0, iter-1 do begin

    ; clean surroundings of reference stars from secondary sources
    ; ------------------------------------------------------------

     for iref = 0, nref_accept-1 do begin
;      loc_bg = bg_stack[*,*,iref]
;      loc_bg = estimate_background(subtracted_stack[*,*,iref],nrad,/POSBG)
;      loc_bg = estimate_background(subtracted_stack[*,*,iref],nrad,/POSBG,SKY_MEDIAN=1)
      slice = rawstack[*,*,iref]
      slice_noise = noise_stack[*,*,iref]
      submask = masks[*,*,iref]
      ; refit stars and background
      starfinder, slice, psf, X_BAD=x_bad, Y_BAD = y_bad, $
            BACKGROUND = loc_bg, BACK_BOX = back_box, $
                     threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
            NOISE_STD = slice_noise, min_correlation, $
            CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
            ESTIMATE_BG = estim_bg, DEBLEND = deblend, DEBLOST = deblost, $
            N_ITER = niter, SILENT=1, $
            GUIDE_X = guide_x, GUIDE_Y = guide_y, $
            SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
            thisx, thisy, thisf, sx, sy, sf, c, STARS = stars, $
            LOGFILE = logfilename, /CUBIC, /POSBG;, /SKY_MEDIAN
 
      print, 'Found ' + strn(n_elements(thisf)) +  ' stars in slice ' + strn(iref) +  '.'

      ref_ind = -1
      compare_lists, thisx, thisy, x_psf_accept[iref]-(xint_accept[iref]-boxhw), y_psf_accept[iref]-(yint_accept[iref]-boxhw), x1c, y1c, x2c ,y2c,  SUBSCRIPTS_1 = ref_ind, SUB1=other_stars, MAX_DISTANCE=mindist
      if (other_stars[0] gt -1) then begin
         secondaries = image_model(thisx[other_stars],thisy[other_stars],thisf[other_stars],boxsize,boxsize,psf,REFERENCE_PIX=[boxhw,boxhw])
          
          
         slice = slice - secondaries - loc_bg
       endif                     ;  (other_stars[0] gt -1)


                                ; skip the following section if the
                                ; reference star has not been
                                ; detectedm for some reason:
      if (ref_ind gt (-1)) then begin
         ; if more than one star is detected at
         ; the reference position, then choose
         ; the brightest one
         reflen = n_elements(ref_ind)
         if (reflen gt 1) then begin
            ord = REVERSE(sort(thisf[ref_ind]))
            ref_ind = ref_ind[ord[0]]
         endif
         xs = boxhw - thisx[ref_ind]
         ys = boxhw - thisy[ref_ind]
         slice = image_shift(slice,xs,ys)
         bg_stack[*,*,iref] = loc_bg
         subtracted_stack[*,*,iref] = (rawstack[*,*,iref] - stars - loc_bg)/thisf[ref_ind]
         stack[*,*,iref] = slice/thisf[ref_ind]
         submask = image_shift(submask,xs,ys)
         psf_weights[iref] = sqrt(thisf[ref_ind])
      endif else begin
         subtracted_stack[*,*,iref] = rawstack[*,*,iref]/thisf[ref_ind]
         stack[*,*,iref] = slice/thisf[ref_ind]
         submask = image_shift(submask,xs,ys)
         psf_weights[iref] = 0
      endelse
    endfor ; end loop over reference stars

    
     psf_weights = round(psf_weights/min(psf_weights))
     if (unweighted eq 1) then psf_weights[*] = 1
     print, 'Weights: '
     print, psf_weights
     psf = stack_median(stack, WEIGHTS=psf_weights, MASK=masks)
     psf_sigma = stack_error(stack, WEIGHTS=psf_weights, MASK=masks)

     ; correct background and make circular mask
     mmm, psf, skymod, skysigma , skyskew
     psf = psf - skymod
     neg = where(psf lt 0)
     psf[neg] = 0
       
    if (debug gt 0) then begin
       writefits, tmpdir + 'stack_' +strtrim(string(it+1),2) + '.fits', stack
       writefits, tmpdir + 'rawstack.fits', rawstack
       writefits, tmpdir + 'subtracted_stack_' +strtrim(string(it+1),2) + '.fits', subtracted_stack
       writefits, tmpdir + 'subtracted_median_' +strtrim(string(it+1),2) + '.fits', MEDIAN(subtracted_stack,DIM=3)
       writefits, tmpdir + 'bg_stack_' +strtrim(string(it+1),2) + '.fits', bg_stack
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'raw.fits', psf
       writefits, tmpdir + 'psf_' +strtrim(string(it+1),2) + 'sigma.fits', psf_sigma
;       STOP
      endif

     psf = circ_mask(psf, boxhw, boxhw, mrad)
     psf = psf/total(psf)  ; normalization of PSF

  endfor                      ; end loop over iter

; if PSF does not yet exist, extract one
; ------------------------------------------
endif else begin

  psf_weights = intarr(nref_accept) ; for optional weighting
  for iref = 0, nref_accept-1 do begin
    subim = stack[*,*,iref]
    submask = masks[*,*,iref]
    ; estimate local sky
    if local_sky then begin
      mmm, subim[bgind], skymod, skysig, skyskew, /SILENT       
      subim = subim - skymod
    endif
    ; sub-pixel shift
    subim = image_shift(subim,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
    submask = image_shift(submask,xint_accept[iref]-x_psf_accept[iref],yint_accept[iref]-y_psf_accept[iref])
    zeroind = where(submask lt 0.99, complement = ones)
    if (use_centroid) then begin
      subim = CENTROIDER(subim,XSHIFT=xs,YSHIFT=ys)
      submask = image_shift(submask,xs,ys)
    endif
    if (zeroind[0] gt -1) then begin
      submask[zeroind] = 0
      submask[ones] = 1
    endif  
    subim = subim * submask
    masks[*,*,iref] = submask

    ; Normalize the reference stars
    ; normalization will underweight stars with saturated cores
    ; should be no problem
    normim = subim
    normim = circ_mask(normim, boxhw, boxhw, nrad)
    fnorm = total(normim)
    saturated = where(subim gt satlevel)
    submask[saturated] = 0
    normim = normim * submask
    subim = subim/fnorm
    psf_weights[iref] = sqrt(fnorm)
    stack[*,*,iref] = subim*submask
    masks[*,*,iref] = submask
   endfor
   psf_weights = round(psf_weights/min(psf_weights))
   if (unweighted eq 1) then psf_weights[*] = 1
   print, 'Weights: '+ strn(psf_weights)
   psf = stack_median(stack, WEIGHTS=psf_weights, MASK=masks)
   psf_sigma = stack_error(stack, WEIGHTS=psf_weights, MASK=masks)

   if (debug gt 0) then begin
     writefits, tmpdir + 'im.fits', psfim
     if withmask then writefits, tmpdir + 'thisfov.fits', thisfov
     writefits, tmpdir + 'bgring.fits', psf*bgring
     writefits, tmpdir + 'psfraw.fits', psf
     writefits, tmpdir + 'rawstack.fits', rawstack
     writefits, tmpdir + 'stackmasks.fits', masks
   endif

   ; circular mask and normalise PSF
   ; -------------------------------
;   MASK_PSF, psf, mrad, PSF_MASKED=psf_masked ; to avoid having negative wings of the PSF
;   psf = psf_masked
;   psf = psf/total(psf)
 
endelse

if (oversamp gt 1) then psf = CREBIN(psf,xax_psf,yax_psf,/TOTAL)


END
