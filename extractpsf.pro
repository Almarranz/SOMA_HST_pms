PRO EXTRACTPSF, zone, filter, epoch

; tmpdir = '../tmp/'
; tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/pruebas/tmp/'
; dir = './'

; zone = 'G028.20-00.05'
; filter = '160w'
; epoch = '1'

path = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'+ zone +'/gaia_alignment/Epoch'+ epoch +'/'
pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/pruebas/'
tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'/epoch'+ epoch +'/tmp/'
results = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'/epoch'+ epoch +'/'


ZP = 25.0 ; random ZP
maskrad = 21
nrad = 4
; path = '/Users/fedriani/Documents/postdoc_iaa/HST_project/HST_data/G028.20-00.05/gaia_alignment/Epoch1/starfinder/'
; path = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/G028.20-00.05/gaia_alignment/Epoch1/'
; results = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/G028.20-00.05/f160w/epoch1/'
; pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/pruebas/'


; nam = zone +'_EP'+ epoch + '_f'+ filter + '_drz_sci'

nam = 'hst_ep'+ epoch + '_f'+ filter + '_drz'
; create tmp directory if necessary
; if not(FILE_TEST(path + 'tmp')) then FILE_MKDIR, path + 'tmp'
if not(FILE_TEST(pruebas + 'tmp')) then FILE_MKDIR, pruebas + 'tmp'


; minimum mag difference and distance for secondary stars near
; reference stars
delta_mag = 5.
delta_r = 10
nref_max = 21 ; max number of PSF reference stars: result not significantly sensitive to this parameter

; Parameters that need to be edited frequently
min_correlation = 0.7 ; high correlation threshold is importante to avoid detecting saturated or corrupted sources
psf_fwhm = 2.0 ; very approximate value
unweighted = 1 ; If 1, then use unweighted median of the stars that are selected to represent the PSF

; Parameters for PSF estimation and StarFinder
; ------------------------------------------------

correl_mag = 4.0 ; as far as I remember this parameter is not very important
deblend = 0     ; deblend close stars?
deblost = 0      ; try to deblend stars that appear to deviate from PSF shape (i.e. they appear to be merged)
niter = 2        ; 2 iterations are standard
rel_thresh = 1   ; use relative threshold for source detection (i.e. sigmas) 
guide_x = ""     ; completely irrelevant, but necessary for our StarFinder version
guide_y = ""     ; completely irrelevant, but necessary for our StarFinder version
back_box = maskrad           ; radius for estimation of background 
psf_size = 2*maskrad+1

;   im = readfits(path + nam + '.fits',header,EXT=0) ;careful here! check the extension!
  im = readfits(path + nam + '.fits',header,EXT=1) ;careful here! check the extension!
  good = where(FINITE(im),complement=isnan)
  im[isnan] = 0
  noise = sqrt(im)
;   writefits, path + 'im'+filter+'.fits', im, /COMPRESS
  writefits, pruebas + 'im'+filter+'.fits', im 
  writefits, tmpdir + 'im.fits', im
  writefits, tmpdir + 'noise.fits', noise
  sz = size(im)
  n1 = sz[1]
  n2 = sz[2]
  support = replicate(0,n1,n2)
  support[good] = 1
  
; 1) First estimate of PSF
; detect PSF reference sources automatically
; ----------------------------
  threshold = 100. * median(noise[where(noise gt 0)])
  background = estimate_background(im,back_box)
  search_objects, im, LOW_SURFACE = background, threshold, $
                  PRE_SMOOTH = 1, MINIF = 2, $ ;THIS WAS CHANGED PRE_SMOOTH AND MINIF. DEFAULTS WERE 1 AND 2.
                  n, x, y, f
  good = where(f gt 0,n)
  x_psf = x[good]
  y_psf = y[good]
  f_psf = f[good]
  print, 'using '+ strn(n) + ' PSF reference stars.'
  debug = 0
  iter = 1                        ; more iterations do not necessarily make it better
  mindist = 10 ; THIS WAS 2, BUT WORKS AT 30. WORKED AT 10 FOR G28
 ; Use Gaussian PSF
  psf = psf_gaussian(NPIXEL=psf_size,FWHM=psf_fwhm,/NORMALIZE,/DOUBLE)
  threshold = 10
  PSFMAKER, x_psf, y_psf, x, y, f, im, noise, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_sigma, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=local_sky, USE_CENTROID=use_centroid, oversamp = oversamp;, THRESHOLD=threshold
;   PSFMAKER, x_psf, y_psf, x, y, f, im, noise, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_sigma, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=local_sky, USE_CENTROID=use_centroid, oversamp = oversamp, THRESHOLD=threshold

 
  mmm, psf, skymod, skysigma , skyskew
  psf = psf - skymod
  neg = where(psf lt 0,count)
  if (count gt 0) then  psf[neg] = 0
  psf = circ_mask(psf, maskrad, maskrad, maskrad)
  psf = psf/total(psf)  ; normalization of PSF
  writefits, tmpdir + 'tmppsf.fits', psf
;STOP
  ; Run StarFinder and iterate search for PSF reference stars and PSF extraction
  ; --------------------------------------------------------------------------

  Threshold = [5, 5] ;THIS WAS 1, WORKS WITH 5.
  estim_bg = 1
  starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
                     threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, DEBLOST = deblost, $
        N_ITER = niter, SILENT=1, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
    writefits, tmpdir + 'stars.fits', stars
    writefits, tmpdir + 'bg.fits', background
    
    sat_th = 15000
    
    unsat = WHERE(f lt sat_th, n_unsat)
    IF n_unsat GT 0 THEN BEGIN
        x = x[unsat]
        y = y[unsat]
        f = f[unsat]
    ENDIF
;     
    ; save stars for use in PSFMAKER
    x_stars = x
    y_stars = y
    f_stars = f
    
  
  
    
 ; 2) PSF extraction
 ; ---------------------

  mag = ZP - 2.5 * alog10(f)
  ord = sort(mag)
  mag = mag[ord]
  x = x[ord]
  y = y[ord]
  x_psf = x
  y_psf = y
 mag_psf = mag
 ISOLATED_STARS, x, y, mag, x_psf, y_psf, mag_psf, delta_mag, delta_r, ind_iso
  x = x[ind_iso]
  y = y[ind_iso]
  f = f[ind_iso]
  mag = mag[ind_iso]
  ord = sort(mag)
  mag = mag[ord]
  x = x[ord]
  y = y[ord]
  f = f[ord]
  n_iso = n_elements(ind_iso)
  print, 'Found ' + strn(n_iso) + ' isolated stars.'
  
  ; Avoid masked sources or sources near edges
  ; Valid reference sources
  ; are those sources that are contained
  ; in field of view (at least to psf_frac part)
  psf_frac = 0.7
  boxhw = psf_size/2
  dummy = replicate(1,psf_size,psf_size)
  dummy = circ_mask(dummy,boxhw,boxhw,maskrad)
  npix_psf = total(dummy)
  valid_ref = replicate(1,n_iso)

  sub_arrays, support, round(x), round(y), psf_size, psffovs, psfmasks
  for iref = 0, n_iso-1 do begin
     tmpmask = psfmasks[*,*,iref]*psffovs[*,*,iref]
     dummy = circ_mask(tmpmask,boxhw,boxhw,maskrad)
     if total(dummy) lt round(psf_frac*npix_psf) then begin
        valid_ref[iref] = 0
     endif
  endfor
  acceptref = where(valid_ref gt 0, nref_accept,complement=reject)
  x = x[acceptref]
  y = y[acceptref]
  f = f[acceptref]
  up = (nref_max < nref_accept) - 1
  
  x_psf = x[0:up]
  y_psf = y[0:up]
  f_psf = f[0:up]

; make map of point sources with Gaussian PSFs (helpful for extremely
; crowded fields)
; -------------------------------------------------------------
  refim = image_model(x_psf,y_psf,f_psf,n1,n2,psf)
  writefits, tmpdir + 'psfstars.fits', refim
  print,'*******************'
  print, 'n_elements(x_psf)', n_elements(x_psf)
  
   iter = 1
   oversamp = 1
   use_centroid = 0
   debug = 0
   PSFMAKER, x_psf, y_psf, x_stars, y_stars, f_stars, im, noise, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_noise, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=1, oversamp = oversamp, USE_CENTROID=use_centroid, THRESHOLD=threshold

  sz = size(psf)
  n1 = sz[1]
  n2 = sz[2]
  mid =  n1/2
;  mmm, psf, skymod, skysigma , skyskew
;  psf = psf - skymod
;  neg = where(psf lt 0)
;  psf[neg] = 0
  psf = circ_mask(psf, mid, mid, maskrad)
  psf = psf/total(psf)  ; normalization of PSF
;   writefits, path + 'psf_'+filter+'.fits', psf
;   writefits, path + 'psf_sigma'+filter+'.fits', psf_noise/total(psf)
  
  writefits, tmpdir + 'psf_'+filter+'.fits', psf
  writefits, tmpdir + 'psf_sigma'+filter+'.fits', psf_noise/total(psf)
  
  print, 'FINITO'
  

END
