PRO EXTRACTPSF_GUI_noAlig, zone, filter, epoch


zone = 'G028.20-00.05'
filter = '160w'
epoch = '1'


; base = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/G028.20-00.05/epoch1/'
base = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'+ zone +'/epoch'+ epoch +'/'




pattern = base + 'hst_*' + filter + '*'
paths = FILE_SEARCH(pattern, /TEST_DIRECTORY)
path = paths + '/'

pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/pruebas/'
tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'_noAlig/epoch'+ epoch +'/tmp/'
results = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'_noAlig/epoch'+ epoch +'/'


ZP = 25.0 ; random ZP

; Keep the scripted PSF model identical to the GUI extraction below.
; XPsf_Extract reported PSF_SIZE = 41 for the successful GUI run.
psf_size = 83L
maskrad = (psf_size - 1L)/2L
nrad = 4

; Set to 1 only after confirming that the GUI-equivalent PSF below is good.
; Leave it at 0 to make this procedure reproduce the GUI extraction directly.
run_starfinder_refinement = 0




nam = FILE_BASENAME(path) + '_drz
print, path
print, nam
; nam = 'hst_ep'+ epoch + '_f'+ filter + '_drz'
; create tmp directory if necessary
; if not(FILE_TEST(path + 'tmp')) then FILE_MKDIR, path + 'tmp'
if not(FILE_TEST(pruebas + 'tmp')) then FILE_MKDIR, pruebas + 'tmp'


; minimum mag difference and distance for secondary stars near
; reference stars
delta_mag = 5.
delta_r = maskrad
nref_max = 10 ; max number of PSF reference stars: result not significantly sensitive to this parameter

; Parameters that need to be edited frequently
min_correlation = 0.9 ; high correlation threshold is importante to avoid detecting saturated or corrupted sources
psf_fwhm = 1.3 ; very approximate value
unweighted = 1 ; If 1, then use unweighted median of the stars that are selected to represent the PSF

; Parameters for PSF estimation and StarFinder
; ------------------------------------------------

correl_mag = 4.0 ; as far as I remember this parameter is not very important
deblend = 0     ; deblend close stars?
deblost = 0      ; try to deblend stars that appear to deviate from PSF shape (i.e. they appear to be merged)
niter = 1        ; 2 iterations are standard
rel_thresh = 1   ; use relative threshold for source detection (i.e. sigmas) 
guide_x = ""     ; completely irrelevant, but necessary for our StarFinder version
guide_y = ""     ; completely irrelevant, but necessary for our StarFinder version
back_box = maskrad           ; radius for estimation of background 

;   im = readfits(path + nam + '.fits',header,EXT=0) ;careful here! check the extension!
  im = readfits(path + nam + '.fits',header,EXT=1) ;careful here! check the extension!
  sz = size(im)
  n1 = sz[1]
  n2 = sz[2]
  
  good = where(FINITE(im),complement=isnan)
  im[isnan] = 0
;   noise = sqrt(im)
; NOISE estimation
    ; SCI is in electrons/s; WHT is effective exposure time in seconds
;       noise = readfits(path + nam + '.fits', header_wht, EXT=2)
    wht = readfits(path + nam + '.fits', header_wht, EXT=2)
;     
;     ; Two input FLT images were drizzled. Header read noise is ~20 e-.
    readnoise = 20.0D
    nexp = 2.0D
;     
    noise = replicate(1.0D, n1, n2)  ; safe initial value
    ok = where(finite(im) AND finite(wht) AND (wht GT 0), n_ok)
;     
    IF n_ok GT 0 THEN BEGIN
       ; Poisson + read-noise uncertainty, expressed in electrons/s.
       ; The max prevents negative sky-subtracted values creating NaNs.
       rate = im[ok] > 0.0D
       noise[ok] = sqrt(rate/wht[ok] + nexp*(readnoise/wht[ok])^2)
   ENDIF
 
;   writefits, path + 'im'+filter+'.fits', im, /COMPRESS
  writefits, pruebas + 'im'+filter+'.fits', im 
  writefits, tmpdir + 'im.fits', im
  writefits, tmpdir + 'noise.fits', noise
 
  support = replicate(0,n1,n2)
  support[good] = 1
  
; 1) GUI-equivalent PSF extraction
; --------------------------------
; These are the ten coordinates printed by XPsf_Extract in the GUI session.
; PSF_EXTRACT performs the GUI's centroiding, median stack, neighbour fitting,
; and saturated-core handling.  Do not replace this result with a Gaussian or
; with PSFMAKER before comparing it with the GUI result.

; x_psf = [464L, 578L, 1146L, 298L, 557L, 452L, 672L, 873L, 464L, 564L]
; y_psf = [814L, 901L, 750L, 712L, 490L, 191L, 649L, 795L, 411L, 1104L]

; Automatic selection of clean PSF stars
seed_sigma = 200.0D
n_auto_max = 40L

background = estimate_background(im, back_box)
seed_threshold = seed_sigma * median(noise[where(noise GT 0)])

search_objects, im, LOW_SURFACE=background, seed_threshold, $
                PRE_SMOOTH=1, MINIF=2, n, x, y, f

; Reject non-positive detections and stars whose full PSF box crosses an edge.
good = where((f GT 0) AND $
             (x GE maskrad) AND (x LT n1-maskrad) AND $
             (y GE maskrad) AND (y LT n2-maskrad), n_good)

IF n_good LT 4 THEN MESSAGE, 'Too few bright PSF candidates found.'

x_cand = x[good]
y_cand = y[good]
f_cand = f[good]
mag_cand = -2.5D * alog10(f_cand)

; With no manually selected secondary-source list, demand strong isolation.
delta_mag_auto = 2.0D
delta_r_auto = 1L * maskrad

ISOLATED_STARS, x_cand, y_cand, mag_cand, $
                x_cand, y_cand, mag_cand, $
                delta_mag_auto, delta_r_auto, ind_iso

n_iso = n_elements(ind_iso)
IF n_iso LT 4 THEN MESSAGE, 'Too few isolated PSF candidates found.'

x_iso = x_cand[ind_iso]
y_iso = y_cand[ind_iso]
f_iso = f_cand[ind_iso]

; Keep only the brightest isolated candidates.
ord = reverse(sort(f_iso))
n_use = n_auto_max
IF n_iso LT n_auto_max THEN n_use = n_iso

x_psf = x_iso[ord[0:n_use-1]]
y_psf = y_iso[ord[0:n_use-1]]

print, 'Using ', n_use, ' automatic isolated PSF stars.'
print, x_psf
print, y_psf
x_secondary = !NULL
y_secondary = !NULL
psf_fwhm_gui = 0.0
background = 0.0

psf_extract, x_psf, y_psf, x_secondary, y_secondary, im, $
             psf_size, psf, psf_fwhm_gui, background, $
             N_FWHM_BACK=9.0, N_FWHM_FIT=2.0, INTERP_TYPE='I'

IF total(psf) LE 0 THEN MESSAGE, 'PSF_EXTRACT returned a non-positive PSF.'
psf = psf / total(psf)
writefits, tmpdir + 'tmppsf_gui_equivalent.fits', psf

print, 'Using ', n_use, ' automatic isolated PSF stars.'
print, x_psf
print, y_psf
STOP
___________________________________
zone = 'G028.20-00.05'
filter = '160w'
epoch = '1'


; base = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/G028.20-00.05/epoch1/'
base = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'+ zone +'/epoch'+ epoch +'/'




pattern = base + 'hst_*' + filter + '*'
paths = FILE_SEARCH(pattern, /TEST_DIRECTORY)
path = paths + '/'

pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/pruebas/'
tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'_noAlig/epoch'+ epoch +'/tmp/'
results = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'_noAlig/epoch'+ epoch +'/'


ZP = 25.0 ; random ZP
maskrad = 41
nrad = 4




nam = FILE_BASENAME(path) + '_drz
print, path
print, nam
; nam = 'hst_ep'+ epoch + '_f'+ filter + '_drz'
; create tmp directory if necessary
; if not(FILE_TEST(path + 'tmp')) then FILE_MKDIR, path + 'tmp'
if not(FILE_TEST(pruebas + 'tmp')) then FILE_MKDIR, pruebas + 'tmp'


; minimum mag difference and distance for secondary stars near
; reference stars
delta_mag = 5.
delta_r = 10
nref_max = 20 ; max number of PSF reference stars: result not significantly sensitive to this parameter

; Parameters that need to be edited frequently
min_correlation = 0.9 ; high correlation threshold is importante to avoid detecting saturated or corrupted sources
psf_fwhm = 2.0 ; very approximate value
unweighted = 1 ; If 1, then use unweighted median of the stars that are selected to represent the PSF

; Parameters for PSF estimation and StarFinder
; ------------------------------------------------

correl_mag = 4.0 ; as far as I remember this parameter is not very important
deblend = 1    ; deblend close stars?
deblost = 0      ; try to deblend stars that appear to deviate from PSF shape (i.e. they appear to be merged)
niter = 1        ; 2 iterations are standard
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
  
   ; This lauch the widget to manually select the psf stars
;    device, decomposed = 0
;    boxsize = 11. ; width of box within which maximum is searched after a click
;    disp_opt = default_display_opt(im)
;    disp_opt.stretch = 'logarithm'
;    disp_opt.range = max(im)*[1.0e-3,1.0]
;    display_image, im, wnum, OPTIONS = disp_opt, MODIFY_OPT = modify_opt
;    click_on_max, im, /MARK, BOXSIZE = boxsize, x_0, y_0
;    measure_centroid, im, x_0, y_0, boxsize
;    confirm_stars, im, wnum, disp_opt, x_0, y_0, psf_size, x_psf, y_psf
;    wdelete, wnum
;    
; ;     x_psf = x_0
; ;     y_psf = y_0
; 
;    print, x_psf, y_psf
   
   
  
  
; 1) First estimate of PSF
; detect PSF reference sources automatically
; ----------------------------
  threshold = 500. * median(noise[where(noise gt 0)])
  background = estimate_background(im,back_box)
  search_objects, im, LOW_SURFACE = background, threshold, $
                  PRE_SMOOTH = 1, MINIF = 2 , $ ;THIS WAS CHANGED PRE_SMOOTH AND MINIF. DEFAULTS WERE 1 AND 2.
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
  
  threshold = 100
  x_secondary = !NULL
  y_secondary = !NULL
  psf_fwhm_gui = 0.0
  background = 0.0
  
  psf_extract, x_psf, y_psf, x_secondary, y_secondary, im, $
             psf_size, psf, psf_fwhm_gui, background, $
             N_FWHM_BACK=9.0, N_FWHM_FIT=2.0, INTERP_TYPE='I'
;   PSFMAKER, x_psf, y_psf, x, y, f, im, noise, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_sigma, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=local_sky, USE_CENTROID=use_centroid, oversamp = oversamp, THRESHOLD=threshold

 
  mmm, psf, skymod, skysigma , skyskew
  psf = psf - skymod
  neg = where(psf lt 0,count)
  if (count gt 0) then  psf[neg] = 0
  psf = circ_mask(psf, maskrad, maskrad, maskrad)
  psf = psf/total(psf)  ; normalization of PSF
  writefits, tmpdir + 'tmppsf.fits', psf
  
 
  
  
STOP
;   Run StarFinder and iterate search for PSF reference stars and PSF extraction
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
    
    sat_th = 50000
;     stop
    
    
;     unsat = WHERE(f lt sat_th, n_unsat)
;     IF n_unsat GT 0 THEN BEGIN
;         x = x[unsat]
;         y = y[unsat]
;         f = f[unsat]
;     ENDIF
; ;     
;     ; save stars for use in PSFMAKER
;     x_stars = x
;     y_stars = y
;     f_stars = f
    
  
  
    
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
  psf_frac = 0.9
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