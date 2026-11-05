PRO ASTROPHOT_orig

; ============================================================
; Loop version of ASTROPHOT_orig.
; Uses the PSFs produced by EXTRACTPSF_orig (loop version) to
; build the StarFinder star lists for filters f160w and f110w
; and epochs 1 and 2.
;
; The HST visit directories are resolved by wildcard:
;   .../SOMA_HST_pms_variability/<zone>/epoch<epoch>/hst_*f<filter>*/
; The drizzled image is <visit>/<visit>_drz.fits (SCI = extension 1)
; and the PSF is read from
;   <visit>/psf_f<filter>_ep<epoch>.fits
; as written by EXTRACTPSF_orig.
;
; The point-source map is built with the extracted PSF
; (not with a Gaussian model), and no FITS file is compressed.
; ============================================================

zone = 'G028.20-00.05'
epochs = ['1', '2']
; epochs = ['1']
filters = ['160w', '110w']
; filters = ['160w']

; Settings for StarFinder you are most likely to 
; want to play with
; these parameters apply to the final StarFinder run
; not to PSF extraction
; Note that the thesholds are set very low, to 0.5 sigma (lower is
; probably even possible)
; This means that our simple noise map overestiamtes the actual uncertainties
; ------------------------------------------------

sf_thresh = [.5,.5]
back_box = 10
deblend = 1
deblost = 0
compbg = 1
posbg = 1

; Apart from the sigma threshold (sf_thresh), the correlation threshold is the most important parameter
; for values close to 1 you may lose real point sources
; for small values (~0.7) you will pick up spurious sources in the
; regionwhere the diffuse emission is high
; With 0.9 you still pick up sources in the region of diffuse emission
min_correlation = 0.9

; General settings for StarFinder
; --------------------------
correl_mag = 4
niter = 2
rel_thresh = 1
guide_x = ""
guide_y = ""

FOR flt = 0, N_ELEMENTS(filters)-1 DO BEGIN

    filter = filters[flt]

    FOR i = 0, N_ELEMENTS(epochs)-1 DO BEGIN
       epoch = epochs[i]
        background = !NULL
       ; --- Resolve the visit directory, the image, and the PSF ---
       base = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/' + zone + '/epoch' + epoch + '/'
       dirs = FILE_SEARCH(base + 'hst_*f' + filter + '*', /TEST_DIRECTORY, COUNT = ndir)

       IF ndir EQ 0 THEN BEGIN
           PRINT, 'ASTROPHOT_orig: no directory found for filter f' + filter + ' epoch ' + epoch
           PRINT, 'ASTROPHOT_orig: pattern was ' + base + 'hst_*f' + filter + '*'
           CONTINUE
       ENDIF

       IF ndir GT 1 THEN BEGIN
           PRINT, 'ASTROPHOT_orig: WARNING, ' + STRTRIM(ndir,2) + ' directories match filter f' + filter + ' epoch ' + epoch + '; using the first:'
           FOR k = 0, ndir-1 DO PRINT, '    ' + dirs[k]
       ENDIF

       path = dirs[0] + '/'
       visit = FILE_BASENAME(dirs[0])
       imfile = path + visit + '_drz.fits'
       psffile = path + 'psf_f' + filter + '_ep' + epoch + '.fits'

       IF NOT FILE_TEST(imfile) THEN BEGIN
           PRINT, 'ASTROPHOT_orig: file not found: ' + imfile
           CONTINUE
       ENDIF

       IF NOT FILE_TEST(psffile) THEN BEGIN
           PRINT, 'ASTROPHOT_orig: PSF not found: ' + psffile
           PRINT, 'ASTROPHOT_orig: run EXTRACTPSF_orig first for f' + filter + ' epoch ' + epoch
           CONTINUE
       ENDIF

       ; create tmp directory if necessary
       tmpdir = path + 'tmp/'
       IF NOT FILE_TEST(tmpdir, /DIRECTORY) THEN FILE_MKDIR, tmpdir

       PRINT, 'ASTROPHOT_orig: epoch ' + epoch + ' filter f' + filter + ' -- image: ' + imfile

       psf = readfits(psffile) ; PSF produced by EXTRACTPSF_orig for this epoch and filter

       ; load image and noise map
       im = readfits(imfile, header, EXT=1) ;careful here! SCI is extension 1 in these drizzled products
       good = where(FINITE(im),complement=isnan)
       im[isnan] = 0
       ; --- Noise map ---
       ; noise = sqrt(im) is NOT the per-pixel noise of a drizzled
       ; count-rate image, and it is exactly zero in the NaN border
       ; regions of the drizzle footprint. With sf_thresh = 0.5 sigma,
       ; a zero noise map makes the detection threshold vanish there,
       ; which produced >16000 spurious suspected stars and the crash
       ; "Illegal subscript range: BACKGROUND" inside StarFinder
       ; (starfinder.pro, line 1133) for epoch 2.
       ; Use the exposure-corrected Poisson noise (same scheme as
       ; EXTRACTPSF_GUI_noAlig) and floor it to a positive value.
       im_ext0 = readfits(imfile, header0, EXT=0)
       exp = sxpar(header0, 'EXPTIME')
       IF (N_ELEMENTS(exp) EQ 1) AND (exp GT 0) THEN BEGIN
           noise = sqrt(abs(im*exp))/exp
       ENDIF ELSE BEGIN
           noise = sqrt(abs(im))
       ENDELSE
       zero = where(noise LE 0.0, nzero)
       IF nzero GT 0 THEN BEGIN
           noise[zero] = median(noise[where(noise GT 0)])
       ENDIF
       writefits, tmpdir + 'im.fits', im
       writefits, tmpdir + 'noise.fits', noise
       sz = size(im)
       n1 = sz[1]
       n2 = sz[2]

; StarFinder run
;####################

       ; Per-combination detection threshold. If the console reports
       ; an excessive number of suspected stars (thousands), raise the
       ; threshold for that combination, e.g.:
       ; IF (epoch EQ '2') AND (filter EQ '160w') THEN sf_thresh_cur = [3.,3.]
       sf_thresh_cur = sf_thresh

       starfinder, im, psf, X_BAD=xbad, Y_BAD = ybad, $
             BACKGROUND = background, BACK_BOX = back_box, $
             sf_thresh_cur, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
             NOISE_STD = noise, min_correlation, $
             CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
             ESTIMATE_BG = compbg, DEBLEND = deblend, DEBLOST = deblost, $
             N_ITER = niter, SILENT=0, $
             GUIDE_X = guide_x, GUIDE_Y = guide_y, $
             SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
             x, y, f, sx, sy, sf, c, STARS = stars, $
             LOGFILE = logfilename, /CUBIC, $
             XBAD = xbad, YBAD = ybad, POSBG = posbg

       IF N_ELEMENTS(f) EQ 0 THEN BEGIN
           PRINT, 'ASTROPHOT_orig: StarFinder found no sources for f' + filter + ' epoch ' + epoch + '; skipping.'
           CONTINUE
           STOP
       ENDIF

       ; --- Output base name for this epoch and filter ---
       basenam = 'hst_' + zone + '_epoch' + epoch

       writefits, path + basenam + '_stars_f'+filter+'.fits', stars, header
       writefits, path + basenam + '_bg_f'+filter+'.fits', background, header
       writefits, path + basenam + '_resid_f'+filter+'.fits', im-stars-background, header
       subtracted = im-stars
       writefits, path + basenam + '_subtracted_f'+filter+'.fits', subtracted, header

       ; save list
       listfile = path + basenam + '_stars_f'+filter+'.txt'
       nstars = n_elements(f)
       openw, outp, listfile, /get_lun
       printf, outp, 'x y f sx sy sf'
       for s = 0, nstars-1 do begin
        xi = round(x[s]) & yi = round(y[s])
        printf, outp, format='(6f13.3)', x[s], y[s], f[s], sx[s], sy[s], sf[s]
       endfor
       free_lun, outp

; make map of point sources with the extracted PSF
; (previously a Gaussian model was used; the map is now
; consistent with the PSF used in the StarFinder fit)
; -------------------------------------------------------------
       readcol, listfile, x, y, f
       im = image_model(x, y, f, n1, n2, psf)
       writefits, tmpdir + basenam + '_map'+filter+'.fits', im

       print, 'Finished ' + basenam + ' filter f' + filter + '.'

    endfor ; this is for the epochs
endfor ; this is for the filters
END