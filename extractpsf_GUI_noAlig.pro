; Reproduce a StarFinder XPsf_Extract GUI PSF extraction from the command line.
; The clicked-star coordinates below are for G028.20-00.05, epoch 1, F160W.

PRO EXTRACTPSF_GUI_noAlig, zone, filter, epoch

  
  ; Command-line defaults for this GUI extraction.
;   zone = 'G028.20-00.05'
;   filter = '160w'
;   filter = '110w'
;   filter = '128n'
;   filter = '164n'
;   epoch = '1'
epochs = ['1', '2']
; epochs = ['1']
filters = ['160w', '110w', '128n', '164n']

; filters = ['160w', '110w']
; filters = ['128n', '164n']
; filters = ['128n']
; filters = ['160w']

zones = [ 'AFGL5180','G028.20-00.05', 'G032.03+00.05', 'G35.2-0.74N', 'G339.88-01.26','IRAS07299-1651', 'IRAS16562-3959']
; zones = [ 'G339.88-01.26']


FOR z = 0, N_ELEMENTS(zones)-1 DO BEGIN
   zone = zones[z]

    For flt = 0, N_ELEMENTS(filters)-1 DO BEGIN
    
        filter = filters[flt]
        
        FOR i = 0, N_ELEMENTS(epochs)-1 DO BEGIN
           epoch = epochs[i]
        
        
          base = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/' + $
                 'SOMA_HST_pms_variability/' + zone + '/epoch' + epoch + '/'
          pattern = base + 'hst_*' + filter + '*'
          paths = FILE_SEARCH(pattern, /TEST_DIRECTORY)
        
          IF (N_ELEMENTS(paths) NE 1) OR (paths[0] EQ '') THEN $
            MESSAGE, 'Expected one image directory matching: ' + pattern
        
          path = paths[0] + '/'
          nam = FILE_BASENAME(paths[0]) + '_drz'
          filename = path + nam + '.fits'
        
          tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/' + $
                   zone + '/F' + filter + '/epoch' + epoch + '/tmp/'
          IF NOT FILE_TEST(tmpdir, /DIRECTORY) THEN FILE_MKDIR, tmpdir
        
          ; HST drizzled image: SCI is extension 1, WHT is extension 2.
          im = READFITS(filename, header_sci, EXT=1)
          noise = READFITS(filename, header_wht, EXT=2)
          
          good = where(FINITE(im),complement=isnan)
          im[isnan] = 0
          noise = sqrt(im)
          ; NOTE: for this product WHT is an effective-exposure/weight map, not a
          ; per-pixel noise standard-deviation map.  PSF_EXTRACT does not use it;
          ; it is read and saved here as requested for downstream use.
        
          bad = WHERE(FINITE(im) EQ 0, n_bad)
          IF n_bad GT 0 THEN im[bad] = 0.0
        
          ; Exact coordinates printed by the successful XPsf_Extract GUI session.
        ;   x_psf = [464L, 578L, 1146L, 298L, 557L, 452L, 672L, 873L, 464L, 564L]
        ;   y_psf = [814L, 901L, 750L, 712L, 490L, 191L, 649L, 795L, 411L, 1104L]
           
;           threshold = 5. * median(noise[where(noise gt 0)])
;           threshold = .1 * median(noise[where(noise gt 0)])

          
          back_box = 21
          background = estimate_background(im,back_box)
          resid = im - background
          sigma = stddev(resid)
          threshold = 3.0*sigma
          background = estimate_background(im,back_box)
          search_objects, im, LOW_SURFACE = background, threshold, $
                          PRE_SMOOTH = 1, MINIF = 2, $ ;THIS WAS CHANGED PRE_SMOOTH AND MINIF. DEFAULTS WERE 1 AND 2.
                          n, x, y, f
          good = where(f gt 0,n)
          x_psf = x[good]
          y_psf = y[good]
          f_psf = f[good]
;           ZP = 25
          ZP = 25
          mag = ZP - 2.5 * alog10(f_psf)
          ord = sort(mag)
          mag = mag[ord]
          x_psf = x_psf[ord]
          y_psf = y_psf[ord]
          f_psf = f_psf[ord]
          mag_psf = mag[ord]
          
          print,N_ELEMENTS(x_psf)
          delta_mag = 5
          delta_r = 10
          ISOLATED_STARS, x_psf, y_psf, mag_psf, x_psf, y_psf, mag_psf, delta_mag, delta_r, ind_iso
           
          x_psf = x_psf[ind_iso]
          y_psf = y_psf[ind_iso]
          f_psf = f_psf[ind_iso]
         
          sz = size(im)
          n1 = sz[1]
          n2 = sz[2]
          psf_fa = psf_gaussian(NPIXEL=41,FWHM=3,/NORMALIZE,/DOUBLE)
        ;   refim = image_model(x,y,f,n1,n2,psf_fa)
        ;   writefits, tmpdir + 'psfstars_iso.fits', refim
          
        ;    x_psf = x_psf[0:30]
        ;    y_psf = y_psf[0:30]
        ;    f_psf = f_psf[0:30]
          
          print, N_ELEMENTS(x_psf)
          
          refim = image_model(x_psf,y_psf,f_psf,n1,n2,psf_fa)
          writefits, tmpdir + 'psfstars.fits', refim
        ;   stop
          
          
          
          ; GUI values reported in the terminal output.
          maskrad = 40
          psf_size = maskrad*2 + 1
          n_fwhm_back = 9.0
          n_fwhm_fit = 2.0
        
          ; No secondary sources were clicked in this GUI extraction.
          ; !NULL is the GUI-compatible representation of an empty selection.
          x_secondary = !NULL
          y_secondary = !NULL
          psf = !NULL
          psf_fwhm = !NULL
          background = !NULL
        
          ; This is the same underlying routine the XPsf_Extract widget calls.
          PSF_EXTRACT, x_psf, y_psf, x_secondary, y_secondary, im, $
                       psf_size, psf, psf_fwhm, background, $
                       N_FWHM_BACK=n_fwhm_back, N_FWHM_FIT=n_fwhm_fit, $
                       INTERP_TYPE='I', /UNWEIGHTED, /NORM_MAX, $
                       N_FWHM_MATCH=1, N_WIDTH=3, MAG_FAC=2
        
          IF TOTAL(psf) LE 0 THEN MESSAGE, 'PSF_EXTRACT returned a non-positive PSF.'
           mid =  n1/2
           print,mid
           
           ; Optional GUI-equivalent halo smoothing
        ; do_halo_smooth = 1
        do_halo_smooth = 0
        
        IF do_halo_smooth THEN BEGIN
           sz_psf = size(psf)
           min_psf_size = min([sz_psf[1], sz_psf[2]])
        
           ; Same defaults used by XPsf_Smooth:
           ; for a 41x41 PSF: r0=10, radial width=5, angular width=22.5 deg.
           r0 = round(float(min_psf_size)/4.0)
           r_width = round(float(r0)/2.0)
           a_width = 22.5 * !pi / 180.0
        
           psf = halo_smooth(psf, r0, $
                             R_WIDTH=r_width, A_WIDTH=a_width, $
                             R_EXP=2, A_EXP=3)
        
           ; Restore unit-flux normalization after smoothing.
        ;    psf = psf / total(psf)
        ENDIF
        
           psf = psf/total(psf)  ; normalization of PSF
          
          WRITEFITS, tmpdir + 'psf_gui_' + filter + '.fits', psf
          WRITEFITS, tmpdir + 'background_gui_' + filter + '.fits', background
          WRITEFITS, tmpdir + 'noise_wht_' + filter + '.fits', noise
        
        ;   PRINT, 'PSF written to: ', tmpdir + 'psf_gui_' + filter + '.fits'
          PRINT, 'PSF written to: ', tmpdir + 'psf_' + filter + '.fits'
          PRINT, 'Measured PSF FWHM (pixels): ', psf_fwhm
        
    endfor ; this is for the epochs
  endfor ; this is for the filters
endfor ; this is for the zone

END
