PRO COMPLETENESS, chip, f0, thismag ;no necesito el zp, necesito f0 :)

;chip = namw
;f0 = vega fluxes
;thismag = the mag. you want to calculate the completenes on
file_path = 'names.txt'

; Read the values from the file
; The FORMAT='A' option reads the values as strings
READCOL, file_path, band, field, FORMAT='A,A' ;Original: survey, band, field FORMAT = 'A,A,A' ; Ahora leemos f0, que cogemos el mismo
;para todos los chips de la banda y lo necesito para calcular el flujo. NO, EL F0 LO METO A MANO 
print, band
print, chip
;PATH
;-------------------------------------
basedir = '/home/data/working/lucia/JWST_latest/' + band + '/' + chip + '_' + field + '/' ;Example: F212N/nrca1_p1/
in_path = basedir                                ;basedir + 'cubeims/chip' + strn(chip) + '/'
out_path = in_path + 'Completeness/'
phot_path =  basedir                                    ;basedir + '/photo/chip' + strn(chip) + '/' 
tmpdir =  basedir + 'tmp/'                                      ;basedir + '/tmp/tmp' + strn(chip) + '/'

; Important parameters
; ----------------------
; use settings from ASTROPHOT.PRO
rebfac = 0                   ; original: 2 (GNS). Para JWST no estamos aplicando rebining
psf_fwhm = 2.0 ;* rebfac      ; approximate PSF FWHM
;psf_box = 20                 ; size of PSF without rebinning -> I would say I don't need this either 
;maskrad = round(12 * rebfac) ; circular mask for PSF -> creo que esto no lo necesito porque no tengo que extraer la PSF

; Create a grid of artificial stars that is sufficiently wide to not increase crowding
; to explore completeness as a funciton of position better, this grid is shifted
; around the image in separate steps
; Ahora no necesito mag_min, nmag ni magstep porque al paralelizar paso las magnitudes por fichero
; ==========================================================================================
step = 26L                ; grid step size in pixels CHANGE -> 10xFWHM, entonces sería 20, porque yo no tengo el rebining factor :). Lo hemos subido un poquito para que sea ~0.8"
gridstep = step/4         ; distance by which to shift the grid in 4 separate runs -> KEEP; shifts de ~2"
randstep = 0.5            ; random jitter of artificial sources around grid positions -> KEEP, doesn't matter much
magstep = 0.5             ; completeness will be probed in steps of 0.5 mag -> KEEP in principle, unless it is too slow; quizás es incluso mucho... pero bueno siempre s epuede interpolar a unas malas 
;mag_min = 19              ; brightest magnitude to examine -> SELECT FROM PRELIMINAR (REDUNDANT) COMPLETENESS!!
;nmag = 6                 ; number of magnitudes to examine -> UNTIL REACHING WHERE THE LF FALLS!!! 
gain = 2.0               ; approximate gain for all detectors (precise value should not matter much)

; parameters to define redeteciton of an artificial star --> ESTO POR AHORA LO MANTENGO, A VER Cómo workea 
; ------------------------------------------------------
dmagmax = magstep/2. ; maximum accepted mteagnitude difference between recovered source and input source
dmax = 1.5           ; maximum position offset between recovered source and input source - USE VALUE FROM COMPUTE_UNCERTAINTIES.PRO

; read image and noise map
; map of number of exposures - exp - is needed to compute photon noise for artificial stars
; ------------------------------------------------------
;im = readfits(in_path + 'B1_chip'+strn(chip)+'_holo_cal.fits.gz')
im = readfits(in_path + 'im.fits', EXT=1, header)
;LEEMOS TAMBIÉN EL PIXAR_SR QUE LO NECESITAMOS PARA CALCULAR EL FLUJO!!
PIXAR_SR = SXPAR(header,'PIXAR_SR','keyword not found')
exp = SXPAR(header,'XPOSURE','keyword not found')
sz = size(im)
n1 = sz[1] ;Antes era n1_0 y n2_0, pero como yo no tengo que recortar la imagen pues directamente n1,n2
n2 = sz[2]

;noise = readfits(in_path + 'B1_chip'+strn(chip)+'_noise_cal.fits.gz')
noise = readfits(in_path + 'err.fits')

;exp = readfits(in_path + 'B1_chip'+strn(chip)+'_exp_cal.fits.gz') ;--> NO TENGO EXP MAPS ASÍ QUE LEEMOS XPOSURE DEL HEADER 

; read list of deteted stars: ONLY NECESSARY FOR GNS
; needed to know which areas of the images contain information
;-------------------------------------------------------------

;readcol, phot_path + 'calibrated_' + band + '_' + strn(field) + '_' + strn(chip) + '_opti.txt', a, d, xreal, yreal, freal, sx, sy, sf, expos, ref_flag
;mm = ZP - 2.5*alog10(freal/texp)
;x_min = round(min(xreal))
;x_max = round(max(xreal))
;y_min = round(min(yreal))
;y_max = round(max(yreal))

; ONLY NEEDED FOR GNS -> LO VOY A COMENTAR THEN 
; In case of GNS1 the images for each chip are much larger than the area that contains any actual data
; extract sub-images to speed up the process
; sub-images are a bit larger than the area containing data in otder to make sure completeness is sampled out to the edges
; ------------------------------------------------------
;x_min = (x_min-step) > 0
;y_min = (y_min-step) > 0
;x_max = (x_max+step) < (n1_0-1)
;y_max = (y_max+step) < (n2_0-1)
;im = im[x_min:x_max,y_min:y_max]
;noise = noise[x_min:x_max,y_min:y_max]
;exp = exp[x_min:x_max,y_min:y_max]
;sz = size(im)
;n1 = sz[1]
;n2 = sz[2]
;writefits, tmpdir + 'im_completeness.fits' , im

; Extract PSF -> ESTO TAMPOCO LO NECESITO -> LEER WebbPSF que haya generado directamente!!!!
; 'mask' is needed but not really used
; ------------------------------------------------------
;mask = replicate(1,n1,n2) ; saturated stars not considered!
;psf = psf_completeness(im, noise, mask, maskrad, psf_fwhm, tmpdir, texp, ZP, /UNWEIGHTED)
;psf = psf/total(psf)
;writefits, in_path + 'psf_holo.fits' , psf


;LECTURA WEBBPSF 
psf = readfits(in_path + 'webbpsf.fits')
;STOP



; loop over magnitudes
; ----------------------------------------------------------------------------------

;for ms = 1, nmag do begin

   
;  thismag = mag_min + ms * magstep      ; magnitude probed in this loop
  ;flux = texp *  10^(0.4*(ZP-thismag))  ; corresponding flux in ADU -> este es el de GNS
  flux = (f0/PIXAR_SR)*10^(-thismag/2.5) ;expresión que utilizamos en los códigos de JWST 
  
  recovered_fluxes = flux               ; recovered_fluxes will save fluxes of measured artificial stars
                                        ; to check whether there is any bias of photometry

  inmap = fltarr(n1,n2)                 ; will contain pixels with value 1 for artificial stars injected into the image
  outmap = fltarr(n1,n2)                ; will contain pixels with value 1 for recovered artificial stars

 
  ; make 4 grids of artificial stars
  ; shifted with respect to each other --> I didn't modify anything here 
  ; ---------------------------------------------------- 
  grid = 0L
  for xg = 0, 3 do begin
    for yg = 0, 3 do begin
      grid++
 
   ; create grid of artificial stars --> did not modify anything 
   ;-------------------------
      xgrid = xg*gridstep + step*findgen(n1/step)
      ygrid = yg*gridstep + step*findgen(n2/step)
      nxgrid = n_elements(xgrid)
      nygrid = n_elements(ygrid)
      xin = fltarr(nxgrid*nygrid)
      yin = fltarr(nxgrid*nygrid)
      npos = 0L
      for i = 0, nxgrid-1 do begin
        for j = 0, nygrid-1 do begin
          xin[npos] = xgrid[i]
          yin[npos] = ygrid[j]
          npos++
        endfor
      endfor
        
     ; add random component of +-0.5 pixel to positions --> did not modify anything 
     ;-------------------------------------------------
     xdiff = randstep*randomu(seedx,npos)
     xin = xin + xdiff
     ydiff = randstep*randomu(seedy,npos)
     yin = yin + ydiff
     fin = replicate(flux,npos)

     ; discard positions outside of fov to avoid errors --> did not modify anything 
     ;----------------------------------
     ind = where((xin ge 0) and (xin le (n1-1)) and (yin ge 0) and (yin le (n2-1)))
     xin = xin[ind]
     yin = yin[ind]
     fin = fin[ind]

    ; create image of artificial stars     --> did not modify anything 
    ;------------------------------------
     print, 'Now creating image of artificial stars...' 
     sources = image_model(xin,yin,fin,n1,n2,psf)
     inmap[round(xin),round(yin)] = 1
;     writefits, tmpdir + 'sources.fits', sources    

    ; add POISSON noise to image of artificial stars
    ; save POISSON noise for each pixel to then add it to the noise map of the original image
    ; readout noise will not change wiht the insertion of aritificial stars
    ;------------------------------------------------------------------------------------------
    
    ;YO NO TENGO EXPOSURE MAPS, ASÍ QUE LO ADAPTO COMO ME HA DICHO CHATGPT, UTILIZANDO EL EXP TIME LEYÉNDOLO DEL HEADER DE LA IMAGEN 
    ;La gain en principio también se lee del header
    print, 'Now creating noise...' 
    ;art_noise = sources
    ;art_noise[*,*] = 0.0
    ;for k = 0, (n1-1) do begin
    ;   for j = 0, (n2-1) do begin
    ;      n_exp = exp[k,j]
          ; use ABS of image because PSF can have minor negativities
    ;      n_photons= abs(sources[k,j]) * n_exp/gain
    ;      if (n_photons gt 0) then $
    ;        sources[k,j] = gain/n_exp * randomu(seed,POISSON=n_photons)
    ;        art_noise[k,j] = gain/n_exp * sqrt(n_photons)
    ;  endfor
    ;endfor
    ;artim = im + sources
    ;art_noise = sqrt(noise^2 + art_noise^2)
    ;writefits, tmpdir + 'noisemap.fits', art_noise
    ;writefits, tmpdir + 'artim.fits', artim
    
    ;COPIO LO DE CHATGPT -> en principio está todo bien a nivel unidades porque las imágenes están en e-/s
    ; inicializa el mapa de ruido de las estrellas artificiales
    art_noise = sources * 0.0

    ; convierte la imagen de estrellas artificiales a electrones totales
    sources_e = sources * exp
    noise_e   = noise   * exp  ; convierte el mapa de errores a electrones

    ;n1 = n_elements(sources_e[*,0])
    ;n2 = n_elements(sources_e[0,*])

    for k = 0, n1-1 do begin
      for j = 0, n2-1 do begin
        ; usa abs() por si la PSF tiene negativos
        n_photons = abs(sources_e[k,j])
        if n_photons gt 0 then begin
          ; genera Poisson
          sources_e[k,j] = randomu(seed, POISSON=n_photons)
          ; guarda la contribución del ruido
          art_noise[k,j] = sqrt(n_photons)
        endif
      endfor
    endfor

    ; convierte de nuevo a e⁻/s para añadir a la imagen original
    artim = im + sources_e / exp

    ; combina el ruido original con el de las estrellas artificiales
    art_noise = sqrt( noise^2 + (art_noise / exp)^2 )

    ; guarda resultados
    writefits, tmpdir + 'artim.fits', artim
    writefits, tmpdir + 'noisemap.fits', art_noise

   ; now run StarFinder to recover sources
   ; Use EXACT parameters as in astrophot.pro
   ; implement any necessary additoinal source selection criteria
   ; -------------------------------------------------------------
   compbg = 1 ;keep
   ;back_box = maskrad ; En astrophot tengo 11 puesto aquí
   back_box = 11
   min_correlation = 0.85 ; El original era 0.7
   correl_mag = 4  ;keep  
   deblend = 1 ;keep
   deblost = 0 ;keep
   niter = 2 ;keep
   rel_thresh = 1 ;keep
   guide_x = "" ;keep
   guide_y = "" ;keep

   threshold = [3.,3.] ;keep
   starfinder, artim, psf, X_BAD=x_bad, Y_BAD = y_bad, $ ;keep
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = (art_noise + noise), min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLOST = deblost, DEBLEND = deblend, N_ITER = niter, /SILENT, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC
    writefits, tmpdir + 'stars'+strn(thismag)+'_'+strtrim(string(grid),2) + '.fits', stars
    writefits, tmpdir + 'bg.fits', background

   ; compare list of detected stars with list of artificial stars
   ; ------------------------------------------------------------
   x1 = xin
   y1 = yin
   x2 = x
   y2 = y
   compare_lists, x1, y1, x2, y2, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2
   nc = n_elements(subc1)
   ; save positions and fluxes of recovered artificial stars
   if (subc1[0] gt -1) then begin
    xout = xin[subc1]
    yout = yin[subc1]
    fout = f[subc1]

    ; Only accept recovered stars if they have the correct magnitude
    ; -------------------------------------------------------------------
    ;dmag = 2.5/alog(10) * abs(fin[subc1] - f[subc2])/fin[subc1]
    dmag = -2.5*alog(10)*(((abs(fin[subc1] - f[subc2])/fin[subc1])*PIXAR_SR)/f0) ; Creo que en realidad da igual aquí poner todas las 
    ;constantes porque es una diferencia, pero bueno por si acaso Im not really sure tbh
    good = where(dmag lt dmagmax,gcount)
    if (gcount gt -1) then begin
     xout = xout[good]
     yout = yout[good]
     fout = fout[good]
     outmap[round(xout),round(yout)] = 1      ; set pixel value of outmap to 1 if a star is recovered
    endif else begin
     print, 'No stars detected:Flux condition. STOP...'
     STOP
    endelse
   endif else begin
    print, 'No stars detected: Position condition. STOP...'
    STOP
   endelse
;   writefits, '../tmp/outmap'+strn(thismag)+'_'+strtrim(string(grid),2) + '_art.fits', outmap
   
   recovered_fluxes = [recovered_fluxes,fout] ; save fluxes of recovered artificial stars
   
  endfor  ; loop over y-grid positions
 endfor   ; loop over x-grid positions

 ; ONLY NEEDED FOR GNS
 ; expand inmap and outmap to size of original image
 ; -----------------------------------------------------------
 ;inmap_large = replicate(1,n1,n2)
 ;outmap_large = replicate(1,n1,n2)
 ;inmap_large[x_min:x_max,y_min:y_max] = inmap
 ;outmap_large[x_min:x_max,y_min:y_max] = outmap
 
 ;CLARO, YO CREO QUE SOLO NECESITO LO DE REPLICATE Y LISTO, O SEA, LO DE INMAP_LARGE Y OUTMAP_LARGE

 ; and save inmap and outmap
 ; -----------------------------------------------------------
 writefits, out_path + 'inmap' + strn(thismag) + '.fits', inmap, /COMPRESS
 writefits, out_path + 'outmap'+ strn(thismag) + '.fits', outmap, /COMPRESS

 print, 'Finished magnitude ' + strn(thismag) + '.'
; STOP
;endfor  ; loop over magnitude steps


END

