PRO COMPLETENESS_hst ;no necesito el zp, necesito f0 :)

t0 = SYSTIME(/SECONDS)
;chip = namw
;f0 = vega fluxes
thismag = 17;the mag. you want to calculate the completenes on
; file_path = 'names.txt'
; 
; Read the values from the file
; The FORMAT='A' option reads the values as strings
; READCOL, file_path, band, field, FORMAT='A,A' ;Original: survey, band, field FORMAT = 'A,A,A' ; Ahora leemos f0, que cogemos el mismo
;para todos los chips de la banda y lo necesito para calcular el flujo. NO, EL F0 LO METO A MANO 
; print, band

; filters = ['160w', '110w', '128n', '164n']

; filters = ['160w', '110w']
; filters = ['128n', '164n']
; filters = ['128n']
filters = ['160w']

; zones = [ 'AFGL5180','G028.20-00.05', 'G032.03+00.05', 'G35.2-0.74N', 'G339.88-01.26','IRAS07299-1651', 'IRAS16562-3959']
zones = [ 'G028.20-00.05']


epoch = '1'
filter = '160w'
zone ='G028.20-00.05'



basedir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/' + zone + '/epoch' + epoch + '/'
in_path = basedir                                ;basedir + 'cubeims/chip' + strn(chip) + '/'
out_path =  '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/' + zone + '/F' + filter + '/epoch' + epoch +  '/completeness/'

tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/' + zone + '/F' + filter + '/epoch' + epoch + '/tmp/'
comdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/' + zone + '/F' + filter + '/epoch' + epoch + '/completeness/'

readcol, tmpdir + 'ZP_' + zone + '_' + filter + '_' + epoch + '.txt', zp, Format ='A'
ZP = float(zp)
print, ZP, zone + '_' + filter + '_' + epoch 

pattern = basedir + 'hst_*' + filter + '*'
paths = FILE_SEARCH(pattern, /TEST_DIRECTORY)
IF (N_ELEMENTS(paths) NE 1) OR (paths[0] EQ '') THEN $
  MESSAGE, 'Expected one image directory matching: ' + pattern

path = paths[0] + '/'
nam = FILE_BASENAME(paths[0]) + '_drz'
filename = path + nam + '.fits'

print,'***** Working on ******************************************************'
print, filename
print,'***********************************************************************'

; Important parameters
; ----------------------



; Create a grid of artificial stars that is sufficiently wide to not increase crowding
; to explore completeness as a funciton of position better, this grid is shifted
; around the image in separate steps
; ==========================================================================================
step = 26L                ; grid step size in pixels CHANGE -> 10xFWHM, entonces sería 20, porque yo no tengo el rebining factor :). Lo hemos subido un poquito para que sea ~0.8"
gridstep = step/4         ; distance by which to shift the grid in 4 separate runs -> KEEP; shifts de ~2"
randstep = 0.5            ; random jitter of artificial sources around grid positions -> KEEP, doesn't matter much
magstep = 1             ; completeness will be probed in steps of 0.5 mag -> KEEP in principle, unless it is too slow; quizás es incluso mucho... pero bueno siempre s epuede interpolar a unas malas 
mag_min = 9              ; brightest magnitude to examine -> SELECT FROM PRELIMINAR (REDUNDANT) COMPLETENESS!!
nmag = 6              ; number of magnitudes to examine -> UNTIL REACHING WHERE THE LF FALLS!!! 
gain = 2.0               ; approximate gain for all detectors (precise value should not matter much)

; parameters to define redeteciton of an artificial star --> ESTO POR AHORA LO MANTENGO, A VER Cómo workea 
; ------------------------------------------------------
dmagmax = magstep/2. ; maximum accepted magnitude difference between recovered source and input source
dmax = 1.5           ; maximum position offset between recovered source and input source - USE VALUE FROM COMPUTE_UNCERTAINTIES.PRO

; noise map
; ------------------------------------------------------
im = readfits(path + nam + '.fits',header,EXT=1) ;careful here! check the extension!
print,sxpar(header,'BUNIT')
im_ext0 = readfits(path + nam + '.fits',header0,EXT=0)
exp = sxpar(header0,'EXPTIME')

finite_mask = FINITE(im)

mask_nan = fltarr(size(im, /DIMENSIONS))
mask_nan[*] = 1.0
bad = WHERE(finite_mask EQ 0, nbad)
IF nbad GT 0 THEN mask_nan[bad] = 0.0
; writefits, tmpdir + 'mask_nan.fits', mask_nan

print, sxpar(header0,'EXPTIME')

sz = size(im)
n1 = sz[1] ;Antes era n1_0 y n2_0, pero como yo no tengo que recortar la imagen pues directamente n1,n2
n2 = sz[2]

print,n1,n2

good = where(FINITE(im),complement=isnan)
im[isnan] = 0
; noise = sqrt(im)
noise = sqrt(abs(im*exp))/exp
;
;LECTURA PSF 
psf = readfits(tmpdir + 'psf_gui_'+filter+'.fits')



; loop over magnitudes
; ----------------------------------------------------------------------------------

; for ms = 1, nmag do begin

   
;  thismag = mag_min + ms * magstep      ; magnitude probed in this loop
  flux = 10^(0.4*(ZP-thismag))
  
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
;   for xg = 0, 1 do begin
;     for yg = 0, 1 do begin
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
     sources = sources*mask_nan
     inmap = inmap*mask_nan
     writefits, comdir + 'sources_mag'+strn(thismag)+'.fits', sources    
     writefits, comdir + 'inmap'+strn(thismag)+'.fits', inmap    
     
    ; add POISSON noise to image of artificial stars
    ; save POISSON noise for each pixel to then add it to the noise map of the original image
    ; readout noise will not change wiht the insertion of aritificial stars
    ;------------------------------------------------------------------------------------------
    
     print, 'Now creating noise...' 
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
    writefits, comdir + 'artim_mag'+strn(thismag)+'.fits', artim
    writefits, comdir + 'noisemap_mag'+strn(thismag)+'.fits', art_noise
;    stop
   ; now run StarFinder to recover sources
   ; Use EXACT parameters as in astrophot.pro
   ; implement any necessary additoinal source selection criteria
   ; -------------------------------------------------------------
   sf_thresh = [.5,.5]
 ;   sf_thresh = [1.,1.]
   back_box = 20
   deblend = 0
   deblost = 0
   compbg = 1
   posbg = 1
   
   compbg = 1 ;keep
   ;back_box = maskrad ; En astrophot tengo 11 puesto aquí
   back_box = 20
   min_correlation = 0.9 ; El original era 0.7
   correl_mag = 4  ;keep  
   deblend = 1 ;keep
   deblost = 0 ;keep
   niter = 2 ;keep
   rel_thresh = 1 ;keep
   guide_x = "" ;keep
   guide_y = "" ;keep
   
  
   
    
   

;    threshold = [3.,3.] ;keep
;    starfinder, artim, psf, X_BAD=x_bad, Y_BAD = y_bad, $ ;keep
;         BACKGROUND = background, BACK_BOX = back_box, $
;         threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
;         NOISE_STD = (art_noise + noise), min_correlation, $
;         CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
; 	ESTIMATE_BG = compbg, DEBLOST = deblost, DEBLEND = deblend, N_ITER = niter, /SILENT, $
; 	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
; 	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
;       	x, y, f, sx, sy, sf, c, STARS = stars, $
;         LOGFILE = logfilename, /CUBIC


    starfinder, artim, psf, X_BAD=xbad, Y_BAD = ybad, $
          BACKGROUND = background, BACK_BOX = back_box, $
          sf_thresh, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
          NOISE_STD = art_noise, min_correlation, $
          CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
  	ESTIMATE_BG = compbg, DEBLEND = deblend, DEBLOST = deblost, $
          N_ITER = niter, SILENT=0, $
  	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
  	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        	x, y, f, sx, sy, sf, c, STARS = stars, $
          LOGFILE = logfilename, /CUBIC, $
          XBAD = xbad, YBAD = ybad, POSBG = posbg

    writefits, comdir + 'stars'+strn(thismag)+'_'+strtrim(string(grid),2) + '.fits', stars
    writefits, comdir + 'bg_mag'+strn(thismag)+'.fits', background


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
    dmag = 2.5/alog(10) * abs(fin[subc1] - f[subc2])/fin[subc1]
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
 writefits, comdir + 'inmap' + strn(thismag) + '.fits', inmap
 writefits, comdir + 'outmap'+ strn(thismag) + '.fits', outmap
;  writefits, out_path + 'inmap' + strn(thismag) + '.fits', inmap, /COMPRESS
;  writefits, out_path + 'outmap'+ strn(thismag) + '.fits', outmap, /COMPRESS

 print,'Injected:',N_ELEMENTS(xin)
 print,'Recovered:',N_ELEMENTS(xout)
 print,'Fraction:',FLOAT(N_ELEMENTS(xout))/FLOAT(N_ELEMENTS(xin))
 
n_in  = total(inmap)
n_out = total(outmap)

print,'Injected total = ',n_in
print,'Recovered total = ',n_out
print,'Fraction = ',float(n_out)/float(n_in)

 print, 'Finished magnitude ' + strn(thismag) + '.'

elapsed = SYSTIME(/SECONDS) - t0
mins = FLOOR(elapsed/60.)
secs = elapsed - 60.*mins
print,'========================================='
print, 'Runtime mag',thismag,':', mins, ' min ', secs, ' s' 
print,'========================================='
; STOP
; endfor  ; loop over magnitude steps


END

