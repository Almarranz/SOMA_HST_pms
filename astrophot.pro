PRO ASTROPHOT, zone, filter, epoch

; zone = 'G028.20-00.05'
; filter = '160w'
; epoch = '1'
; path = '/Users/fedriani/Documents/postdoc_iaa/HST_project/HST_data/G028.20-00.05/gaia_alignment/Epoch1/starfinder/'
path = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'+ zone +'/gaia_alignment/Epoch'+ epoch +'/'
pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/pruebas/'
tmpdir = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'/epoch'+ epoch +'/tmp/'
results = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/'+ zone +'/f'+ filter +'/epoch'+ epoch +'/'



; nam = zone +'_EP'+ epoch + '_f'+ filter + '_drz_sci'
nam = 'hst_ep'+ epoch + '_f'+ filter + '_drz'


print, nam

  ; create tmp directory if necessary
;   if not(FILE_TEST(path + 'tmp')) then FILE_MKDIR, path + 'tmp'
  if not(FILE_TEST(results + 'tmp')) then FILE_MKDIR, results + 'tmp'
  

;   psf = readfits(path + 'psf_'+filter+'.fits') ;use Rainer's PSF worst case scenerio. psf = readfits(path + 'psf_'+filter+'_rainer.fits')
  psf = readfits(tmpdir + 'psf_'+filter+'.fits') ;use Rainer's PSF worst case scenerio. psf = readfits(path + 'psf_'+filter+'_rainer.fits')

  ; load image and noise map
;   im = readfits(path + nam + '.fits',header,EXT=0) ;careful here! check the extension!
  im = readfits(path + nam + '.fits',header,EXT=1) ;careful here! check the extension!
  good = where(FINITE(im),complement=isnan)
  im[isnan] = 0
  noise = sqrt(im)
;   writefits, path + 'im.fits', im, /COMPRESS
;   writefits, path + 'noise.fits', noise, /COMPRESS
 
  writefits, tmpdir + 'im.fits', im
  writefits, tmpdir + 'noise.fits', noise
 
  sz = size(im)
  n1 = sz[1]
  n2 = sz[2]

; Settings for StarFinder you are most likely to 
; want to play with
; these parameters apply to the final StarFinder run
; not to PSF extraction
; Note that the thesholds are set very low, to 0.5 sigma (lower is
; probably even possible)
; This means that our simple noise map overestiamtes the actual uncertainties
; ------------------------------------------------

  sf_thresh = [.5,.5]
;   sf_thresh = [1.,1.]
  back_box = 20
  deblend = 0
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


; StarFinder run
;####################

  starfinder, im, psf, X_BAD=xbad, Y_BAD = ybad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        sf_thresh, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLEND = deblend, DEBLOST = deblost, $
        N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC, $
        XBAD = xbad, YBAD = ybad, POSBG = posbg

  writefits, tmpdir + STRMID(nam,0,100) + '_stars'+filter+'.fits', stars, header
  writefits, tmpdir + STRMID(nam,0,100) + '_bg'+filter+'.fits', background, header
  writefits, tmpdir + STRMID(nam,0,100) + '_resid'+filter+'.fits', im-stars-background, header
  subtracted = im-stars
  writefits, tmpdir + STRMID(nam,0,100) + '_subtracted'+filter+'.fits', subtracted, header
  
;   writefits, path + STRMID(nam,3,100) + '_stars'+filter+'.fits', stars, header, /COMPRESS
;   writefits, path + STRMID(nam,3,100) + '_bg'+filter+'.fits', background, header, /COMPRESS
;   writefits, path + STRMID(nam,3,100) + '_resid'+filter+'.fits', im-stars-background, header, /COMPRESS
;   subtracted = im-stars
;   writefits, path + STRMID(nam,3,100) + '_subtracted'+filter+'.fits', subtracted, header, /COMPRESS
;   
  ; save list
  ; select stars in region with more than covfrac coverage
  nstars = n_elements(f)
;   openw, outp, path + STRMID(nam,3,100) + '_stars'+filter+'.txt', /get_lun
  openw, outp, results + STRMID(nam,0,100) + '_stars'+filter+'.txt', /get_lun
  printf, outp, 'x y f sx sy sf'
  for s = 0, nstars-1 do begin
   xi = round(x[s]) & yi = round(y[s])
   printf, outp, format='(6f13.3)', x[s], y[s], f[s], sx[s], sy[s], sf[s]
  endfor
  free_lun, outp

; make map of point sources with Gaussian PSFs (helpful for extremely
; crowded fields)
; -------------------------------------------------------------
;   readcol, path + STRMID(nam,3,100) + '_stars'+filter+'.txt', x, y, f
  readcol, results + STRMID(nam,0,100) + '_stars'+filter+'.txt', x, y, f,SKIPLINE =1, Format ='A,A,A,A,A,A'
  x=float(x)
  y=float(y)
  f=float(f)
  dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 1.0, Sigma_y: 1.0, Angle: 0.0})
  im = image_model(x,y,f,n1,n2,'gaussian', dat)
;   writefits, path + STRMID(nam,3,100) + '_map'+filter+'.fits', im, /COMPRESS
  writefits, tmpdir + STRMID(nam,0,100) + '_map'+filter+'.fits', im 


  print, 'Finished ' + nam + '.'

END
