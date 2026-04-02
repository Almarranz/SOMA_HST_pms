import sys
sys.path.append("/Users/amartinez/Desktop/pythons_imports/")
import numpy as np
from astropy import units as u
from astropy.time import Time
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import sys
from matplotlib import rcParams
from astroquery.gaia import Gaia
import IPython
import os
# import cluster_finder
from filters import filter_gaia_data
import Polywarp as pw
import skimage as ski
from astropy.table import Table
from compare_lists import compare_lists
from astropy.stats import sigma_clip
from alignator_gaia import alig_gaia
from astropy.table import unique
import gns_cluster_finder
from filters import filter_hst_data
from astropy.modeling.models import Polynomial2D
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling import models, fitting
from astropy.coordinates import search_around_sky
from collections import Counter
from sys import exit as stop
from astropy.coordinates import Longitude
from astropy.io import fits
from astropy.wcs import WCS
from pyplots import plot_two_hists_sigma
from pyplots import plot_two_pm_hists
from alignator_looping import alg_loop
from ds9_region import region_vectors
from ds9_region import region
from matplotlib.colors import LogNorm
import cluster_finder
from hst_irZP import get_vegazp
import astroalign as aa
from pyplots import plot_two_hists_sigma
from pyplots import sig_cl
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams

# plt.rcParams["mathtext.fontset"] = 'dejavuserif'
# rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
# # Enable automatic plotting mode
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# # IPython.get_ipython().run_line_magic('matplotlib', 'inline')


folder = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'
# results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
# results = ''
# zone = 'G032.03+00.05'
zone = 'G028.20-00.05'
band = '160w'
# band = '128n'
# epoch = 2

pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'


# =============================================================================
# Gaia parametres
# =============================================================================
radius = 200*u.arcsec
max_sep = 50*u.mas
mag_gaia = [10,18]
e_pm_gaia = 0.5
e_pos_gaia = 0.5
# =============================================================================
# HST observation
# =============================================================================

pixSca = 0.12825 #arcsec/pixel
# pixSca = 0.00001 #arcsec/pixel#

e_pos_cat  = 0.05# in arcsec. The position errors from starfinder lists are largely overstatimated!!!

red_techn = 'Gaia'
# red_techn = 'Original'
# =============================================================================
# Aligment paremeters
# =============================================================================
# pre_transf = 'polynomial'
pre_transf = None

# transf = 'polynomial'
# transf = 'affine' # Transformation before loopings
transf = 'similarity'
order_trans = 1
align_loop = 'yes'
# align_loop = 'no'
align = 'Polywarp'
max_deg = 3# If this is <2 it does not enter the alignment loop. 
max_loop = 3
# gaia_clipping = 'one_one'# Clipp the Gaia outlayer one by one
gaia_clipping = 'all'# Clipp the Gaia outlayer all at once
centered_in = 1
destination = 1
sig_cl_H = 3
d_fa = 50*u.mas
align_by = 'Polywarp'#!!!
# align_by = '2DPoly'#!!!
f_mode = 'W' # f_mode only useful for 2Dpoly
# f_mode = 'WnC'
# f_mode = 'NW'
# f_mode = 'NWnC'
# =============================================================================
# GRID PARAMS
# =============================================================================
grid_s = None
# grid_s = 2# si1ce of the grid cell in arcsec 
grid_Hmin = 12
grid_Hmax = 18
isolation_radius = 0.7#arcsec isolation of the grid stars 
# =============================================================================
# Proper motions param
# ===========================================================================
max_dis_pm = 0.150#in arcsec
sig_H = 3# discrd pm for stars with delta H over sig_H
e_pm_cat = 20# im mas/yr
d_pm = 150*u.mas
# =============================================================================
# CLUSTERS
# =============================================================================
# look_for_cluster = 'no'
look_for_cluster = 'yes'

# =============================================================================
# Dictionaries
# =============================================================================
cat_dic = {}
obst_dic = {}
gaia_dic = {}
# center = SkyCoord(ra = 282.41000248, dec = -0.78212655, unit = 'degree', frame = 'icrs')


bad = []
lopping = 1
# for loop in range(1):
wloop_counter = 0
# ZP = 25 # INVENTED


for epoch in range(1,3):
    results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
    
    if red_techn == 'Gaia':
        cat = Table.read(results + f'{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii')
        ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_f{band}_drz_sci.fits')
        wcs = WCS(ima[0].header)
        scale_pix = np.sqrt(ima[0].header['CD1_1']**2 + ima[0].header['CD2_1']**2 )*3600
        
    if red_techn == 'Original':
        cat = Table.read(results + f'hst_ep{epoch}_f{band}_drz_stars{band}.txt', format = 'ascii')
        ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/hst_ep{epoch}_f{band}_drz.fits')
        wcs = WCS(ima[1].header)
        scale_pix = np.sqrt(ima[1].header['CD1_1']**2 + ima[1].header['CD2_1']**2 )*3600
    
    print(pixSca, scale_pix)
# =============================================================================
#         ZP calculation
# =============================================================================
    mjd = ima[0].header['EXPSTART']
    ZP = get_vegazp(mjd, band = f'f{band}')
    print(ZP)


    cat['H'] =  (-2.5*np.log10(cat['f']) + ZP).round(4)
    cat['dH'] = ((2.5 / np.log(10)) * (cat['sf'] / cat['f'])).round(4)
    
    cat.write(results + f'calib_{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii', overwrite = True )
    

    
    cat_rd = wcs.pixel_to_world(cat['x'], cat['y'])
    
    if epoch == 1:
        center = SkyCoord(ra = np.mean(cat_rd.ra.value), dec = np.mean(cat_rd.dec.value),
                          unit = 'degree', frame = 'icrs')
    
    cat['ra'] = cat_rd.ra
    cat['dec'] = cat_rd.dec
    
    cat['x'] = cat['x']+1 # This +1 is just to make the region file.
    cat['y'] = cat['y']+1
    
    region(cat, 'x', 'y',
           name = f'{zone}_Ep{epoch}_stars_xy_0.5s',
           save_in = pruebas,
           wcs = 'physical',
           color = 'cyan',
           marker = 'cross')
    
    cat['x'] = cat['x']-1 
    cat['y'] = cat['y']-1
    
    cat['x'] = cat['x']*pixSca # These are arcsec
    cat['y'] = cat['y']*pixSca
    
    cat['sx'] = cat['sx']*pixSca # These are arcsec
    cat['sy'] = cat['sy']*pixSca
    
    cat['sxy'] = np.sqrt(cat['sx']**2 + cat['sy']**2)
    cat = cat[cat['sxy'] > 0]
    
    fig, (ax, ax2) = plt.subplots(1,2, figsize = (12,6))
    ax2.set_title(f'ZP = {ZP: .3f}')
    ax.set_title(f'Epoch{epoch}')
    h = ax.hexbin(cat['H'], cat['sxy'], norm = LogNorm())
    cbar = plt.colorbar(h, ax=ax)
    ax.set_ylabel('$\sqrt{\sigma{x}^2 + \sigma{y}^2}$ [arcsec]', fontsize = 12)
    ax.set_xlabel('[H]', fontsize = 12)
    ax.axhline(e_pos_cat, ls = 'dashed', color = 'red')
    
    h2 = ax2.hexbin(cat['H'], cat['dH'], norm = LogNorm(), cmap = 'inferno')
    cbar = plt.colorbar(h2, ax=ax2)
    ax2.set_ylabel('[dH]', fontsize = 12)
    ax2.set_xlabel('[H]', fontsize = 12)
    # ax2.axhline(e_pos_cat, ls = 'dashed', color = 'red')
    fig.tight_layout()
    mjd = ima[0].header['EXPSTART']   # e.g., 57610.4918468
    obst = Time(mjd, format='mjd', scale='utc')
    obst_dic[f'tyr{epoch}'] = obst.decimalyear
    

    cat = filter_hst_data(cat, max_e_pos = e_pos_cat)
   
    obst_dic[f't{epoch}'] = obst
    cat_dic[f'cat{epoch}'] = cat

cat1 = cat_dic['cat1'] 
cat2 = cat_dic['cat2'] 
# %%
fig, ax = plt.subplots(1,1)
ax.scatter(cat1['ra'], cat1['dec'])
ax.scatter(cat2['ra'], cat2['dec'])
# ax.scatter(cat1['x'], cat1['y'])
# ax.scatter(cat2['x'], cat2['y'])
ax.axis('scaled')

centered_in = 1
if centered_in == 1:
    center = SkyCoord(ra = np.mean(cat1['ra']), dec = np.mean(cat1['dec']), unit = 'degree')
elif centered_in == 2:
    center = SkyCoord(ra = np.mean(cat2['ra']), dec = np.mean(cat2['dec']), unit = 'degree')
    
# ax.scatter(center.ra, center.dec, s = 400)

cat1_rd = SkyCoord(ra = cat1['ra'], dec = cat1['dec'], unit = 'degree')
cat2_rd = SkyCoord(ra = cat2['ra'], dec = cat2['dec'], unit = 'degree')



xg_1, yg_1 = center.spherical_offsets_to(cat1_rd)
xg_2, yg_2 = center.spherical_offsets_to(cat2_rd)

tag = center.skyoffset_frame()

cat1_t = cat1_rd.transform_to(tag)
cat2_t = cat2_rd.transform_to(tag)


cat1['xp'] = cat1_t.lon.to(u.arcsec)
cat1['yp'] = cat1_t.lat.to(u.arcsec)
cat2['xp'] = cat2_t.lon.to(u.arcsec)
cat2['yp'] = cat2_t.lat.to(u.arcsec)

# %%

idx1, idx2, sep2d, _ = search_around_sky(cat1_rd, cat2_rd, max_sep)

count1 = Counter(idx1)
count2 = Counter(idx2)

# Step 3: Create mask for one-to-one matches only
mask_unique = np.array([
    count1[i1] == 1 and count2[i2] == 1
    for i1, i2 in zip(idx1, idx2)
])

# Step 4: Apply the mask
idx1_clean = idx1[mask_unique]
idx2_clean = idx2[mask_unique]

cat1_match= cat1[idx1_clean]
cat2_match = cat2[idx2_clean]


diff_H = cat2_match['H'] - cat1_match['H']

mask_H, l_lim,h_lim = sigma_clip(diff_H, sigma=sig_cl_H, masked = True, return_bounds= True, maxiters= 50)

cat2_match = cat2_match[np.logical_not(mask_H.mask)]
cat1_match = cat1_match[np.logical_not(mask_H.mask)]


# %%
fig,(ax,ax2) = plt.subplots(1,2)
# ax.set_title(f'Matching starts {len(cat2_match)}')
ax.hist(diff_H[np.logical_not(mask_H.mask)], bins = 'auto',histtype = 'step', label = f'\u2605 = {len(cat2_match)}')
ax.hist(diff_H, bins = 'auto', color = 'grey', alpha = 0.3)
ax.axvline(np.mean(diff_H), color = 'k', ls = 'dashed', label = '$\overline{\Delta H}$= %.2f$\pm$%.2f'%(np.mean(diff_H),np.std(diff_H)), lw = 0.5)
ax.axvline(l_lim, ls = 'dashed', color ='r', label ='%s$\sigma$'%(sig_cl_H))
ax.set_xlim(-1,1)
ax2.hist2d(cat2_match['H'],diff_H[np.logical_not(mask_H.mask)], bins = 20, norm = LogNorm())
ax2.set_ylim(-1,1)
ax.axvline(h_lim, ls = 'dashed', color ='r')
ax.set_xlabel('$\Delta H$')
ax.legend() 

c1_m = np.array([cat1_match['xp'],cat1_match['yp']]).T
c2_m = np.array([cat2_match['xp'],cat2_match['yp']]).T

sig_cl_H_aligment = sig_cl_H
if destination == 1:
    # Time lapse to move Gaia Stars.
    
    dt = (obst_dic['tyr1'] - obst_dic['tyr2'])*u.yr
    
    p,(_,_)= aa.find_transform(c2_m,c1_m,max_control_points=100)
    
    # p = ski.transform.estimate_transform('similarity',
    #                                 g2_m, 
    #                                 g1_m)
    # 
    print("Translation: (x, y) = (%.2f, %.2f)"%(p.translation[0],p.translation[1]))
    print("Rotation: %.2f deg"%(p.rotation * 180.0/np.pi)) 
    print("Rotation: %.0f arcmin"%(p.rotation * 180.0/np.pi*60)) 
    print("Rotation: %.0f arcsec"%(p.rotation * 180.0/np.pi*3600)) 
    
    
    
    
    loop = 0
    comom_ls = []
    dic_xy = {}
    dic_Kx ={}
    dic_xy_final = {}
    
        
    cat2_xy = np.array((cat2['xp'],cat2['yp'])).T
    cat2_xyt = p(cat2_xy)
    
    s_ls = compare_lists(cat2_xyt, np.array([cat1['xp'],cat1['yp']]).T, d_fa.to(u.arcsec).value)
    print(f'Common stars after astroaling similaryty:{len(s_ls)}')
    
    cat2['xp'] = cat2_xyt[:,0]
    cat2['yp'] = cat2_xyt[:,1]
    
    
    # gns2 = alg_rel(gns2, gns1,'xp', 'yp', align_by,use_grid,max_deg = max_deg, d_m = d_m,f_mode = f_mode,grid_s = grid_s )
    # def alg_loop(gns_A, gns_B,col1, col2, align_by, max_deg, d_m,                        max_loop,  use_grid,grid_s= None, f_mode = None  ) :
    cat2 = alg_loop(cat2, cat1, 'xp', 'yp', align_by, max_deg, d_fa.to(u.arcsec).value, max_loop,sig_cl_H = sig_cl_H_aligment, 
                    grid_s = grid_s, grid_Hmin = grid_Hmin, grid_Hmax = grid_Hmax ,isolation_radius = isolation_radius,  f_mode = f_mode, mag_lim_alig= None)


l1_xy = np.array([cat1['xp'], cat1['yp']]).T
l2_xy = np.array([cat2['xp'], cat2['yp']]).T
l_12 = compare_lists(l1_xy,l2_xy,d_pm.to(u.arcsec).value)


print(30*'*'+'\nComon stars to be use for pm calculation :%s\n'%(len(l_12))+30*'*')
cat1_mi = cat1[l_12['ind_1']]
cat2_mi = cat2[l_12['ind_2']]



dx = (cat1_mi['xp'].value - cat2_mi['xp'].value)*1000
dy = (cat1_mi['yp'].value - cat2_mi['yp'].value)*1000

cat1_mi['sx'].unit = u.mas
cat1_mi['sy'].unit = u.mas
cat2_mi['sx'].unit = u.mas
cat2_mi['sy'].unit = u.mas

dpm_x = np.sqrt((cat2_mi['sx'].to(u.mas))**2 + (cat1_mi['sx'].to(u.mas))**2)/dt.to(u.year)
dpm_y = np.sqrt((cat2_mi['sy'].to(u.mas))**2 + (cat1_mi['sy'].to(u.mas))**2)/dt.to(u.year)


pm_x = (dx*u.mas)/dt.to(u.year)
pm_y = (dy*u.mas)/dt.to(u.year)

cat1_mi['pm_x']  = pm_x
cat1_mi['pm_y']  = pm_y

cat2_mi['pm_x']  = pm_x
cat2_mi['pm_y']  = pm_y

cat1_mi['dpm_x']  = dpm_x
cat1_mi['dpm_y']  = dpm_y

cat2_mi['dpm_x']  = dpm_x
cat2_mi['dpm_y']  = dpm_y

# %%


fig, (ax,ax2) = plt.subplots(1,2, figsize = (7,3.5))
ax.hist(cat1_mi['pm_x'], histtype = 'step', bins = 'auto', label = f'$\mu_{{x}} = {np.mean(cat1_mi["pm_x"]): .2f}$\n$\sigma_{{x}} = {np.std(cat1_mi["pm_x"]): .2f}$')
ax.set_xlabel('$\mu_{x}$ [mas/yr]')
ax.set_ylabel('#')
ax.legend()
ax2.hist(cat1_mi['pm_y'], histtype = 'step', bins = 'auto', label = f'$\mu_{{y}} = {np.mean(cat1_mi["pm_y"]): .2f}$\n$\sigma_{{x}} = {np.std(cat1_mi["pm_y"]): .2f}$')
ax2.set_xlabel('$\mu_{y}$ [mas/yr]')
ax2.legend()
fig.tight_layout()
# %%
# =============================================================================
# GAIA
# =============================================================================


try:
    
    gaia = Table.read(pruebas  +f'gaia_{zone}_{radius.value: .0f}.ecsv', format = 'ascii.ecsv')
    # gaia = Table.read(pruebas1  + 'gaia_f1%s_f2%s_r400.ecsv'%(field_one,field_two))
    
    # gaia = Table.read('/Users/amartinez/Desktop/PhD/HAWK/GNS_1relative_SUPER/pruebas/gaia_f1D19_f216_r269.ecsv')
    print('Gaia from table')
except:
    print('Gaia from web')
    # center = SkyCoord(l = np.mean(gns1['l']), b = np.mean(gns1['b']), unit = 'degree', frame = 'galactic').icrs

    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
    Gaia.ROW_LIMIT = -1  # it not especifty, Default rows are limited to 50. 
    # j = Gaia.cone_search_async(center, radius = abs(radius))
    j = Gaia.cone_search_async(center, radius = radius)
    gaia = j.get_results()
    os.makedirs(pruebas, exist_ok=True)
    gaia.write(pruebas  + f'gaia_{zone}_{radius.value: .0f}.ecsv', format = 'ascii.ecsv', overwrite = True)
    gaia = Table.read(pruebas  +f'gaia_{zone}_{radius.value: .0f}.ecsv', format = 'ascii.ecsv')

gaia = filter_gaia_data(
    gaia_table=gaia,
    astrometric_params_solved=31,
    duplicated_source= False,
    parallax_over_error_min=-10,
    astrometric_excess_noise_sig_max=2,
    phot_g_mean_mag_min= mag_gaia[1],
    phot_g_mean_mag_max= mag_gaia[0],
    pm_min=None,
    pmra_error_max=e_pm_gaia,
    pmdec_error_max=e_pm_gaia,
    ra_error_max = e_pos_gaia,
    dec_error_max = e_pos_gaia,
    # min_angular_separation_arcsec = 0.1*u.arcsec
    )


gaia['id'] = np.arange(len(gaia))

gaia.sort('phot_g_mean_mag')
gfig, (gax, gax2) = plt.subplots(1,2)
gax.scatter(gaia['phot_g_mean_mag'],gaia['pmra_error'],color = 'k')
gax.scatter(gaia['phot_g_mean_mag'],gaia['pmdec_error'], color = 'grey')
# gax2.scatter(gaia['phot_g_mean_mag'],gaia['ra_error'],color = 'k')
# gax2.scatter(gaia['phot_g_mean_mag'],gaia['dec_error'],color = 'grey')
gax.set_ylabel('$\sigma \mu$')
gax2.set_ylabel('$\sigma \mu (gns)$')
gax.set_xlabel('[G]')
gax2.set_xlabel('[H]')
gfig.tight_layout()

# %%
fig, ax = plt.subplots(1,1)
ax.scatter(cat1_mi['ra'], cat1_mi['dec'], label = 'HST',s =1 )
ax.scatter(gaia['ra'], gaia['dec'],zorder = 2, label = 'Gaia', s =1 )
ax.legend()





gaia_rd = SkyCoord(ra = gaia['ra'], dec = gaia['dec'],
                   pm_ra_cosdec = gaia['pmra'],
                   pm_dec = gaia['pmdec'], 
                   frame = 'icrs', obstime = 'J2016')

gaia_rdt = gaia_rd.apply_space_motion(new_obstime = Time(obst_dic['tyr1'], format='decimalyear'))


gaia['ra1'] = gaia_rdt.ra
gaia['dec1'] = gaia_rdt.dec

xp_g, yp_g = center.spherical_offsets_to(gaia_rdt.frame)

gaia['x'] = xp_g.to(u.arcsec) 
gaia['y'] = yp_g.to(u.arcsec) 

# This proyect Gaia proper motions on the same tangentical plas of that of the proyected x,y coordenates
offset_frame = center.skyoffset_frame()

c_proj = gaia_rdt.transform_to(offset_frame)

# Step 3: Extract Gaia PM in tangent plane (same as your XY frame, in mas/yr)
pm_x_gaia = c_proj.pm_lon_coslat  # mas/yr
pm_y_gaia = c_proj.pm_lat        # mas/yr

gaia['pm_x'] = pm_x_gaia
gaia['pm_y'] = pm_y_gaia

gaia_c = SkyCoord(ra = gaia['ra1'], dec = gaia['dec1'], frame = 'icrs')
cat_c = SkyCoord(ra = cat1_mi['ra'], dec = cat1_mi['dec'], frame = 'icrs')

idx, d2d, _ = gaia_c.match_to_catalog_sky(cat_c, nthneighbor=1)
match_mask = d2d < max_sep
gaia_m = gaia[match_mask]
cat_m = cat1_mi[idx[match_mask]]
print(40*'+')
unicos = unique(cat_m, keep = 'first')
print(len(cat_m),len(unicos))
print(40*'+')



fig, (ax, ax1, ax2) = plt.subplots(1,3, figsize = (15,5))

# ax.scatter(gns_w[::100], gns_cat['b'][::100], alpha =0.1, color = 'k')
# ax.scatter(ga_w, gaia['b'], label = 'Gaia', s= 10)

ax.set_title(f'Matches = {len(cat_m)}\nMin dist = {max_sep} ')
ax.scatter(cat1_mi['ra'], cat1_mi['dec'])
ax.scatter(gaia_m['ra'], gaia_m['dec'],s =10, label = 'Gaia Match')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.legend()


diff_pmx_all = gaia_m['pmra'] - cat_m['pm_x']
diff_pmy_all = gaia_m['pmdec'] - cat_m['pm_y']
# diff_pmx_all = gaia_m['pm_x'] - cat_m['pm_x']
# diff_pmy_all = gaia_m['pm_y'] - cat_m['pm_y']


mask, l = sig_cl(diff_pmx_all, diff_pmy_all, 3)

diff_pmx = diff_pmx_all[mask]
diff_pmy = diff_pmy_all[mask]

# %%
fig, (ax, ax2) = plt.subplots(1,2, figsize = (10,5))
ax.hist(diff_pmx_all, bins = 'auto', color = 'gray',alpha = 0.3)
ax2.hist(diff_pmy_all, bins = 'auto', color = 'gray',alpha = 0.3)
ax.hist(diff_pmx, histtype = 'step', label = f'$\Delta\overline{{\mu}}_{{x}} = {np.mean(diff_pmx): .2f}$\n$\sigma_{{x}} = {np.std(diff_pmx): .2f}$')
ax2.hist(diff_pmy, histtype = 'step', label = f'$\Delta\overline{{\mu}}_{{y}} = {np.mean(diff_pmy): .2f}$\n$\sigma_{{x}} = {np.std(diff_pmy): .2f}$')
ax.set_xlabel('$\Delta \mu_{xp}$ [mas/yr]')
ax2.set_xlabel('$\Delta \mu_{yp}$ [mas/yr]')
ax.legend()
ax2.legend()

# %%

if look_for_cluster == 'yes':
    
   
    # modes = ['pm_xy_color']
    modes = ['pm_xy']
    knn = 30
    gen_sim = 'kernnel'
    sim_lim ='minimun'
    # sim_lim ='mean'

    clus_dic = cluster_finder.finder(cat1_mi, 'pm_x', 'pm_y',
                                 'x', 'y', 
                                 'ra', 'dec',
                                 modes[0],
                                 'H','H',
                                 knn,gen_sim,sim_lim, save_reg = None)
   
#     elif destination == 1:
#         clus_dic = gns_cluster_finder.finder(gns1_mpm['pm_x'], gns1_mpm['pm_y'],
#                                      gns1_mpm['xp'], gns1_mpm['yp'], 
#                                      gns1_mpm['l'].value, gns1_mpm['b'].value,
#                                      modes[0],
#                                      gns1_mpm['H'],gns1_mpm['H'],
#                                      knn,gen_sim,sim_lim, save_reg = pruebas1)
# # %%


# =============================================================================
# # Hyper-velocity star
# # =============================================================================

# hv = 25# We cosider hv if pm > 25 mas/yr (1000 km/s)

# hv_mask = np.sqrt(gns1_mi['pm_x']**2 + gns1_mi['pm_y']**2) > hv

# gns1_hv = gns1_mi[hv_mask]
# gns2_hv = gns2_mi[hv_mask]

# gns1_hv['H2'] = gns2_hv['H']
# gns1_hv.write(pruebas1 +  f'HV_gns1_pmSuper_F1_{field_one}_F2_{field_two}.ecsv', format = 'ascii.ecsv', overwrite = True)
# print(gns1_hv['pm_x', 'pm_y','H', 'H2'])


















































