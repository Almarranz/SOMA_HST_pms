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
from glob import glob
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
# results = ''
zone = 'G032.03+00.05'
# zone = 'G028.20-00.05'
zones = [ 'AFGL5180','G028.20-00.05', 'G032.03+00.05', 'G35.2-0.74N', 'G339.88-01.26','IRAS07299-1651', 'IRAS16562-3959']
# zones = [ 'IRAS16562-3959']
# zones = ['G032.03+00.05']

band_ls = ['160w', '110w']
names_ls = ['H', 'J']

# band_ls = ['160w', '128n']
# names_ls = ['H', 'Paβ']






# =============================================================================
# HST observation
# =============================================================================

pixSca = 0.12825 #arcsec/pixel
# pixSca = 0.00001 #arcsec/pixel#


# e_pos_cat  = 0.05# in arcsec. The position errors from starfinder lists are largely overstatimated!!!


# =============================================================================
# Aligment paremeters
# =============================================================================
max_sep = 80*u.mas


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
sig_cl_H = 30
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
look_for_cluster = 'no'
# look_for_cluster = 'yes'


# center = SkyCoord(ra = 282.41000248, dec = -0.78212655, unit = 'degree', frame = 'icrs')


bad = []
lopping = 1
# for loop in range(1):
wloop_counter = 0
# ZP = 25 # INVENTED
# for zone in zones[4:5]:
for zone in zones:
    
    for epoch in range(1,3):
        
        # =============================================================================
        # Dictionaries
        # =============================================================================
        cat_dic = {}
        obst_dic = {}

        for band in band_ls:
            print(band)
            
            # zone = zones[1]
            
            if band == '160w':
                b_name = 'H'
            elif band == '110w':
                b_name = 'J'
            elif band == '128n':
                b_name = 'Paβ'
            elif band == '164N':
                b_name = 'Fe II'
        
            pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'
            tmp1 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch1/tmp/'
            tmp2 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch2/tmp/'
            pm_results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/'
        
           
                
            
            results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
            cat = Table.read(results + f'hst_{zone}_epoch{epoch}_stars_f{band}.txt', format = 'ascii')
            cat['id'] = np.arange(1,len(cat)+1)
        
            pattern = folder + f"{zone}/epoch{epoch}/hst_*f{band}*"
            
            ima_name = os.path.basename(glob(pattern)[0])
            ima = fits.open(glob(pattern)[0] + f'/{ima_name}_drz.fits')
            wcs = WCS(ima[1].header)
            scale_pix = np.sqrt(ima[1].header['CD1_1']**2 + ima[1].header['CD2_1']**2 )*3600
            
            print(pixSca, scale_pix)
        # =============================================================================
        #         ZP calculation
        # =============================================================================
            mjd = ima[0].header['EXPSTART']
            ZP = get_vegazp(mjd, band = f'f{band}')
            print(ZP)
        
        
            cat[f'{b_name}'] =  (-2.5*np.log10(cat['f']) + ZP).round(4)
            cat[f'd{b_name}'] = ((2.5 / np.log(10)) * (cat['sf'] / cat['f'])).round(4)
            
            
            print(30*"=")
            print(results)
            print(30*"=")
            
            cat_rd = wcs.pixel_to_world(cat['x'], cat['y'])
            
            
            cat['ra'] = cat_rd.ra
            cat['dec'] = cat_rd.dec
            
           
            cat['x'] = cat['x']*pixSca # These are arcsec
            cat['y'] = cat['y']*pixSca
            
            cat['sx'] = cat['sx']*pixSca # These are arcsec
            cat['sy'] = cat['sy']*pixSca
            
            cat['sxy'] = np.sqrt(cat['sx']**2 + cat['sy']**2)
            # cat = cat[cat['sxy'] > 0]
            
            cat.write(results + f'calib_{zone}_EP{epoch}_f{band}.txt', format = 'ascii', overwrite = True )

            
            
            fig, (ax, ax2) = plt.subplots(1,2, figsize = (8,4))
            ax2.set_title(f'Epoch{epoch} ZP = {ZP: .3f}')
            ax.set_title(f'{zone} {b_name}')
            h = ax.hexbin(cat[f'{b_name}'], cat['sxy'], norm = LogNorm())
            cbar = plt.colorbar(h, ax=ax, label = f'\u2605 ={len(cat)}')
            ax.set_ylabel('$\sqrt{\sigma{x}^2 + \sigma{y}^2}$ [arcsec]', fontsize = 12)
            ax.set_xlabel(f'[{b_name}]', fontsize = 12)
            # ax.axhline(e_pos_cat, ls = 'dashed', color = 'red')
            
            h2 = ax2.hexbin(cat[f'{b_name}'], cat[f'd{b_name}'], norm = LogNorm(), cmap = 'inferno')
            cbar = plt.colorbar(h2, ax=ax2, label = f'\u2605 ={len(cat)}' )  
            ax2.set_ylabel(f'[d{b_name}]', fontsize = 12)
            ax2.set_xlabel(f'[{b_name}]', fontsize = 12)
            # ax2.axhline(e_pos_cat, ls = 'dashed', color = 'red')
            fig.tight_layout()
            mjd = ima[0].header['EXPSTART']   # e.g., 57610.4918468
            obst = Time(mjd, format='mjd', scale='utc')
            obst_dic[f'tyr{epoch}'] = obst.decimalyear
            
        
            # cat = filter_hst_data(cat, max_e_pos = e_pos_cat)
           
            # obst_dic[f't{epoch}'] = obst
            # cat_dic[f'cat{epoch}'] = cat
        
            
            
            if band == '160w':
                center = SkyCoord(ra = np.mean(cat_rd.ra.value), dec = np.mean(cat_rd.dec.value),
                                  unit = 'degree', frame = 'icrs')
           
            
            xg_1, yg_1 = center.spherical_offsets_to(cat_rd)
            
            
            tag = center.skyoffset_frame()
            
            cat_t = cat_rd.transform_to(tag)
            
            
            
            cat['xp'] = cat_t.lon.to(u.arcsec)
            cat['yp'] = cat_t.lat.to(u.arcsec)
            
            
            cat_dic[f'{b_name}'] = cat
        
# %%
            

        cat1 = cat_dic[f'{names_ls[0]}']
        cat1[f'{names_ls[1]}'] = 99.9999
        cat1[f'd{names_ls[1]}'] = 99.9999
        
        cat2 = cat_dic[f'{names_ls[1]}']
        
        crd1 = SkyCoord(ra = cat1['ra'], dec = cat1['dec'], unit = 'degree', frame = 'fk5')
        crd2 = SkyCoord(ra = cat2['ra'], dec = cat2['dec'], unit = 'degree', frame = 'fk5')
        
        
        idx, d2d, _ = crd1.match_to_catalog_sky(crd2, nthneighbor=1)
        match_mask = d2d < max_sep
        cat1_m = cat1[match_mask]
        cat2_m = cat2[idx[match_mask]]
        print(40*'+')
        unicos = unique(cat1_m, keep = 'first')
        print(len(cat1_m),len(unicos))
        print(40*'+')
        cat1[f'{names_ls[1]}'][match_mask] = cat2[idx[match_mask]][f'{names_ls[1]}']
        cat1[f'd{names_ls[1]}'][match_mask] = cat2[idx[match_mask]][f'd{names_ls[1]}']
        
        
        
        color_l = 2
        b_name1 = f'{names_ls[0]}'
        b_name2 = f'{names_ls[1]}'
        fig, ax = plt.subplots(1,1, figsize =(4,4)) 
        ax.set_title(f'{zone} \u2605 = {len(cat1_m)}')
        # h = ax.hexbin(cat2_m[b_name2] - cat1_m[b_name1], cat1_m[b_name1], norm = LogNorm(),cmap = 'viridis' )
        h = ax.hexbin(cat1[f'{names_ls[1]}'][cat1[f'{names_ls[1]}']<99] - cat1[f'{names_ls[0]}'][cat1[f'{names_ls[1]}']<99], cat1[f'{names_ls[0]}'][cat1[f'{names_ls[1]}']<99], norm = LogNorm(),cmap = 'viridis' )
        cbar = plt.colorbar(h, ax=ax, aspect = 40)
        ax.axvline(color_l, color = 'red', ls = 'dashed',label = f'{b_name2} - {b_name1} = {color_l}')
        ax.invert_yaxis()
        ax.set_xlabel(f'{b_name2} - {b_name1}')
        ax.set_ylabel(f'{b_name1}')
        ax.legend()
        
        c1_m = np.array([cat1_m['xp'],cat1_m['yp']]).T
        c2_m = np.array([cat2_m['xp'],cat2_m['yp']]).T
        p = ski.transform.estimate_transform('similarity',c2_m, c1_m)
        print(p)
        
        print("Translation: xp, yp [mas]  = (%.4f, %.4f)"%(p.translation[0]*1000,p.translation[1]*1000))
        print("Rotation: %.2f deg"%(p.rotation * 180.0/np.pi)) 
        print("Rotation: %.0f arcmin"%(p.rotation * 180.0/np.pi*60)) 
        print("Rotation: %.3f arcsec"%(p.rotation * 180.0/np.pi*3600)) 
        
        comb_f = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band_ls[0]}/epoch{epoch}/'
        cat1.write(comb_f + f'hst_{zone}_epoch{epoch}_{names_ls[0]}{names_ls[1]}.txt', format = 'ascii', overwrite = True)
# %%

ep1_f = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band_ls[0]}/epoch1/'
ep1 = Table.read(ep1_f + f'hst_{zone}_epoch1_{names_ls[0]}{names_ls[1]}.txt', format = 'ascii')
ep2_f = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band_ls[0]}/epoch2/'
ep2 = Table.read(ep2_f + f'hst_{zone}_epoch2_{names_ls[0]}{names_ls[1]}.txt', format = 'ascii')


fig, ax = plt.subplots(1,1)
ax.scatter(ep1['xp'], ep1['yp'],)
ax.scatter(ep2['xp'], ep2['yp'],s=0.2)
# ax.scatter(ep1['ra'], ep1['dec'],)
# ax.scatter(ep2['ra'], ep2['dec'],s=0.2)























    