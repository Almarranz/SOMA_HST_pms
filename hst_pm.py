#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 12:45:15 2026

@author: amartinez
"""

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
from filters import filter_gns_data
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
folder = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'
zone = 'G032.03+00.05'
band = '160w'
# epoch = 2

pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'


# =============================================================================
# Gaia parametres
# =============================================================================
radius = 200*u.arcsec
max_sep = 50*u.mas
mag_gaia = [12,18]
e_pm_gaia = 0.5
e_pos_gaia = 1
# =============================================================================
# HST observation
# =============================================================================

pixSca = 0.12825 #arcsec/pixel
# pixSca = 0.00001 #arcsec/pixel

# =============================================================================
# Aligment paremeters
# =============================================================================
# pre_transf = 'similarity'
transf = 'polynomial'
# transf = 'affine'
# transf = 'similarity'
order_trans = 1
align_loop = 'yes'
# align_loop = 'no'
align = 'Polywarp'
max_deg = 3# If this is <2 it does not enter the alignment loop. 
max_loop = 10

# =============================================================================
# Dictionaries
# =============================================================================
cat_dic = {}
obst_dic = {}
gaia_dic = {}
center = SkyCoord(ra = 282.41000248, dec = -0.78212655, unit = 'degree', frame = 'icrs')
for epoch in range(1,3):
    cat = Table.read(folder + f'{zone}/gaia_alignment/Epoch{epoch}/starfinder/{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii')
    # ima = fits.open(f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_{band}_drz_sci.fits')
    ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/G032.03+00.05_EP{epoch}_f{band}_drz_sci.fits')
    wcs = WCS(ima[0].header)
    
    cat_rd = wcs.pixel_to_world(cat['x'], cat['y'])
    
    # center = SkyCoord(ra = np.mean(cat_rd.ra.value), dec = np.mean(cat_rd.dec.value),
    #                   unit = 'degree', frame = 'icrs')
    
    cat['ra'] = cat_rd.ra
    cat['dec'] = cat_rd.dec
    
    cat['x'] = cat['x']*pixSca # These are arcsec
    cat['y'] = cat['y']*pixSca
    
    mjd = ima[0].header['EXPSTART']   # e.g., 57610.4918468
    obst = Time(mjd, format='mjd', scale='utc')
    obst_dic[f't{epoch}'] = obst.decimalyear
    
    # plt.scatter(ra_dec.ra, ra_dec.dec)
    # plt.scatter(center.ra, center.dec)
    # stop()
    
    
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
        
    gaia['id'] = np.arange(len(gaia))
    
    
    fig, ax2 = plt.subplots(1,1)
    ax2.scatter(gaia['phot_g_mean_mag'],gaia['pmra_error'], s= 2, label = 'Gaia $\delta \mu_{ra}$')
    ax2.scatter(gaia['phot_g_mean_mag'],gaia['pmdec_error'], s= 2, label = 'Gaia $\delta \mu_{dec}$')
    ax2.axvline(mag_gaia[0], color = 'r', ls = 'dashed', label = 'cuts')
    ax2.axvline(mag_gaia[1], color = 'r', ls = 'dashed', label = 'cuts')
    ax2.axhline(e_pm_gaia, color = 'r', ls = 'dashed')
    
    ax2.set_xlabel('[G]')
    ax2.set_ylabel('$\delta \mu$ [mas/yr]')
    
    fig.tight_layout(pad=1.0)
    lg = ax2.legend(markerscale=3.0)
    
    
    
    
    gaia = filter_gaia_data(
        gaia_table=gaia,
        astrometric_params_solved=31,
        duplicated_source= False,
        parallax_over_error_min=-10,
        astrometric_excess_noise_sig_max=2,
        phot_g_mean_mag_min = mag_gaia[1],
        phot_g_mean_mag_max = mag_gaia[0],
        pm_min=None,
        pmra_error_max= e_pm_gaia,
        pmdec_error_max= e_pm_gaia,
        ra_error_max= e_pos_gaia,
        dec_error_max= e_pos_gaia,
        ruwe = 1.4
        )
    
    # =============================================================================
    # This elimanates random Gaia stars
    # =============================================================================
    # =============================================================================
    #     cut_gaia = int(len(gaia)*0.70)
    #     idx_cut= np.random.choice(len(gaia), size=cut_gaia, replace=False)
    # 
    #     # Subtable with the same structure and column names
    #     gaia= gaia[idx_cut]
    # =============================================================================
    # =============================================================================
    # This elimanates random Gaia stars
    # =============================================================================
      
    gaia_rd = SkyCoord(ra = gaia['ra'], dec = gaia['dec'],
                       pm_ra_cosdec = gaia['pmra'],
                       pm_dec = gaia['pmdec'], 
                       frame = 'icrs', obstime = 'J2016')
    
    gaia_rdt = gaia_rd.apply_space_motion(new_obstime = obst)
    
    
    xp_g, yp_g = center.spherical_offsets_to(gaia_rdt.frame)
    
    gaia['rdt'] = gaia_rdt
    
    gaia['x'] = xp_g.to(u.arcsec)
    gaia['y'] = yp_g.to(u.arcsec)
    
    # gaia['pm_l'] = gaia_lb.pm_l_cosb
    # gaia['pm_b'] = gaia_lb.pm_b
    
    tg = Time(['2016-01-01T00:00:00'],scale='utc')
    
    fig, (ax, ax2) = plt.subplots(1,2)
    ax.scatter(cat['x'], cat['y'])
    ax2.scatter(gaia['x'], gaia['y'])
    
    # stop(211)
    
    gaia_dic[f'gaia{epoch}'] = gaia
    # %
    
    def sig_cl(x, y,s):
        mx, lx, hx = sigma_clip(x , sigma = s, masked = True, return_bounds= True)
        my, ly, hy = sigma_clip(y , sigma = s, masked = True, return_bounds= True)
        m_xy = np.logical_and(np.logical_not(mx.mask),np.logical_not(my.mask))
        
        return m_xy, [lx,hx,ly,hy]
    
    idx, d2d, _ = gaia_rdt.match_to_catalog_sky(cat_rd, nthneighbor=1)
    match_mask = d2d < max_sep
    gaia_m = gaia[match_mask]
    cat_m = cat[idx[match_mask]]
    print(40*'+')
    unicos = unique(cat_m, keep = 'first')
    print(len(cat_m),len(unicos))
    print(40*'+')
    
    # %%
    
    dra = cat_m['ra'] - gaia_m['ra'] 
    ddec = cat_m['dec'] - gaia_m['dec'] 
    cat_m['diff_ra'] = dra.to(u.mas)
    cat_m['diff_dec'] = ddec.to(u.mas)
    
    plot_two_hists_sigma(cat_m,'diff_ra', 'diff_dec', 'dRA', 'dDec', bins = 'auto', sig_clip=3, title1 = 'dRA, dDec')
    
    m_rd, lim = sig_cl(dra, ddec, s=3)
    
    cat_m = cat_m[m_rd]
    gaia_m = gaia_m[m_rd]
    
    
    # Fit transform
    if transf == 'polynomial':
        p = ski.transform.estimate_transform(
            transf, np.array([cat_m['x'], cat_m['y']]).T,
            np.array([gaia_m['x'], gaia_m['y']]).T, order=order_trans
        )
    # elif transf == 'Weight':
    #     model_x = Polynomial2D(degree=order_trans)
    #     model_y = Polynomial2D(degree=order_trans)
        
    #     # Linear least-squares fitter
    #     fitter = LinearLSQFitter()
        
    #     fit_xw = fitter(model_x, gns_m['xp'], gns_m['yp'],  gaia_m['xp'], weights= 1/np.sqrt( gns_m['sl']**2 +  gns_m['sb']**2))  # Fit x-coordinates
    #     fit_yw = fitter(model_y, gns_m['xp'], gns_m['yp'],  gaia_m['yp'],weights= 1/np.sqrt( gns_m['sl']**2 +  gns_m['sb']**2)) 
        
        
    else:
        p = ski.transform.estimate_transform(
            transf, np.array([cat_m['x'], cat_m['y']]).T,
            np.array([gaia_m['x'], gaia_m['y']]).T
        )
       
    
    cat_trans = p(np.array([cat['x'], cat['y']]).T)
    cat['x'] = cat_trans[:, 0]
    cat['y'] = cat_trans[:, 1]
    
    cat_xy = np.array([cat['x'],cat['y']]).T
    gaia_xy = np.array([gaia['x'],gaia['y']]).T
    xy_mat = compare_lists(cat_xy, gaia_xy, max_sep.to(u.arcsec).value*2)
    
    
    cat_m = cat[xy_mat['ind_1']]
    gaia_m = gaia[xy_mat['ind_2']]
    
    
    dx = cat_m['x'] - gaia_m['x'] 
    dy = cat_m['y'] - gaia_m['y'] 
    cat_m['diff_x'] = dx*1000
    cat_m['diff_y'] = dy*1000
    
    plot_two_hists_sigma(cat_m,'diff_x', 'diff_y', 'dx', 'dy', bins = 'auto', sig_clip=3, title1 = 'dx, dy')
    
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax1.set_title(f'Matches = {len(xy_mat)}\nMin dist = {max_sep} ')
    # ax1.scatter(cat['x'][::100], cat['y'][::100], alpha =0.1, color = 'k')
    ax1.scatter(cat['x'], cat['y'], alpha =0.1, color = 'k')
    ax1.scatter(gaia['x'], gaia['y'],s =10, label = 'Gaia')
    ax1.scatter(cat['x'][xy_mat['ind_1']], cat['y'][xy_mat['ind_1']], label = 'GNS1 match')
    ax1.scatter(gaia['x'][xy_mat['ind_2']], gaia['y'][xy_mat['ind_2']],s =10, label = 'Gaia match')
    ax1.set_xlabel('x[arcsec]')
    ax1.set_xlabel('y[arcsec]')
    ax1.axis('equal')
    ax1.legend()
    
    
    if align_loop  == 'yes':
        # Optional final refinement
        # gns_cat = alg_rel(gns_cat, gaia, 'xp', 'yp', align, max_deg, max_sep.to(u.arcsec).value)
        cat = alig_gaia(cat, gaia, 'x', 'y', align, max_deg, max_sep.to(u.arcsec).value, max_loop)
      # alg_loop(gns_cat, gaia, 'xp', 'yp', align, max_deg, max_sep.to(u.arcsec).value, max_loop)
    elif align_loop == 'no':
        p2 = ski.transform.estimate_transform('polynomial', 
            np.array([cat['x'][xy_mat['ind_1']], cat['y'][xy_mat['ind_1']]]).T,
            np.array([gaia['x'][xy_mat['ind_2']], gaia['y'][xy_mat['ind_2']]]).T, order=2)
        
        cat_trans = p2(np.array([cat['x'], cat['y']]).T)
        
        cat['x'] = cat_trans[:,0]
        cat['y'] = cat_trans[:,1]
    
    
    
    cat_xy = np.array([cat['x'],cat['y']]).T
    gaia_xy = np.array([gaia['x'],gaia['y']]).T
    xy_mat = compare_lists(cat_xy, gaia_xy, max_sep.to(u.arcsec).value)
    
    ax2.set_title(f'Matches = {len(xy_mat)}\nMin dist = {max_sep} ')
    ax2.scatter(cat['x'][::100], cat['y'][::100], alpha =0.1, color = 'k')
    ax2.scatter(gaia['x'], gaia['y'],s =10, label = 'Gaia')
    ax2.scatter(cat['x'][xy_mat['ind_1']], cat['y'][xy_mat['ind_1']], label = 'GNS1 match')
    ax2.scatter(gaia['x'][xy_mat['ind_2']], gaia['y'][xy_mat['ind_2']],s =10, label = 'Gaia match')
    ax2.set_xlabel('x[arcsec]')
    ax2.yaxis.set_label_position("right")
    ax2.set_ylabel('y [arcsec]')
    ax2.legend()
    fig.tight_layout()
    plt.show()
    
    
    
    cat_xy = np.array([cat['x'],cat['y']]).T
    gaia_xy = np.array([gaia['x'],gaia['y']]).T
    xy_mat = compare_lists(cat_xy, gaia_xy, max_sep.to(u.arcsec).value*2)
    
    
    cat_m = cat[xy_mat['ind_1']]
    gaia_m = gaia[xy_mat['ind_2']]
    
    
    dx = cat_m['x'] - gaia_m['x'] 
    dy = cat_m['y'] - gaia_m['y'] 
    cat_m['diff_x'] = dx*1000
    cat_m['diff_y'] = dy*1000
    
    plot_two_hists_sigma(cat_m,'diff_x', 'diff_y', 'dx', 'dy', bins = 'auto', sig_clip=3, title1 = f'After alignment Epoch {epoch} ')
    
    cat_dic[f'cat{epoch}'] = cat
    
cat1 = cat_dic['cat1']
cat2 = cat_dic['cat2']

cat1_xy = np.array([cat1['x'],cat1['y']]).T
cat2_xy = np.array([cat2['x'],cat2['y']]).T
# %%


xy_cat = compare_lists(cat1_xy, cat2_xy, 0.150)

cat1_c = cat1[xy_cat['ind_1']]
cat2_c = cat2[xy_cat['ind_2']]

dt = obst_dic['t2'] - obst_dic['t1']
pm_x = ((cat2_c['x']*u.arcsec - cat1_c['x']*u.arcsec)/(dt*u.yr)).to(u.mas/u.yr)
pm_y = ((cat2_c['y']*u.arcsec - cat1_c['y']*u.arcsec)/(dt*u.yr)).to(u.mas/u.yr)

cat1_c['pmx'] = pm_x
cat1_c['pmy'] = pm_y

plot_two_pm_hists(cat1_c, 'pmx', 'pmy', r'\mu_{x}', r'\mu_{y}')



gaia1 = gaia_dic['gaia1']
plot_two_pm_hists(gaia1, 'pmra', 'pmdec', r'\mu_{RA}', r'\mu_{Dec}', title1 = 'Gaia pm')

cat1_rd = SkyCoord(ra = cat1_c['ra'], dec = cat1_c['dec'], unit = 'degree', frame = 'icrs', obstime = f'J{obst_dic["t1"]}')
idx, d2d, _ = gaia1['rdt'] .match_to_catalog_sky(cat1_rd, nthneighbor=1)
match_mask = d2d < max_sep*3
gaia1_m = gaia1[match_mask]
cat1_m = cat1_c[idx[match_mask]]

dpmx = gaia1_m['pmra'] - cat1_m['pmx']
dpmy = gaia1_m['pmdec'] - cat1_m['pmy']

cat1_m['dpmx'] = dpmx
cat1_m['dpmy'] = dpmy


# plot_two_hists_sigma(cat1_m, 'dpmx', 'dpmy', r'\Delta\,pm RA', r'\Delta\,pm Dec', bins = 'auto', title1 = f'Matches = {len(cat1_m)}')
plot_two_hists_sigma(cat1_m, 'dpmx', 'dpmy', r'\Delta\,pm RA', r'\Delta\,pm Dec', bins = 'auto', title1 = f'Matches = {len(cat1_m)}')







