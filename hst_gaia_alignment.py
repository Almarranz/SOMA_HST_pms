
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 10:30:27 2022

@author: amartinez
"""

# Generates offsets for Gaia stars and astroaligned with xy coordinates in GNS1

# Here we are going to align GNS (1 and 2) to Gaia reference frame for the 
# each of tha gns epochs
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
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})
rcParams.update({'figure.figsize':(10,5)})
rcParams.update({
    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
# Enable automatic plotting mode
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')


folder = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'
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
mag_gaia = [12,18]
e_pm_gaia = 0.3
e_pos_gaia = 0.5
# =============================================================================
# HST observation
# =============================================================================

pixSca = 0.12825 #arcsec/pixel
# pixSca = 0.00001 #arcsec/pixel

# =============================================================================
# Aligment paremeters
# =============================================================================
# pre_transf = 'polynomial'
pre_transf = None

transf = 'polynomial'
# transf = 'affine'
# transf = 'similarity'
order_trans = 1
# align_loop = 'yes'
align_loop = 'no'
align = 'Polywarp'
max_deg = 3# If this is <2 it does not enter the alignment loop. 
max_loop = 3
# gaia_clipping = 'one_one'# Clipp the Gaia outlayer one by one
gaia_clipping = 'all'# Clipp the Gaia outlayer all at once

# =============================================================================
# Proper motions param
# ===========================================================================
max_dis_pm = 0.150#in arcsec
sig_H = 3# discrd pm for stars with delta H over sig_H
e_pm_cat = 1# im mas/yr


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
while lopping > 0:
    

    for epoch in range(1,3):
        cat = Table.read(folder + f'{zone}/gaia_alignment/Epoch{epoch}/starfinder/{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii')
        # ima = fits.open(f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_{band}_drz_sci.fits')
        ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_f{band}_drz_sci.fits')
        wcs = WCS(ima[0].header)
        
        cat_rd = wcs.pixel_to_world(cat['x'], cat['y'])
        
        if epoch == 1:
            center = SkyCoord(ra = np.mean(cat_rd.ra.value), dec = np.mean(cat_rd.dec.value),
                              unit = 'degree', frame = 'icrs')
        
        cat['ra'] = cat_rd.ra
        cat['dec'] = cat_rd.dec
        
        cat['x'] = cat['x']*pixSca # These are arcsec
        cat['y'] = cat['y']*pixSca
        
        cat['sx'] = cat['sx']*pixSca # These are arcsec
        cat['sy'] = cat['sy']*pixSca
        
        mjd = ima[0].header['EXPSTART']   # e.g., 57610.4918468
        obst = Time(mjd, format='mjd', scale='utc')
        # obst_dic[f't{epoch}'] = obst.decimalyear
        obst_dic[f't{epoch}'] = obst
    
        cat_dic[f'cat{epoch}'] = cat
        
        region(cat, 'ra', 'dec',
               name = f'{zone}_Ep{epoch}_stars',
               save_in = pruebas)
    
    stop(179)
    # m_mask1 = (gns1['H'] > m_lim[0]) & (gns1['H'] < m_lim[1])
    # gns1 = gns1[m_mask1]
    
    # unc_cut1 = (gns1['sl']< max_sig) & (gns1['sb'] < max_sig)
    # gns1 = gns1[unc_cut1]
    
    # m_mask2 = (gns2['H'] > m_lim[0]) & (gns2['H'] < m_lim[1])
    # gns2 = gns2[m_mask2]
    
    
    # unc_cut2 = (gns2['sl']<max_sig) & (gns2['sb']<max_sig)
    # gns2 = gns2[unc_cut2]
    
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
       

   
    

     
    if len(bad)>0:#!!!
        del_1 = np.isin(gaia['id'], bad)#!!!
        gaia = gaia[np.logical_not(del_1)]#!!!
    
    
    
    
    fig, ax2 = plt.subplots(1,1)
    ax2.scatter(gaia['phot_g_mean_mag'],gaia['pmra_error'], s= 2, label = 'Gaia $\delta \mu_{ra}$')
    ax2.scatter(gaia['phot_g_mean_mag'],gaia['pmdec_error'], s= 2, label = 'Gaia $\delta \mu_{dec}$')
    ax2.axvline(mag_gaia[0], color = 'r', ls = 'dashed', label = 'quality cuts')
    ax2.axvline(mag_gaia[1], color = 'r', ls = 'dashed')
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

    
    
  
    gaia_rd = SkyCoord(ra = gaia['ra'], dec = gaia['dec'],
                       pm_ra_cosdec = gaia['pmra'],
                       pm_dec = gaia['pmdec'], 
                       frame = 'icrs', obstime = 'J2016')
    
   
    tg = Time(['2016-01-01T00:00:00'],scale='utc')
    
    
    
    # %
    
    def sig_cl(x, y,s):
        mx, lx, hx = sigma_clip(x , sigma = s, masked = True, return_bounds= True)
        my, ly, hy = sigma_clip(y , sigma = s, masked = True, return_bounds= True)
        m_xy = np.logical_and(np.logical_not(mx.mask),np.logical_not(my.mask))
        
        return m_xy, [lx,hx,ly,hy]
    
    # Define catalogs and times as a list of tuplesçç
    cat1 = cat_dic['cat1']
    t1_cat = obst_dic['t1']
    
    cat2 = cat_dic['cat2']
    t2_cat = obst_dic['t2']
    
    catalogs = [
        {'name': 'cat1', 'cat': cat1, 'time': t1_cat, 'tag': '1'},
        {'name': 'cat2', 'cat': cat2, 'time': t2_cat, 'tag': '2'}
    ]
    
    for c,cats in enumerate(catalogs):
        # max_sep = max_sep
        print(f"\n===== Aligning {cats['name']} =====")
    
       
        # Recalculate Gaia projected positions at catalog's epoch
        
    
        # gaia_lb = SkyCoord(
        #     ra=gaia['ra'], dec=gaia['dec'],
        #     pm_ra_cosdec=gaia['pmra'], pm_dec=gaia['pmdec'],
        #     frame='icrs', obstime='J2016'
        # ).galactic
        
      
        # gaia_l = gaia_lb.l + gaia_lb.pm_l_cosb*dt.to(u.yr)
        # gaia_b = gaia_lb.b + gaia_lb.pm_b*dt.to(u.yr)
        
        # gaia_lbt = SkyCoord(l = gaia_l, b = gaia_b, frame = 'galactic')
        
        gaia_rd = SkyCoord(
            ra=gaia['ra'], dec=gaia['dec'],
            pm_ra_cosdec=gaia['pmra'], pm_dec=gaia['pmdec'],
            frame='icrs', obstime='J2016'
        )
        
        gaia_rdt = gaia_rd.apply_space_motion(new_obstime = cats['time'])
        
        gaia[f'ra{cats["tag"]}'] = gaia_rdt.ra
        gaia[f'dec{cats["tag"]}'] = gaia_rdt.dec
        
        xp_g, yp_g = center.spherical_offsets_to(gaia_rdt.frame)
        
        
        
        gaia['x'] = xp_g.to(u.arcsec) 
        gaia['y'] = yp_g.to(u.arcsec) 
        
        gaia[f'x_{c+1}'] = xp_g.to(u.arcsec) 
        gaia[f'y_{c+1}'] = yp_g.to(u.arcsec) 
         
    
        hst_cat = cats['cat']
        gaia_c = SkyCoord(ra = gaia[f'ra{cats["tag"]}'], dec = gaia[f'dec{cats["tag"]}'], frame='icrs')
        cat_c = SkyCoord(ra = hst_cat['ra'], dec = hst_cat['dec'], frame='icrs')
        
        
# =============================================================================
#         # This proyect Gaia proper motions on the same tangentical plas of that of the proyected x,y coordenates
#         # The projected pm are exactly the same that the Gaia propermotions...
#         offset_frame = center.skyoffset_frame()
# 
#         c_proj = gaia_rdt.transform_to(offset_frame)
# 
#         # Step 3: Extract Gaia PM in tangent plane (same as your XY frame, in mas/yr)
#         pm_x_gaia = c_proj.pm_lon_coslat  # mas/yr
#         pm_y_gaia = c_proj.pm_lat        # mas/yr
# 
#         gaia['pm_x'] = pm_x_gaia
#         gaia['pm_y'] = pm_y_gaia
# =============================================================================
        
        
        
    
# =============================================================================
#         idx1, idx2, sep2d, _ = search_around_sky(gaia_c, cat_c, max_sep)
# 
#         count1 = Counter(idx1)
#         count2 = Counter(idx2)
# 
#         # Step 3: Create mask for one-to-one matches only
#         mask_unique = np.array([
#             count1[i1] == 1 and count2[i2] == 1
#             for i1, i2 in zip(idx1, idx2)
#         ])
# 
#         # Step 4: Apply the mask
#         idx1_clean = idx1[mask_unique]
#         idx2_clean = idx2[mask_unique]
# 
#         gaia_m= gaia[idx1_clean]
#         cat_m = hst_cat[idx2_clean]
#         
#        
#         
#         
#         print(40*'+')
#         unicos = unique(cat_m, keep = 'first')
#         print(len(cat_m),len(unicos))
#         print(40*'+')
#         
# =============================================================================
        idx, d2d, _ = gaia_c.match_to_catalog_sky(cat_c, nthneighbor=1)
        match_mask = d2d < max_sep
        gaia_m = gaia[match_mask]
        cat_m = hst_cat[idx[match_mask]]
        print(40*'+')
        unicos = unique(cat_m, keep = 'first')
        print(len(cat_m),len(unicos))
        print(40*'+')
    
       
        
        fig, (ax, ax1, ax2) = plt.subplots(1,3, figsize = (15,5))
        fig.suptitle(f'GNS{c+1}')
        # ax.scatter(gns_w[::100], gns_cat['b'][::100], alpha =0.1, color = 'k')
        # ax.scatter(ga_w, gaia['b'], label = 'Gaia', s= 10)
        ax.scatter(hst_cat['ra'][::10], hst_cat['dec'][::10], alpha =0.1, color = 'k')
        ax.set_title(f'Matches = {len(cat_m)}\nMin dist = {max_sep} ')
        ax.scatter(cat_m['ra'], cat_m['dec'], label = f'cat{c} Match')
        ax.scatter(gaia_m['ra'], gaia_m['dec'],s =10, label = 'Gaia Match')
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        ax.legend()
        
# =============================================================================
#         # Apply a similiraty trasnformation fisrt
#         Does not make much of a difference
# =============================================================================
        
        if pre_transf is not None:
            print(30*'-'+'\n'+f'Applyin a pre transformation of {pre_transf}\n'+ 30*'-')
            if pre_transf == 'polynomial':
                psim = ski.transform.estimate_transform(
                pre_transf, np.array([cat_m['ra'], cat_m['dec']]).T,
                    np.array([gaia_m['ra'], gaia_m['dec']]).T, order = 1)
                
    
            else:
                psim = ski.transform.estimate_transform(
                pre_transf, np.array([cat_m['ra'], cat_m['dec']]).T,
                    np.array([gaia_m['ra'], gaia_m['dec']]).T)
            
            # if pre_transf == 'affine':
            #     psim = ski.transform.estimate_transform(
            #     pre_transf, np.array([gns_m['l'], gns_m['b']]).T,
            #         np.array([gaia_m['l'], gaia_m['b']]).T)
            
            
            cat_sim = psim(np.array([hst_cat['ra'], hst_cat['dec']]).T)
            hst_cat['ra'] = cat_sim[:, 0]*u.deg
            hst_cat['dec'] = cat_sim[:, 1]*u.deg



        gaia_c = SkyCoord(ra = gaia[f'ra{cats["tag"]}'], dec = gaia[f'dec{cats["tag"]}'], frame='icrs')
        gns_c = SkyCoord(ra = hst_cat['ra'], dec = hst_cat['dec'], frame='icrs')
        
        idx1, idx2, sep2d, _ = search_around_sky(gaia_c, gns_c, max_sep)

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

        gaia_m= gaia[idx1_clean]
        cat_m = hst_cat[idx2_clean]
        
        # Number of rows to select
       


        print(40*'+')
        unicos = unique(cat_m, keep = 'first')
        print(len(cat_m),len(unicos))
        print(40*'+')
        
# %%
        

# =============================================================================
#         This are the position residuals BEFORE any transformation
# =============================================================================
        if c == 0:        
            rcParams.update({
            "figure.figsize": (10, 5),
            "font.size": 18,
            "axes.labelsize": 18,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 16
        })
            dl_pre = cat_m['ra'].to(u.mas) - gaia_m['ra'].to(u.mas)
            db_pre = cat_m['dec'].to(u.mas) - gaia_m['dec'].to(u.mas)
            m_pre, pre_lims = sig_cl(dl_pre, db_pre, 3)
            dl_prem = dl_pre[m_pre]
            db_prem = db_pre[m_pre]
            
            fig_pre, (ax_pre, ax1_pre) = plt.subplots(1,2,figsize = (11,5.5))
            #
            ax_pre.set_title(f'Gaia vs HST{c+1} (Stars = {len(gaia_m)})')
            ax1_pre.set_title(f'Residuas before any transformation')
            # ax1.set_title(f'Matching stars  = {len(d_xm)}')
            ax_pre.set_ylabel('# stars')
            ax_pre.hist(dl_pre, bins = 10,  color = 'grey', alpha = 0.5)
            ax1_pre.hist(db_pre, bins = 'auto',  color = 'grey', alpha = 0.5)
            ax_pre.hist(dl_prem, bins = 10, histtype = 'step',lw = 2,label = '$\overline{\Delta RA}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(dl_prem.value),np.std(dl_prem.value)))
            ax1_pre.hist(db_prem,bins = 'auto',histtype = 'step', lw = 2, label = '$\overline{\Delta Dec}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(db_prem.value),np.std(db_prem.value)))
            ax_pre.legend(loc = 1)
            ax1_pre.legend()
            ax_pre.set_xlabel('$\Delta$RA [mas]')
            ax1_pre.set_xlabel('$\Delta$Dec [mas]')
            ax_pre.axvline(pre_lims[0].value, color = 'r', ls = 'dashed')
            ax_pre.axvline(pre_lims[1].value, color = 'r', ls = 'dashed')
            ax1_pre.axvline(pre_lims[2].value, color = 'r', ls = 'dashed')
            ax1_pre.axvline(pre_lims[3].value, color = 'r', ls = 'dashed')
            # ax_pre.set_xlim(-20,20)
            # ax1_pre.set_xlim(-20,20)
            fig_pre.tight_layout()
            
            meta = {'Script': '/Users/amartinez/Desktop/PhD/HAWK/GNS_pm_scripts/GNS_pm_absolute_SUPER/gns_gaia_alignment.py'}
            # plt.savefig(f'/Users/amartinez/Desktop/PhD/My_papers/GNS_pm_catalog/images/ABS_F1_{field_one}_gaia_resi_pos_prealign.png', dpi = 150, transparent = True, metadata = meta)
        
        # continue
        
# %%
        

        # Fit transform
        if transf == 'polynomial':
            p = ski.transform.estimate_transform(
                transf, np.array([cat_m['x'], cat_m['y']]).T,
                np.array([gaia_m['x'], gaia_m['y']]).T, order=order_trans
            )
            
        else:
            p = ski.transform.estimate_transform(
                transf, np.array([cat_m['x'], cat_m['y']]).T,
                np.array([gaia_m['x'], gaia_m['y']]).T
            )
       
    
        
        
        cat_trans = p(np.array([hst_cat['x'], hst_cat['y']]).T)
        hst_cat['x'] = cat_trans[:, 0]
        hst_cat['y'] = cat_trans[:, 1]
        print(p.params)
        
        cat_xy = np.array([hst_cat['x'],hst_cat['y']]).T
        gaia_xy = np.array([gaia['x'],gaia['y']]).T
        xy_mat = compare_lists(cat_xy, gaia_xy, max_sep.to(u.arcsec).value)
        
        # for i in range(15):
            
            
        #     # p = ski.transform.estimate_transform(
        #     #     transf, np.array([gns_cat['xp'][xy_mat['ind_1']] , gns_cat['yp'][xy_mat['ind_1']] ]).T,
        #     #     np.array([gaia['xp'][xy_mat['ind_2']] , gaia['yp'][xy_mat['ind_2']] ]).T)
            
        #     p = ski.transform.estimate_transform(
        #         'polynomial', np.array([gns_cat['xp'][xy_mat['ind_1']] , gns_cat['yp'][xy_mat['ind_1']] ]).T,
        #         np.array([gaia['xp'][xy_mat['ind_2']] , gaia['yp'][xy_mat['ind_2']] ]).T, order = 2)
            
        #     gns_trans = p(np.array([gns_cat['xp'], gns_cat['yp']]).T)
        #     gns_cat['xp'] = gns_trans[:, 0]
        #     gns_cat['yp'] = gns_trans[:, 1]
        
        #     gns_xy = np.array([gns_cat['xp'],gns_cat['yp']]).T
        #     gaia_xy = np.array([gaia['xp'],gaia['yp']]).T
        #     xy_mat = compare_lists(gns_xy, gaia_xy, max_sep.to(u.arcsec).value)
            
        #     print(len(xy_mat))
        #     print(p.params)
    
            
    
        # sys.exit()
        ax1.set_title(f'Matches = {len(xy_mat)}\nMin dist = {max_sep} ')
        ax1.scatter(hst_cat['x'][::100], hst_cat['y'][::100], alpha =0.1, color = 'k')
        ax1.scatter(gaia['x'], gaia['y'],s =10, label = 'Gaia')
        ax1.scatter(hst_cat['x'][xy_mat['ind_1']], hst_cat['y'][xy_mat['ind_1']], label = f'cat{c+1} match')
        ax1.scatter(gaia['x'][xy_mat['ind_2']], gaia['y'][xy_mat['ind_2']],s =10, label = 'Gaia match')
        ax1.set_xlabel('xp[arcsec]')
        ax1.legend()
        
        
        
        if align_loop  == 'yes':
            # Optional final refinement
            
            # hst_cat = alg_loop(hst_cat, gaia, 'x', 'y', align, max_deg, max_sep.to(u.arcsec).value, max_loop)
            hst_cat = alig_gaia(hst_cat, gaia, 'x', 'y', align, max_deg, max_sep.to(u.arcsec).value, max_loop)
       
        elif align_loop == 'no':
            p2 = ski.transform.estimate_transform('polynomial', 
                np.array([hst_cat['x'][xy_mat['ind_1']], hst_cat['y'][xy_mat['ind_1']]]).T,
                np.array([gaia['x'][xy_mat['ind_2']], gaia['y'][xy_mat['ind_2']]]).T, order=2)
            
            hst_trans = p2(np.array([hst_cat['x'], hst_cat['y']]).T)
            
            hst_cat['x'] = hst_trans[:,0]
            hst_cat['y'] = hst_trans[:,1]
        
        
        
        hst_xy = np.array([hst_cat['x'],hst_cat['y']]).T
        gaia_xy = np.array([gaia['x'],gaia['y']]).T
        xy_mat = compare_lists(hst_xy, gaia_xy, max_sep.to(u.arcsec).value)
        
        ax2.set_title(f'Matches = {len(xy_mat)}\nMin dist = {max_sep} ')
        ax2.scatter(hst_cat['x'][::100], hst_cat['y'][::100], alpha =0.1, color = 'k')
        ax2.scatter(gaia['x'], gaia['y'],s =10, label = 'Gaia')
        ax2.scatter(hst_cat['x'][xy_mat['ind_1']], hst_cat['y'][xy_mat['ind_1']], label = 'GNS1 match')
        ax2.scatter(gaia['x'][xy_mat['ind_2']], gaia['y'][xy_mat['ind_2']],s =10, label = 'Gaia match')
        ax2.set_xlabel('xp[arcsec]')
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel('yp [arcsec]')
        ax2.legend()
        fig.tight_layout()
        plt.show()
        # Residuals
        gns_xy = np.array([hst_cat['x'], hst_cat['y']]).T
        gaia_xy = np.array([gaia['x'], gaia['y']]).T
        xy_al = compare_lists(hst_xy, gaia_xy, max_sep.to(u.arcsec).value)
    
        d_x = (gaia['x'][xy_al['ind_2']] - hst_cat['x'][xy_al['ind_1']]).to(u.mas) 
        d_y = (gaia['y'][xy_al['ind_2']] - hst_cat['y'][xy_al['ind_1']] ).to(u.mas)
        
        sig_pm = 3
        # m_dx, l_dx, h_dx = sigma_clip(d_x, sigma = sig_pm, masked = True, return_bounds= True)
        # m_dy, l_dy, h_dy = sigma_clip(d_y, sigma = sig_pm, masked = True, return_bounds= True)
        # m_dxy = np.logical_and(np.logical_not(m_dx.mask),np.logical_not(m_dy.mask))
        
        m_dxy, lims = sig_cl(d_x,d_y, sig_pm)
        
        d_xm = d_x[m_dxy]
        d_ym = d_y[m_dxy]
        
        fig, (ax, ax1) = plt.subplots(1,2)
        # ax.set_title(f'Gaia vs {cat['name']} (proyected)')
        ax1.set_title(f'Matching stars  = {len(d_xm)}')
        ax.hist(d_x,  color = 'grey', alpha = 0.5)
        ax1.hist(d_y, color = 'grey', alpha = 0.5)
        ax.hist(d_xm, histtype = 'step',label = '$\overline{\Delta x}$ = %.2f mas\n$\sigma$ = %.2f'%(np.mean(d_xm.value),np.std(d_xm.value)))
        ax1.hist(d_ym,histtype = 'step', label = '$\overline{\Delta y}$ = %.2f mas\n$\sigma$ = %.2f'%(np.mean(d_ym.value),np.std(d_ym.value)))
        ax.legend(); ax1.legend()
        ax.set_xlabel('$\Delta$xp [mas]')
        ax1.set_xlabel('$\Delta$yp [mas]')
        ax.axvline(lims[0].value, color = 'r', ls = 'dashed')
        ax.axvline(lims[1].value, color = 'r', ls = 'dashed')
        ax1.axvline(lims[2].value, color = 'r', ls = 'dashed')
        ax1.axvline(lims[3].value, color = 'r', ls = 'dashed')
        
        cats['cat']['x'] = hst_cat['x']
        cats['cat']['y'] = hst_cat['y']
        plt.show()
    
    
    
    
        # print('Pasannado')
    # %%
    
    # =============================================================================
    #HST proper motions
    # =============================================================================
    cat1_al = catalogs[0]['cat']
    cat2_al = catalogs[1]['cat']
    
    dt_cat = catalogs[1]['time'].decimalyear - catalogs[0]['time'].decimalyear 
    
    
    cat1_gxy  = np.array([cat1_al['x'], cat1_al['y']]).T 
    cat2_gxy  = np.array([cat2_al['x'], cat2_al['y']]).T 
    
    cat_com = compare_lists(cat1_gxy, cat2_gxy,max_dis_pm )
    
    cat1_gxy  = cat1_gxy[cat_com['ind_1']]  
    cat2_gxy  = cat2_gxy[cat_com['ind_2']]  
    
    cat1 = cat1_al[cat_com['ind_1']]
    cat2 = cat2_al[cat_com['ind_2']]
    
    # dH = gns1['H'] - gns2['H']
    # mask_H, lH, hH = sigma_clip(dH , sigma = sig_H, masked = True, return_bounds= True)
    
    # fig, ax = plt.subplots(1,1)
    # ax.hist(dH, bins = 'auto', label = '$\overline{\Delta H}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(dH), np.std(dH)))
    # ax.set_xlabel('$\Delta$H GNS')
    # ax.axvline(lH, color = 'red', ls = 'dashed', label = f'{sig_H}$\sigma$')
    # ax.axvline(hH, color = 'red', ls = 'dashed')
    # ax.legend()

    
    # sys.exit(517)
    
    # gns1 = gns1[np.logical_not(mask_H.mask)]
    # gns2 = gns2[np.logical_not(mask_H.mask)]
    # gns_com = gns_com[np.logical_not(mask_H.mask)]
    # gns1_gxy = gns1_gxy[np.logical_not(mask_H.mask)]
    # gns2_gxy = gns2_gxy[np.logical_not(mask_H.mask)]
    
    
    pm_x = (cat_com['l2_x']*u.arcsec - cat_com['l1_x']*u.arcsec).to(u.mas)/(dt_cat*u.yr)
    pm_y = (cat_com['l2_y']*u.arcsec - cat_com['l1_y']*u.arcsec).to(u.mas)/(dt_cat*u.yr)
    
    
   
    err1_x = (cat1['sx'])*u.arcsec
    err1_y = (cat1['sy'])*u.arcsec
    err2_x = (cat2['sx'])*u.arcsec
    err2_y = (cat2['sy'])*u.arcsec
    
    dpm_x = np.sqrt((err1_x.to(u.mas))**2 + (err2_x.to(u.mas))**2)/(dt_cat*u.yr)
    dpm_y = np.sqrt((err1_y.to(u.mas))**2 + (err2_y.to(u.mas))**2)/(dt_cat*u.yr)

    
    cat1['pm_x'] = pm_x
    cat1['pm_y'] = pm_y
    cat2['pm_x'] = pm_x
    cat2['pm_y'] = pm_y
    
    cat1['dpm_x'] = dpm_x
    cat1['dpm_y'] = dpm_y
    cat2['dpm_x'] = dpm_x
    cat2['dpm_y'] = dpm_y
    
    
    
    
    # cat1.meta['m_lim_gns'] = m_lim
    cat1.meta['max_dis_pm'] = max_dis_pm
    cat1.meta['e_pm_gaia'] = e_pm_gaia
    # cat1.meta['mag_min_gaia'] = mag_min_gaia
    
# =============================================================================
#     Save the catlogs with the pm
# =============================================================================
    # cat1.write(pruebas1 + f'gns1_pmSuper_F1_{field_one}_F2_{field_two}.ecvs',format = 'ascii.ecsv', overwrite = True)
    # cat2.write(pruebas2 + f'gns2_pmSuper_F1_{field_one}_F2_{field_two}.ecvs',format = 'ascii.ecsv', overwrite = True)
    
    


    
    # gns2.meta['f_mode'] = f_mode


    # stop()
    cat1 = filter_hst_data(cat1, max_e_pm = e_pm_cat)
    cat2 = filter_hst_data(cat2, max_e_pm = e_pm_cat)

    # 
    # %%
    fig, (ax,ax2) = plt.subplots(1,2)
    
    
    bins = 30
    ax.set_title(f'{zone}, f{band}')
    ax.hist(pm_x, bins = bins, color = 'grey', alpha = 0.3)
    ax2.hist(pm_y, bins = bins, color = 'grey', alpha = 0.3)
    
    ax.hist(cat1['pm_x'], bins = bins, histtype = 'step', label = '$\overline{\mu}_{x}$ = %.2f\n$\sigma$ =%.2f'%(np.mean(cat1['pm_x']),np.std(cat1['pm_x'])))
    # ax.axvline(lpm[0] , ls = 'dashed', color = 'r')
    # ax.axvline(lpm[1] , ls = 'dashed', color = 'r')
    
    ax2.hist(cat1['pm_y'], bins = bins, histtype = 'step', label = '$\overline{\mu}_{y}$ = %.2f\n$\sigma$ =%.2f'%(np.mean(cat1['pm_y']),np.std(cat1['pm_y'])))
    # ax2.axvline(lpm[2] , ls = 'dashed', color = 'r')
    # ax2.axvline(lpm[3] , ls = 'dashed', color = 'r')
    
    ax.legend()
    ax2.legend()
    
    ax.set_xlabel('$\mu_{xp}$ [mas/yr]')
    ax2.set_xlabel('$\mu_{yp}$ [mas/yr]')
    fig.tight_layout()
    plt.show()
    
    
    g_fac = 1# make the min distance 3 times bigger when comrin with Gaia
    
    cat1_xy = np.array([cat1['x'], cat1['y']]).T
    gaia1_xy = np.array([gaia['x_1'], gaia['y_1']]).T
    cat1_ga = compare_lists(cat1_xy, gaia1_xy, max_sep.to(u.arcsec).value*g_fac)
    
    cat2_xy = np.array([cat2['x'], cat2['y']]).T
    gaia2_xy = np.array([gaia['x_2'], gaia['y_2']]).T
    cat2_ga = compare_lists(cat2_xy, gaia2_xy, max_sep.to(u.arcsec).value*g_fac)
    
    fig, axm = plt.subplots(1,1)
    axm.set_title('Loop %s. Gaia matches %s'%(wloop_counter, len(gaia['ra'][cat1_ga['ind_2']])))
    axm.scatter(cat1['ra'], cat1['dec'], color = 'k', alpha = 0.01)
    axm.scatter(gaia['ra'][cat1_ga['ind_2']], gaia['dec'][cat1_ga['ind_2']])
    axm.invert_xaxis()
    axm.axis('equal')
    
    
    
   
    dpm_x = (gaia['pmra'][cat1_ga['ind_2']] - cat1['pm_x'][cat1_ga['ind_1']]) 
    dpm_y = (gaia['pmdec'][cat1_ga['ind_2']] - cat1['pm_y'][cat1_ga['ind_1']])
    m_pm, lims = sig_cl(dpm_x, dpm_y, sig_pm)
    
    bad_pm = np.logical_not(m_pm)
    dpm_xm_bad = dpm_x[bad_pm]
    dpm_ym_bad = dpm_y[bad_pm]
    bad_gaia_pm = gaia[cat1_ga['ind_2']][bad_pm]
    
    bad_loop = len(bad)
    
    if len(bad_gaia_pm) > 0:
        if gaia_clipping == 'one_one':
            max_bad_pm = np.argmax(np.sqrt(dpm_xm_bad**2 + dpm_ym_bad**2))
            bad_gaia_pm = bad_gaia_pm[max_bad_pm]
            axm.annotate(bad_gaia_pm['id'], (bad_gaia_pm['ra'], bad_gaia_pm['dec']), xytext=(1, 1), textcoords='offset points', fontsize=15, color='black')
        
        axm.scatter(bad_gaia_pm['ra'], bad_gaia_pm['dec'], facecolor = 'none', s = 200, color = 'red', label = f'cat1 pm {sig_pm}$\sigma$')

        bad.append(bad_gaia_pm['id'])
        
    dx = (gaia['x_1'][cat1_ga['ind_2']] - cat1['x'][cat1_ga['ind_1']])*1e3
    dy = (gaia['y_1'][cat1_ga['ind_2']] - cat1['y'][cat1_ga['ind_1']])*1e3    
    m_xy, limx =  sig_cl(dx, dy, sig_pm)
    bad_xy = np.logical_not(m_xy)
    dx_bad = dx[bad_xy]
    dy_bad = dy[bad_xy]
    bad_gaia_xy = gaia[cat1_ga['ind_2']][bad_xy]
    
    if len(bad_gaia_xy) > 0:
        if gaia_clipping == 'one_one':
            max_bad_xy = np.argmax(np.sqrt(dx_bad**2 + dy_bad**2))
            bad_gaia_xy = bad_gaia_xy[max_bad_xy]
            axm.annotate(bad_gaia_xy['id'], (bad_gaia_xy['ra'], bad_gaia_xy['dec']), xytext=(1, 1), textcoords='offset points', fontsize=15, color='black')
        axm.scatter(bad_gaia_xy['ra'], bad_gaia_xy['dec'], marker = 'x', s = 200, color = 'red', label = f'cat1 xy {sig_pm}$\sigma$')

        bad.append(bad_gaia_xy['id'])
       
    
    dx2 = (gaia['x_2'][cat2_ga['ind_2']] - cat2['x'][cat2_ga['ind_1']])*1e3
    dy2 = (gaia['y_2'][cat2_ga['ind_2']] - cat2['y'][cat2_ga['ind_1']])*1e3
    m_xy2, limx2 =  sig_cl(dx2, dy2, sig_pm)
    bad_xy2 = np.logical_not(m_xy2)
    dx2_bad = dx2[bad_xy2]
    dy2_bad = dy2[bad_xy2]
    bad_gaia_xy2 = gaia[cat2_ga['ind_2']][bad_xy2]
    
    if len(bad_gaia_xy2) > 0:
        if gaia_clipping == 'one_one':
            max_bad_xy2 = np.argmax(np.sqrt(dx2_bad**2 + dy2_bad**2))
            bad_gaia_xy2 = bad_gaia_xy2[max_bad_xy2]
            axm.annotate(bad_gaia_xy2['id'], (bad_gaia_xy2['ra'], bad_gaia_xy2['dec']), xytext=(1, 1), textcoords='offset points', fontsize=15, color='black')
        axm.scatter(bad_gaia_xy2['ra'], bad_gaia_xy2['dec'], marker = '+', s = 200, color = 'red', label = f'cat2 xy {sig_pm}$\sigma$')

        bad.append(bad_gaia_xy2['id'])

    
    axm.legend()
    
    
    
    # if len(dpm_xm_bad) > 0:
    #     max_bad = np.argmax(np.sqrt(dpm_xm_bad**2 + dpm_ym_bad**2))
    #     ax.scatter(bad_gaia_pm['l'], bad_gaia_pm['b'], facecolor = 'none', s = 200, color = 'red')
    #     ax.scatter(bad_gaia_xy['l'], bad_gaia_xy['b'], marker = 'x', s = 200, color = 'red')
    #     ax.scatter(bad_gaia_xy2['l'], bad_gaia_xy2['b'], marker = '+', s = 200, color = 'red')
    #     del_1 = np.isin(gaia['id'], bad_gaia_pm[max_bad]['id'])#!!!
    #     gaia = gaia[np.logical_not(del_1)]#!!!
     
    
    fig, (ax1,ax2, ax3) = plt.subplots(1,3, figsize =(20,7))
    fig.suptitle(f'Aling degree with Gaia = {max_deg-1}, Max sep = {max_sep}, Stars = {len(dpm_x)}')
    ax1.set_title('HST-Gaia pm residuals')
    ax1.scatter(dpm_x,dpm_y, color = 'k', alpha = 0.3)
    # ax.scatter(dpm_xm,dpm_ym,edgecolor = 'k', label = '$\overline{\Delta \mu_{x}}$ = %.2f, $\sigma$ = %.2f\n''$\overline{\Delta \mu_{y}}$ = %.2f, $\sigma$ = %.2f'%(np.mean(dpm_xm),np.std(dpm_xm),np.mean(dpm_ym),np.std(dpm_ym)))
    ax1.axvline(lims[0], color = 'r', ls = 'dashed', label = f'{sig_pm}$\sigma$')
    ax1.axvline(lims[1], color = 'r', ls = 'dashed')
    ax1.axhline(lims[2], color = 'r', ls = 'dashed')
    ax1.axhline(lims[3], color = 'r', ls = 'dashed')
    
    
    for x, y, label in zip(dpm_x[bad_pm],
                           dpm_y[bad_pm],
                           gaia['id'][cat1_ga['ind_2']][bad_pm]):
        print(x, y, label )
        ax1.annotate(str(label), xy=(x, y), xytext=(1,1), textcoords='offset points',
                    fontsize=15, color='black')
    
    # ax.annotate(36, (5,5))
    ax1.legend(fontsize =15)
    ax1.axis('equal')
    ax1.set_xlabel('$\Delta \mu_{x}$ [mas/yr]')
    ax1.set_ylabel('$\Delta \mu_{y}$ [mas/yr]')
    
    for x, y, label in zip(dx[bad_xy],
                           dy[bad_xy],
                           gaia['id'][cat1_ga['ind_2']][bad_xy]):
        print(x, y, label )
        ax2.annotate(str(label), xy=(x, y), xytext=(1, 1), textcoords='offset points',
                    fontsize=15, color='black')
    
    
    ax2.set_title('cat1-Gaia pos. residuals')
    ax2.scatter(dx,dy, color = 'k', alpha = 0.3)
    # ax2.scatter(dx_m,dy_m, edgecolor = 'k',label = '$\overline{\Delta x}$ = %.2f, $\sigma$ = %.2f\n''$\overline{\Delta y}$ = %.2f, $\sigma$ = %.2f'%(np.mean(dx_m),np.std(dx_m),np.mean(dy_m),np.std(dy_m)))
    ax2.axvline(limx[0], color = 'r', ls = 'dashed', label = f'{sig_pm}$\sigma$')
    ax2.axvline(limx[1], color = 'r', ls = 'dashed')
    ax2.axhline(limx[2], color = 'r', ls = 'dashed')
    ax2.axhline(limx[3], color = 'r', ls = 'dashed')
    ax2.legend()
    ax2.axis('equal')
    ax2.set_xlabel('$\Delta$x [mas]')
    ax2.set_ylabel('$\Delta$x [mas]')
    
    for x, y, label in zip(dx2[bad_xy2],
                           dy2[bad_xy2],
                           gaia['id'][cat2_ga['ind_2']][bad_xy2]):
        print(x, y, label )
        ax3.annotate(str(label), xy=(x, y), xytext=(1, 1), textcoords='offset points',
                    fontsize=15, color='black')
    
    ax3.set_title('cat2-Gaia pos. residuals')
    ax3.scatter(dx2,dy2, color = 'k', alpha = 0.3)
    # ax3.scatter(dx_m2,dy_m2,edgecolor = 'k', label = '$\overline{\Delta x}$ = %.2f, $\sigma$ = %.2f\n''$\overline{\Delta y}$ = %.2f, $\sigma$ = %.2f'%(np.mean(dx_m2),np.std(dx_m2),np.mean(dy_m2),np.std(dy_m2)))
    ax3.axvline(limx2[0], color = 'r', ls = 'dashed', label = f'{sig_pm}$\sigma$')
    ax3.axvline(limx2[1], color = 'r', ls = 'dashed')
    ax3.axhline(limx2[2], color = 'r', ls = 'dashed')
    ax3.axhline(limx2[3], color = 'r', ls = 'dashed')
    ax3.legend()
    ax3.axis('equal')
    ax3.set_xlabel('$\Delta$x [mas]')
    ax3.set_ylabel('$\Delta$x [mas]')
    
    fig.tight_layout()
    plt.show()  
    
        
 
   
    # %%
    rcParams.update({
    "figure.figsize": (10, 5),
    "font.size": 18,
    "axes.labelsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16
})
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    
    ax2.set_title(f'# Gaia = {len(dpm_x)} ')
    ax.set_title(f'{zone}, f{band}')
    ax.hist(dpm_x, histtype='step', bins='auto', lw=2,
            label='$\overline{\Delta \mu_{RA}}$ = %.2f'
                  '\n$\sigma_{\mu Ra}$ = %.2f' % 
                  (np.mean(dpm_x), np.std(dpm_x)))
    
    ax2.hist(dpm_y, histtype='step', bins='auto', lw=2,
             label='$\overline{\Delta \mu_{Dec}}$ = %.2f'
                   '\n$\sigma_{\mu Dec}$ = %.2f' % 
                   (np.mean(dpm_y), np.std(dpm_y)))
    
    ax.set_xlabel(r'$\Delta \mu_{RA}$ [mas/yr]')
    ax2.set_xlabel(r'$\Delta \mu_{Dec}$ [mas/yr]')
    ax.set_ylabel('# stars')
    ax.set_xlim(-3,3)
    ax2.set_xlim(-3,3)
    ax.legend(loc=1)
    ax2.legend(loc=1)
    
    fig.tight_layout()

    meta = {'Script': '/Users/amartinez/Desktop/PhD/HAWK/GNS_pm_scripts/GNS_pm_absolute_SUPER/gns_gaia_alignment.py'}
    # plt.savefig(f'/Users/amartinez/Desktop/PhD/My_papers/GNS_pm_catalog/images/ABS_F1_{field_one}_gaia_resi_pm.png', dpi = 150, transparent = True, metadata = meta)


        
    # %%
    # bad = []
    
    
     
    print(30*'☠️')
    print(bad_loop, len(bad))
    bad = np.hstack(bad)
    print('bad = ',[x.tolist() for x in np.unique(bad)])
    # print('bad = ',[x.tolist() for x in np.unique(bad)])
    print(30*'☠️')
        
    bad = list(np.unique(bad))
    
    if len(bad) == bad_loop:
        
        lopping = 0
        print('no  more 3sigmas')
        # gaia[cat1_ga['ind_2']].write(pruebas1  + f'gaia_refstars_F{field_one}_F{field_two}.txt', format = 'ascii', overwrite = True)
        # sys.exit('no  more 3sigmas')

    
    wloop_counter += 1
    
# %%


region_vectors(
    table=cat1[cat1_ga['ind_1']],
    ra_col='ra',
    dec_col='dec',
    pmra_col='pm_x',
    pmdec_col='pm_y',
    name=f'{zone}_{band}_f{band}_vec',
    save_in = pruebas,
    color='cyan',
    wcs='fk5',
    scale=1
)
# 
region_vectors(
    table=gaia[cat1_ga['ind_2']],
    ra_col='ra',
    dec_col='dec',
    pmra_col='pmra',
    pmdec_col='pmdec',
    name=f'gaia_{zone}_{band}_f{band}_vec',
    save_in = pruebas,
    color='red',
    wcs='fk5',
    scale=1,
    width = 5
)




# values = [len(gns1_ga),
#       np.mean(dpm_x), np.std(dpm_x),
#       np.mean(dpm_y), np.std(dpm_y)]

# filepath = pruebas1 + 'alignment_stimation.txt'

# # Check if file already exists
# file_exists = os.path.isfile(filepath)

# # Open in append mode
# with open(filepath, 'a') as f:
#     # If file didn't exist, write header
#     if not file_exists:
#         header = "num_gaia mean_dpm_x std_dpm_x mean_dpm_y std_dpm_y\n"
#         f.write(header)

#     # Write the new line
#     line = " ".join(map(str, values)) + "\n"
#     f.write(line)
# # %%


# # %%
# if look_for_cluster == 'yes':
    
   
#     # modes = ['pm_xy_color']
#     modes = ['pm_xy']
#     knn = 20
#     gen_sim = 'kernnel'
#     sim_lim ='minimun'
#     # sim_lim ='mean'

# # clus_dic = gns_cluster_finder.finder(gns1['pm_xp'], gns1['pm_yp'],
# #                              # gns1['xp'], gns1['yp'], 
# #                              gns1['xp'], gns1['yp'], 
# #                              gns1['l'].value, gns1['b'].value,
# #                              modes[0],
# #                              gns1['H'],gns1['H'],
# #                              knn,gen_sim,sim_lim, save_reg = pruebas1)

#     clus_dic = gns_cluster_finder.finder(gns2['pm_xp'], gns2['pm_yp'],
#                                  # gns1['xp'], gns1['yp'], 
#                                  gns2['xp'], gns2['yp'], 
#                                  gns2['l'].value, gns2['b'].value,
#                                  modes[0],
#                                  gns2['H'],gns2['H'],
#                                  knn,gen_sim,sim_lim, save_reg = pruebas1)






