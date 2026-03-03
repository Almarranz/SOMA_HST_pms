
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
from alignator import alignator
from alignator_relative import alg_rel
import skimage as ski
from astropy.table import Table
from compare_lists import compare_lists
from astropy.stats import sigma_clip
from alignator_looping import alg_loop
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
# %%
# field_one = 10
# chip_one = 0
# field_two = 4
# chip_two = 0

field_one = 'B6'
chip_one = 0
field_two = 20
chip_two = 0

# # field_one = 'D12'
# field_one = 'D13'
# chip_one = 0
# field_two = 1
# chip_two = 0

# field_one = 16
# chip_one = 0
# field_two = 7
# chip_two = 0


# field_one = 'D19'
# chip_one = 0
# field_two = 16
# chip_two = 0

# field_one = 19
# chip_one = 0
# field_two = 9
# chip_two = 0

# field_one = 24
# chip_one = 0
# field_two = 9
# chip_two = 0

# field_one = 1
# chip_one = 0
# field_two = 9
# chip_two = 0
# %%



if field_one == 7 or field_one == 12 or field_one == 10 or field_one == 16 or field_one == 2 or field_one == 24 :
    t1_gns = Time(['2015-06-07T00:00:00'],scale='utc')

elif field_one == 19:
    t1_gns = Time(['2015-06-26T00:00:00'],scale='utc')
    
elif field_one == 'B6' or field_one == 1 :
    t1_gns = Time(['2016-06-13T00:00:00'],scale='utc')
elif field_one ==  'B1':
    t1_gns = Time(['2016-05-20T00:00:00'],scale='utc')
elif field_one ==  'D12':
    t1_gns = Time(['2017-06-03T00:00:00'],scale='utc')
elif field_one ==  'D13':
    t1_gns = Time(['2017-06-24T00:00:00'],scale='utc')
elif field_one ==  'D19':
    t1_gns = Time(['2018-05-21T00:00:00'],scale='utc')
else:
    print(f'NO time detected for this field_one = {field_one}')
    sys.exit()

if field_two == 7 or field_two == 5:
    t2_gns = Time(['2022-05-27T00:00:00'],scale='utc')
elif field_two == 4:
    t2_gns = Time(['2022-04-05T00:00:00'],scale='utc')
elif field_two == 20:
    t2_gns = Time(['2022-07-25T00:00:00'],scale='utc')
elif field_two == 1 or field_two == 9 :
    t2_gns = Time(['2021-09-17T00:00:00'],scale='utc')
elif field_two == 16:
    t2_gns = Time(['2022-08-14T00:00:00'],scale='utc')
else:
    print(f'NO time detected for this field_two = {field_two}')
    sys.exit()

dt_gns = t2_gns - t1_gns

# %%
# ===============================Constants=====================================

# =============================================================================
# Quality cuts
# =============================================================================
max_sig = 0.05


# =============================================================================
# Alignment params
# =============================================================================
rebfacI = 2
rebfacII = 2
# gaia_clipping = 'one_one'# Clipp the Gaia outlayer one by one
gaia_clipping = 'all'# Clipp the Gaia outlayer all at once

# %%
align_loop = 'yes'
# align_loop = 'no'
max_loop = 5
align = 'Polywarp'
# align = '2DPoly'
# f_mode = 'WnC'
sep_both = [50,50]
max_sep_ls = [sep_both[0]*u.mas,sep_both[1]*u.mas]#!!!
max_deg = 3# If this is <2 it does not enter the alignment loop. 

trans_ls = ['polynomial','affine','similarity','Weight']
transf = trans_ls[0]
# pre_transf = trans_ls[2]# Pre transformation
order_trans = 1
clip_in_alig = 'yes' # Clipps the 3sigmas in position during the alignment
# clip_in_alig = None

m_lim = [12 ,18]

# =============================================================================
# Proper motions param
# ===========================================================================
max_dis_pm = 0.150#in arcsec
sig_H = 3# discrd pm for stars with delta H over sig_H
e_pm_gns = 1# im mas/yr

# =============================================================================
# Gaia Settings 
# ============================================================================+
e_pm_gaia = 0.5#!!!
mag_min_gaia = 18#!!!
e_pos_gaia = 0.5
# =============================================================================
# Clustering 
# =============================================================================
# look_for_cluster = 'yes'
look_for_cluster = 'no'

# =============================================================================

GNS_1='/Users/amartinez/Desktop/PhD/HAWK/GNS_1/lists/%s/chip%s/'%(field_one, chip_one)
GNS_2='/Users/amartinez/Desktop/PhD/HAWK/GNS_2/lists/%s/chip%s/'%(field_two, chip_two)

pruebas1 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_1absolute_SUPER/pruebas/'
pruebas2 = '/Users/amartinez/Desktop/PhD/HAWK/GNS_2absolute_SUPER/pruebas/'

bad = []
lopping = 1
# for loop in range(1):
wloop_counter = 0
while lopping > 0:
    
    # gns1 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/pruebas/F{field_one}/{field_one}_H_chips_opti.ecsv',  format = 'ascii.ecsv')
    # gns2 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/pruebas/F{field_two}/{field_two}_H_chips_opti.ecsv', format = 'ascii.ecsv')
    
    
    
    # gns1 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS1/F{field_one}/{field_one}_H_chips_opti.ecsv',  format = 'ascii.ecsv')
    gns2 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS2/F{field_two}/{field_two}_H_chips_opti.ecsv', format = 'ascii.ecsv')
    
    # gns1 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS1/F{field_one}/{field_one}_H_chips_opti_rebfac{rebfacI}.ecsv',  format = 'ascii.ecsv')
    # gns2 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS2/F{field_two}/{field_two}_H_chips_opti_rebfac{rebfacII}.ecsv', format = 'ascii.ecsv')
    
    # gns1 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS1/F{field_one}/{field_one}_H_chips_opti_rebfac{rebfacI}_VVV.ecsv',  format = 'ascii.ecsv')
    # gns2 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS2/F{field_two}/{field_two}_H_chips_opti_rebfac{rebfacII}_VVV.ecsv', format = 'ascii.ecsv')
    
    # gns1 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS1/F{field_one}/{field_one}_H_chips_opti_noDup_rebfac{rebfacI}.ecsv',  format = 'ascii.ecsv')
    # gns2 = Table.read(f'/Users/amartinez/Desktop/Projects/GNS_gd/superlists/GNS2/F{field_two}/{field_two}_H_chips_opti_noDup_rebfac{rebfacII}.ecsv', format = 'ascii.ecsv')
 
    gns1 = Table.read('/Users/amartinez/Desktop/Projects/GNS_gd/pruebas/FB1/B1_and_B6_comb.ecsv',  format = 'ascii.ecsv')

    
    gns1['l'] = Longitude(gns1['l']).wrap_at('180d')
    gns2['l'] = Longitude(gns2['l']).wrap_at('180d')
    
    
    buenos1 = (gns1['l']>min(gns2['l'])) & (gns1['l']<max(gns2['l'])) & (gns1['b']>min(gns2['b'])) & (gns1['b']<max(gns2['b']))
    gns1 = gns1[buenos1]
    
    buenos2 = (gns2['l']>min(gns1['l'])) & (gns2['l']<max(gns1['l'])) & (gns2['b']>min(gns1['b'])) & (gns2['b']<max(gns1['b']))
    gns2 = gns2[buenos2]
    
    
    
    m_mask1 = (gns1['H'] > m_lim[0]) & (gns1['H'] < m_lim[1])
    gns1 = gns1[m_mask1]
    
    unc_cut1 = (gns1['sl']< max_sig) & (gns1['sb'] < max_sig)
    gns1 = gns1[unc_cut1]
    
    m_mask2 = (gns2['H'] > m_lim[0]) & (gns2['H'] < m_lim[1])
    gns2 = gns2[m_mask2]
    
    
    unc_cut2 = (gns2['sl']<max_sig) & (gns2['sb']<max_sig)
    gns2 = gns2[unc_cut2]
    
    
    
    center = SkyCoord(l = np.mean(gns1['l']), b = np.mean(gns1['b']), unit = 'degree', frame = 'galactic')
    # center_1 = SkyCoord(l = np.mean(gns1['l']), b = np.mean(gns1['b']), unit = 'degree', frame = 'galactic')
    # center_2 = SkyCoord(l = np.mean(gns2['l']), b = np.mean(gns2['b']), unit = 'degree', frame = 'galactic')
    
    gns1_lb = SkyCoord(l = gns1['l'], b = gns1['b'], unit ='deg', frame = 'galactic')
    gns2_lb = SkyCoord(l = gns2['l'], b = gns2['b'], unit ='deg', frame = 'galactic')
    
    xg_1, yg_1 = center.spherical_offsets_to(gns1_lb)
    xg_2, yg_2 = center.spherical_offsets_to(gns2_lb)
    
    
    gns1['xp'] = xg_1.to(u.arcsec)
    gns1['yp'] = yg_1.to(u.arcsec)
    gns2['xp'] = xg_2.to(u.arcsec)
    gns2['yp'] = yg_2.to(u.arcsec)
    
    # radius = abs(np.min(gns1['l'])-np.max(gns1['l']))*1*u.degree
    radius = 300*u.arcsec
    # %
    
    gaia = Table.read('/Users/amartinez/Desktop/PhD/Catalogs/Gaia/gaia_over_GNS.txt', format = 'ascii.ecsv')
# =============================================================================
#     try:
#         
#         gaia = Table.read(pruebas1  + 'Noooo_gaia_f1%s_f2%s_r%.0f.ecsv'%(field_one,field_two,radius.to(u.arcsec).value))
#         # gaia = Table.read(pruebas1  + 'gaia_f1%s_f2%s_r400.ecsv'%(field_one,field_two))
#         
#         # gaia = Table.read('/Users/amartinez/Desktop/PhD/HAWK/GNS_1relative_SUPER/pruebas/gaia_f1D19_f216_r269.ecsv')
#         print('Gaia from table')
#     except:
#         print('Gaia from web')
#         # center = SkyCoord(l = np.mean(gns1['l']), b = np.mean(gns1['b']), unit = 'degree', frame = 'galactic').icrs
#     
#         Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select early Data Release 3
#         Gaia.ROW_LIMIT = -1  # it not especifty, Default rows are limited to 50. 
#         # j = Gaia.cone_search_async(center, radius = abs(radius))
#         j = Gaia.cone_search_async(center, radius = radius)
#         gaia = j.get_results()
#         os.makedirs(pruebas1, exist_ok=True)
#         gaia.write(pruebas1  + 'gaia_f1%s_f2%s_r%.0f.ecsv'%(field_one,field_two,radius.to(u.arcsec).value), overwrite = True)
# =============================================================================
    
    gaia['id'] = np.arange(len(gaia))
    gaia['l'] = Longitude(gaia['l']).wrap_at('180d')

   
    

     
    if len(bad)>0:#!!!
        del_1 = np.isin(gaia['id'], bad)#!!!
        gaia = gaia[np.logical_not(del_1)]#!!!
    
    
    
    
    fig, ax2 = plt.subplots(1,1)
    ax2.scatter(gaia['phot_g_mean_mag'],gaia['pmra_error'], s= 2, label = 'Gaia $\delta \mu_{ra}$')
    ax2.scatter(gaia['phot_g_mean_mag'],gaia['pmdec_error'], s= 2, label = 'Gaia $\delta \mu_{dec}$')
    ax2.axvline(mag_min_gaia, color = 'r', ls = 'dashed', label = 'pm cuts')
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
        phot_g_mean_mag_min= mag_min_gaia,
        phot_g_mean_mag_max=None,
        pm_min=0,
        pmra_error_max=e_pm_gaia,
        pmdec_error_max=e_pm_gaia,
        ra_error_max=e_pos_gaia,
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
  
    gaia_lb = SkyCoord(ra = gaia['ra'], dec = gaia['dec'],
                       pm_ra_cosdec = gaia['pmra'],
                       pm_dec = gaia['pmdec'], 
                       frame = 'icrs', obstime = 'J2016').galactic
    xp_g, yp_g = center.spherical_offsets_to(gaia_lb.frame)
    gaia['xp'] = xp_g.to(u.arcsec)
    gaia['yp'] = yp_g.to(u.arcsec)
    
    gaia['pm_l'] = gaia_lb.pm_l_cosb
    gaia['pm_b'] = gaia_lb.pm_b
    
    tg = Time(['2016-01-01T00:00:00'],scale='utc')
    
    
    
    # %
    
    def sig_cl(x, y,s):
        mx, lx, hx = sigma_clip(x , sigma = s, masked = True, return_bounds= True)
        my, ly, hy = sigma_clip(y , sigma = s, masked = True, return_bounds= True)
        m_xy = np.logical_and(np.logical_not(mx.mask),np.logical_not(my.mask))
        
        return m_xy, [lx,hx,ly,hy]
    
    # Define catalogs and times as a list of tuples
    catalogs = [
        {'name': 'GNS1', 'gns': gns1, 'time': t1_gns, 'tag': '1'},
        {'name': 'GNS2', 'gns': gns2, 'time': t2_gns, 'tag': '2'}
    ]
    
    for c,cat in enumerate(catalogs):
        max_sep = max_sep_ls[c]
        print(f"\n===== Aligning {cat['name']} =====")
    
        dt = cat['time'] - Time('2016-01-01T00:00:00', scale='utc')
    
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
        
        gaia_rdt = gaia_rd.apply_space_motion(new_obstime = cat['time'])
        
        gaia_lbt = gaia_rdt.galactic 
        
        gaia[f'l{cat["tag"]}'] = gaia_lbt.l
        gaia[f'b{cat["tag"]}'] = gaia_lbt.b
        
        xp_g, yp_g = center.spherical_offsets_to(gaia_lbt.frame)
        
        gaia['xp'] = xp_g.to(u.arcsec) 
        gaia['yp'] = yp_g.to(u.arcsec) 
        
        gaia[f'xp_{c+1}'] = xp_g.to(u.arcsec) 
        gaia[f'yp_{c+1}'] = yp_g.to(u.arcsec) 
         
    
        # xp_g, yp_g = center.spherical_offsets_to(gaia_lb.frame)
        # gaia['xp'] = xp_g.to(u.arcsec) + gaia_lb.pm_l_cosb * dt.to(u.yr)
        # gaia['yp'] = yp_g.to(u.arcsec) + gaia_lb.pm_b * dt.to(u.yr)
        
        # gaia[f'xp_{c+1}'] = xp_g.to(u.arcsec) + gaia_lb.pm_l_cosb * dt.to(u.yr)
        # gaia[f'yp_{c+1}'] = yp_g.to(u.arcsec) + gaia_lb.pm_b * dt.to(u.yr)
    
        # ga_gtc = center.spherical_offsets_by(gaia['xp'], gaia['yp'])
        # gaia[f'l{cat["tag"]}'] = ga_gtc.l
        # gaia[f'b{cat["tag"]}'] = ga_gtc.b
    
        gns_cat = cat['gns']
        gaia_c = SkyCoord(l=gaia[f'l{cat["tag"]}'], b=gaia[f'b{cat["tag"]}'], frame='galactic')
        gns_c = SkyCoord(l=gns_cat['l'], b=gns_cat['b'], frame='galactic')
        
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
        gns_m = gns_cat[idx2_clean]
        
       
        
        
        print(40*'+')
        unicos = unique(gns_m, keep = 'first')
        print(len(gns_m),len(unicos))
        print(40*'+')
        
# =============================================================================
#         idx, d2d, _ = gaia_c.match_to_catalog_sky(gns_c, nthneighbor=1)
#         match_mask = d2d < max_sep
#         gaia_m = gaia[match_mask]
#         gns_m = gns_cat[idx[match_mask]]
#         print(40*'+')
#         unicos = unique(gns_m, keep = 'first')
#         print(len(gns_m),len(unicos))
#         print(40*'+')
#     
# =============================================================================
       
        
        fig, (ax, ax1, ax2) = plt.subplots(1,3, figsize = (15,5))
        fig.suptitle(f'GNS{c+1}')
        # ax.scatter(gns_w[::100], gns_cat['b'][::100], alpha =0.1, color = 'k')
        # ax.scatter(ga_w, gaia['b'], label = 'Gaia', s= 10)
        ax.scatter(gns_cat['l'][::10], gns_cat['b'][::10], alpha =0.1, color = 'k')
        ax.set_title(f'Matches = {len(gns_m)}\nMin dist = {max_sep} ')
        ax.scatter(gns_m['l'], gns_m['b'], label = f'GNS{c} Match')
        ax.scatter(gaia_m['l'], gaia_m['b'],s =10, label = 'Gaia Match')
        ax.set_xlabel('l')
        ax.set_ylabel('b')
        ax.legend()
        
# =============================================================================
#         # Apply a similiraty trasnformation fisrt
#         Does not make much of a difference
# =============================================================================
        
# =============================================================================
#         if pre_transf == 'similarity':
#             psim = ski.transform.estimate_transform(
#             pre_transf, np.array([gns_m['l'], gns_m['b']]).T,
#                 np.array([gaia_m['l'], gaia_m['b']]).T)
#         
#         # if pre_transf == 'affine':
#         #     psim = ski.transform.estimate_transform(
#         #     pre_transf, np.array([gns_m['l'], gns_m['b']]).T,
#         #         np.array([gaia_m['l'], gaia_m['b']]).T)
#         
#         else:
#             psim = ski.transform.estimate_transform(
#             pre_transf, np.array([gns_m['l'], gns_m['b']]).T,
#                 np.array([gaia_m['l'], gaia_m['b']]).T, order = 1)
#             
#         gns_sim = psim(np.array([gns_cat['l'], gns_cat['b']]).T)
#         gns_cat['l'] = gns_sim[:, 0]*u.deg
#         gns_cat['b'] = gns_sim[:, 1]*u.deg
# =============================================================================
# =============================================================================
#         # Apply a similiraty trasnformation fisrt
#         Does not make much of a difference
# =============================================================================
        ns_cat = cat['gns']
        gaia_c = SkyCoord(l=gaia[f'l{cat["tag"]}'], b=gaia[f'b{cat["tag"]}'], frame='galactic')
        gns_c = SkyCoord(l=gns_cat['l'], b=gns_cat['b'], frame='galactic')
        
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
        gns_m = gns_cat[idx2_clean]
        
        # Number of rows to select
       


        print(40*'+')
        unicos = unique(gns_m, keep = 'first')
        print(len(gns_m),len(unicos))
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
            dl_pre = gns_m['l'].to(u.mas) - gaia_m['l'].to(u.mas)
            db_pre = gns_m['b'].to(u.mas) - gaia_m['b'].to(u.mas)
            m_pre, pre_lims = sig_cl(dl_pre, db_pre, 3)
            dl_prem = dl_pre[m_pre]
            db_prem = db_pre[m_pre]
            
            fig_pre, (ax_pre, ax1_pre) = plt.subplots(1,2,figsize = (11,5.5))
            #
            ax_pre.set_title(f'Gaia vs GNS{c+1} (Stars = {len(gaia_m)})')
            ax1_pre.set_title(f'Residuas before any transformation')
            # ax1.set_title(f'Matching stars  = {len(d_xm)}')
            ax_pre.set_ylabel('# stars')
            ax_pre.hist(dl_pre, bins = 10,  color = 'grey', alpha = 0.5)
            ax1_pre.hist(db_pre, bins = 'auto',  color = 'grey', alpha = 0.5)
            ax_pre.hist(dl_prem, bins = 10, histtype = 'step',lw = 2,label = '$\overline{\Delta l}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(dl_prem.value),np.std(dl_prem.value)))
            ax1_pre.hist(db_prem,bins = 'auto',histtype = 'step', lw = 2, label = '$\overline{\Delta b}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(db_prem.value),np.std(db_prem.value)))
            ax_pre.legend(loc = 1)
            ax1_pre.legend()
            ax_pre.set_xlabel('$\Delta$l [mas]')
            ax1_pre.set_xlabel('$\Delta$b [mas]')
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
                transf, np.array([gns_m['xp'], gns_m['yp']]).T,
                np.array([gaia_m['xp'], gaia_m['yp']]).T, order=order_trans
            )
        elif transf == 'Weight':
            model_x = Polynomial2D(degree=order_trans)
            model_y = Polynomial2D(degree=order_trans)
            
            # Linear least-squares fitter
            fitter = LinearLSQFitter()
            
            fit_xw = fitter(model_x, gns_m['xp'], gns_m['yp'],  gaia_m['xp'], weights= 1/np.sqrt( gns_m['sl']**2 +  gns_m['sb']**2))  # Fit x-coordinates
            fit_yw = fitter(model_y, gns_m['xp'], gns_m['yp'],  gaia_m['yp'],weights= 1/np.sqrt( gns_m['sl']**2 +  gns_m['sb']**2)) 
            
            
        else:
            p = ski.transform.estimate_transform(
                transf, np.array([gns_m['xp'], gns_m['yp']]).T,
                np.array([gaia_m['xp'], gaia_m['yp']]).T
            )
       
    
        
        if transf == 'Weight':
            
            gns_cat['xp']  = fit_xw(gns_cat['xp'], gns_cat['yp'])
            gns_cat['yp']  = fit_yw(gns_cat['xp'], gns_cat['yp'])
        
        else:
            gns_trans = p(np.array([gns_cat['xp'], gns_cat['yp']]).T)
            gns_cat['xp'] = gns_trans[:, 0]
            gns_cat['yp'] = gns_trans[:, 1]
            print(p.params)
        
        gns_xy = np.array([gns_cat['xp'],gns_cat['yp']]).T
        gaia_xy = np.array([gaia['xp'],gaia['yp']]).T
        xy_mat = compare_lists(gns_xy, gaia_xy, max_sep.to(u.arcsec).value)
        
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
        ax1.scatter(gns_cat['xp'][::100], gns_cat['yp'][::100], alpha =0.1, color = 'k')
        ax1.scatter(gaia['xp'], gaia['yp'],s =10, label = 'Gaia')
        ax1.scatter(gns_cat['xp'][xy_mat['ind_1']], gns_cat['yp'][xy_mat['ind_1']], label = 'GNS1 match')
        ax1.scatter(gaia['xp'][xy_mat['ind_2']], gaia['yp'][xy_mat['ind_2']],s =10, label = 'Gaia match')
        ax1.set_xlabel('xp[arcsec]')
        ax1.legend()
        
        if align_loop  == 'yes':
            # Optional final refinement
            # gns_cat = alg_rel(gns_cat, gaia, 'xp', 'yp', align, max_deg, max_sep.to(u.arcsec).value)
            gns_cat = alg_loop(gns_cat, gaia, 'xp', 'yp', align, max_deg, max_sep.to(u.arcsec).value, max_loop)
            
        elif align_loop == 'no':
            p2 = ski.transform.estimate_transform('polynomial', 
                np.array([gns_cat['xp'][xy_mat['ind_1']], gns_cat['yp'][xy_mat['ind_1']]]).T,
                np.array([gaia['xp'][xy_mat['ind_2']], gaia['yp'][xy_mat['ind_2']]]).T, order=2)
            
            gns_trans = p2(np.array([gns_cat['xp'], gns_cat['yp']]).T)
            
            gns_cat['xp'] = gns_trans[:,0]
            gns_cat['yp'] = gns_trans[:,1]
        
        
        
        gns_xy = np.array([gns_cat['xp'],gns_cat['yp']]).T
        gaia_xy = np.array([gaia['xp'],gaia['yp']]).T
        xy_mat = compare_lists(gns_xy, gaia_xy, max_sep.to(u.arcsec).value)
        
        ax2.set_title(f'Matches = {len(xy_mat)}\nMin dist = {max_sep} ')
        ax2.scatter(gns_cat['xp'][::100], gns_cat['yp'][::100], alpha =0.1, color = 'k')
        ax2.scatter(gaia['xp'], gaia['yp'],s =10, label = 'Gaia')
        ax2.scatter(gns_cat['xp'][xy_mat['ind_1']], gns_cat['yp'][xy_mat['ind_1']], label = 'GNS1 match')
        ax2.scatter(gaia['xp'][xy_mat['ind_2']], gaia['yp'][xy_mat['ind_2']],s =10, label = 'Gaia match')
        ax2.set_xlabel('xp[arcsec]')
        ax2.yaxis.set_label_position("right")
        ax2.set_ylabel('yp [arcsec]')
        ax2.legend()
        fig.tight_layout()
        plt.show()
        # Residuals
        gns_xy = np.array([gns_cat['xp'], gns_cat['yp']]).T
        gaia_xy = np.array([gaia['xp'], gaia['yp']]).T
        xy_al = compare_lists(gns_xy, gaia_xy, max_sep.to(u.arcsec).value)
    
        d_x = (gaia['xp'][xy_al['ind_2']] - gns_cat['xp'][xy_al['ind_1']]).to(u.mas) 
        d_y = (gaia['yp'][xy_al['ind_2']] -gns_cat['yp'][xy_al['ind_1']] ).to(u.mas)
        
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
        
        cat['gns']['xp'] =  gns_cat['xp']
        cat['gns']['yp'] =  gns_cat['yp']
        plt.show()
    
        # stop(724)
    
    
    
        # print('Pasannado')
    # %%
    
    # =============================================================================
    # GNS proper motions
    # =============================================================================
    gns1_al = catalogs[0]['gns']
    gns2_al = catalogs[1]['gns']
    
    
    
    gns1_gxy  = np.array([gns1_al['xp'], gns1_al['yp']]).T 
    gns2_gxy  = np.array([gns2_al['xp'], gns2_al['yp']]).T 
    
    gns_com = compare_lists(gns1_gxy, gns2_gxy,max_dis_pm )
    
    gns1_gxy  = gns1_gxy[gns_com['ind_1']]  
    gns2_gxy  = gns2_gxy[gns_com['ind_2']]  
    
    gns1 = gns1_al[gns_com['ind_1']]
    gns2 = gns2_al[gns_com['ind_2']]
    
    dH = gns1['H'] - gns2['H']
    mask_H, lH, hH = sigma_clip(dH , sigma = sig_H, masked = True, return_bounds= True)
    
    fig, ax = plt.subplots(1,1)
    ax.hist(dH, bins = 'auto', label = '$\overline{\Delta H}$ = %.2f\n$\sigma$ = %.2f'%(np.mean(dH), np.std(dH)))
    ax.set_xlabel('$\Delta$H GNS')
    ax.axvline(lH, color = 'red', ls = 'dashed', label = f'{sig_H}$\sigma$')
    ax.axvline(hH, color = 'red', ls = 'dashed')
    ax.legend()

    
    # sys.exit(517)
    
    gns1 = gns1[np.logical_not(mask_H.mask)]
    gns2 = gns2[np.logical_not(mask_H.mask)]
    gns_com = gns_com[np.logical_not(mask_H.mask)]
    gns1_gxy = gns1_gxy[np.logical_not(mask_H.mask)]
    gns2_gxy = gns2_gxy[np.logical_not(mask_H.mask)]
    
    
    pm_x = (gns_com['l2_x']*u.arcsec - gns_com['l1_x']*u.arcsec).to(u.mas)/dt_gns.to(u.yr)
    pm_y = (gns_com['l2_y']*u.arcsec - gns_com['l1_y']*u.arcsec).to(u.mas)/dt_gns.to(u.yr)
    
    dpm_x = np.sqrt((gns1['sl'].to(u.mas))**2 + (gns1['sl'].to(u.mas))**2)/dt_gns.to(u.year)
    dpm_y = np.sqrt((gns2['sb'].to(u.mas))**2 + (gns1['sb'].to(u.mas))**2)/dt_gns.to(u.year)

    
    gns1['pm_xp'] = pm_x
    gns1['pm_yp'] = pm_y
    gns2['pm_xp'] = pm_x
    gns2['pm_yp'] = pm_y
    
    gns1['dpm_x'] = dpm_x
    gns1['dpm_y'] = dpm_y
    gns2['dpm_x'] = dpm_x
    gns2['dpm_y'] = dpm_y
    
    
    
    gns1.meta['Gaia_sep'] = sep_both*u.mas
    gns1.meta['m_lim_gns'] = m_lim
    gns1.meta['max_dis_pm'] = max_dis_pm
    gns1.meta['e_pm_gaia'] = e_pm_gaia
    gns1.meta['mag_min_gaia'] = mag_min_gaia
    
# =============================================================================
#     Save the catlogs with the pm
# =============================================================================
    gns1.write(pruebas1 + f'gns1_pmSuper_F1_{field_one}_F2_{field_two}.ecvs',format = 'ascii.ecsv', overwrite = True)
    gns2.write(pruebas2 + f'gns2_pmSuper_F1_{field_one}_F2_{field_two}.ecvs',format = 'ascii.ecsv', overwrite = True)
    
    


    
    # gns2.meta['f_mode'] = f_mode


    
    gns1 = filter_gns_data(gns1, max_e_pm = e_pm_gns)
    gns2 = filter_gns_data(gns2, max_e_pm = e_pm_gns)

    # 
    # %%
    fig, (ax,ax2) = plt.subplots(1,2)
    
    
    bins = 30
    ax.hist(pm_x, bins = bins, color = 'grey', alpha = 0.3)
    ax2.hist(pm_y, bins = bins, color = 'grey', alpha = 0.3)
    
    ax.hist(gns1['pm_xp'], bins = bins, histtype = 'step', label = '$\overline{\mu}_{xp}$ = %.2f\n$\sigma$ =%.2f'%(np.mean(gns1['pm_xp']),np.std(gns1['pm_xp'])))
    # ax.axvline(lpm[0] , ls = 'dashed', color = 'r')
    # ax.axvline(lpm[1] , ls = 'dashed', color = 'r')
    
    ax2.hist(gns1['pm_yp'], bins = bins, histtype = 'step', label = '$\overline{\mu}_{yp}$ = %.2f\n$\sigma$ =%.2f'%(np.mean(gns1['pm_yp']),np.std(gns1['pm_yp'])))
    # ax2.axvline(lpm[2] , ls = 'dashed', color = 'r')
    # ax2.axvline(lpm[3] , ls = 'dashed', color = 'r')
    
    ax.legend()
    ax2.legend()
    
    ax.set_xlabel('$\mu_{xp}$ [mas/yr]')
    ax2.set_xlabel('$\mu_{yp}$ [mas/yr]')
    fig.tight_layout()
    plt.show()
    
    g_fac = 1# make the min distance 3 times bigger when comrin with Gaia
    
    gns1_xy = np.array([gns1['xp'], gns1['yp']]).T
    gaia1_xy = np.array([gaia['xp_1'], gaia['yp_1']]).T
    gns1_ga = compare_lists(gns1_xy, gaia1_xy, max_sep_ls[0].to(u.arcsec).value*g_fac)
    
    gns2_xy = np.array([gns2['xp'], gns2['yp']]).T
    gaia2_xy = np.array([gaia['xp_2'], gaia['yp_2']]).T
    gns2_ga = compare_lists(gns2_xy, gaia2_xy, max_sep_ls[1].to(u.arcsec).value*g_fac)
    
    fig, axm = plt.subplots(1,1)
    axm.set_title('Loop %s. Gaia matches %s'%(wloop_counter, len(gaia['l'][gns1_ga['ind_2']])))
    axm.scatter(gns1['l'], gns1['b'], color = 'k', alpha = 0.01)
    axm.scatter(gaia['l'][gns1_ga['ind_2']], gaia['b'][gns1_ga['ind_2']])
    axm.invert_xaxis()
    axm.axis('equal')
    
   
    
   
    dpm_x = (gaia['pm_l'][gns1_ga['ind_2']] - gns1['pm_xp'][gns1_ga['ind_1']]) 
    dpm_y = (gaia['pm_b'][gns1_ga['ind_2']] - gns1['pm_yp'][gns1_ga['ind_1']])
    m_pm, lims = sig_cl(dpm_x, dpm_y, sig_pm)
    
    bad_pm = np.logical_not(m_pm)
    dpm_xm_bad = dpm_x[bad_pm]
    dpm_ym_bad = dpm_y[bad_pm]
    bad_gaia_pm = gaia[gns1_ga['ind_2']][bad_pm]
    
    bad_loop = len(bad)
    
    if len(bad_gaia_pm) > 0:
        if gaia_clipping == 'one_one':
            max_bad_pm = np.argmax(np.sqrt(dpm_xm_bad**2 + dpm_ym_bad**2))
            bad_gaia_pm = bad_gaia_pm[max_bad_pm]
            axm.annotate(bad_gaia_pm['id'], (bad_gaia_pm['l'], bad_gaia_pm['b']), xytext=(1, 1), textcoords='offset points', fontsize=15, color='black')
        
        axm.scatter(bad_gaia_pm['l'], bad_gaia_pm['b'], facecolor = 'none', s = 200, color = 'red', label = f'GNS1 pm {sig_pm}$\sigma$')

        bad.append(bad_gaia_pm['id'])
        
    dx = (gaia['xp_1'][gns1_ga['ind_2']] - gns1['xp'][gns1_ga['ind_1']])*1e3
    dy = (gaia['yp_1'][gns1_ga['ind_2']] - gns1['yp'][gns1_ga['ind_1']])*1e3    
    m_xy, limx =  sig_cl(dx, dy, sig_pm)
    bad_xy = np.logical_not(m_xy)
    dx_bad = dx[bad_xy]
    dy_bad = dy[bad_xy]
    bad_gaia_xy = gaia[gns1_ga['ind_2']][bad_xy]
    
    if len(bad_gaia_xy) > 0:
        if gaia_clipping == 'one_one':
            max_bad_xy = np.argmax(np.sqrt(dx_bad**2 + dy_bad**2))
            bad_gaia_xy = bad_gaia_xy[max_bad_xy]
            axm.annotate(bad_gaia_xy['id'], (bad_gaia_xy['l'], bad_gaia_xy['b']), xytext=(1, 1), textcoords='offset points', fontsize=15, color='black')
        axm.scatter(bad_gaia_xy['l'], bad_gaia_xy['b'], marker = 'x', s = 200, color = 'red', label = f'GNS1 xy {sig_pm}$\sigma$')

        bad.append(bad_gaia_xy['id'])
       
    
    dx2 = (gaia['xp_2'][gns2_ga['ind_2']] - gns2['xp'][gns2_ga['ind_1']])*1e3
    dy2 = (gaia['yp_2'][gns2_ga['ind_2']] - gns2['yp'][gns2_ga['ind_1']])*1e3
    m_xy2, limx2 =  sig_cl(dx2, dy2, sig_pm)
    bad_xy2 = np.logical_not(m_xy2)
    dx2_bad = dx2[bad_xy2]
    dy2_bad = dy2[bad_xy2]
    bad_gaia_xy2 = gaia[gns2_ga['ind_2']][bad_xy2]
    
    if len(bad_gaia_xy2) > 0:
        if gaia_clipping == 'one_one':
            max_bad_xy2 = np.argmax(np.sqrt(dx2_bad**2 + dy2_bad**2))
            bad_gaia_xy2 = bad_gaia_xy2[max_bad_xy2]
            axm.annotate(bad_gaia_xy2['id'], (bad_gaia_xy2['l'], bad_gaia_xy2['b']), xytext=(1, 1), textcoords='offset points', fontsize=15, color='black')
        axm.scatter(bad_gaia_xy2['l'], bad_gaia_xy2['b'], marker = '+', s = 200, color = 'red', label = f'GNS2 xy {sig_pm}$\sigma$')

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
    fig.suptitle(f'Aling degree with Gaia = {max_deg-1}, Max sep = {max_sep_ls}, Stars = {len(dpm_x)}')
    ax1.set_title('GNS-Gaia pm residuals')
    ax1.scatter(dpm_x,dpm_y, color = 'k', alpha = 0.3)
    # ax.scatter(dpm_xm,dpm_ym,edgecolor = 'k', label = '$\overline{\Delta \mu_{x}}$ = %.2f, $\sigma$ = %.2f\n''$\overline{\Delta \mu_{y}}$ = %.2f, $\sigma$ = %.2f'%(np.mean(dpm_xm),np.std(dpm_xm),np.mean(dpm_ym),np.std(dpm_ym)))
    ax1.axvline(lims[0], color = 'r', ls = 'dashed', label = f'{sig_pm}$\sigma$')
    ax1.axvline(lims[1], color = 'r', ls = 'dashed')
    ax1.axhline(lims[2], color = 'r', ls = 'dashed')
    ax1.axhline(lims[3], color = 'r', ls = 'dashed')
    
    
    for x, y, label in zip(dpm_x[bad_pm],
                           dpm_y[bad_pm],
                           gaia['id'][gns1_ga['ind_2']][bad_pm]):
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
                           gaia['id'][gns1_ga['ind_2']][bad_xy]):
        print(x, y, label )
        ax2.annotate(str(label), xy=(x, y), xytext=(1, 1), textcoords='offset points',
                    fontsize=15, color='black')
    
    
    ax2.set_title('GNS1-Gaia pos. residuals')
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
                           gaia['id'][gns2_ga['ind_2']][bad_xy2]):
        print(x, y, label )
        ax3.annotate(str(label), xy=(x, y), xytext=(1, 1), textcoords='offset points',
                    fontsize=15, color='black')
    
    ax3.set_title('GNS2-Gaia pos. residuals')
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
    
    ax.set_title(f'# Gaia = {len(dpm_x)} ')
    ax.hist(dpm_x, histtype='step', bins='auto', lw=2,
            label='$\overline{\Delta \mu_{l}}$ = %.2f'
                  '\n$\sigma_b$ = %.2f' % 
                  (np.mean(dpm_x), np.std(dpm_x)))
    
    ax2.hist(dpm_y, histtype='step', bins='auto', lw=2,
             label='$\overline{\Delta \mu_{b}}$ = %.2f'
                   '\n$\sigma_l$ = %.2f' % 
                   (np.mean(dpm_y), np.std(dpm_y)))
    
    ax.set_xlabel(r'$\Delta \mu_{l}$ [mas/yr]')
    ax2.set_xlabel(r'$\Delta \mu_{b}$ [mas/yr]')
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
        gaia[gns1_ga['ind_2']].write(pruebas1  + f'gaia_refstars_F{field_one}_F{field_two}.txt', format = 'ascii', overwrite = True)
        # sys.exit('no  more 3sigmas')
    
    wloop_counter += 1
    
# stop(973)
values = [len(gns1_ga),
      np.mean(dpm_x), np.std(dpm_x),
      np.mean(dpm_y), np.std(dpm_y)]

filepath = pruebas1 + 'alignment_stimation.txt'

# Check if file already exists
file_exists = os.path.isfile(filepath)

# Open in append mode
with open(filepath, 'a') as f:
    # If file didn't exist, write header
    if not file_exists:
        header = "num_gaia mean_dpm_x std_dpm_x mean_dpm_y std_dpm_y\n"
        f.write(header)

    # Write the new line
    line = " ".join(map(str, values)) + "\n"
    f.write(line)
# %%


# %%
if look_for_cluster == 'yes':
    
   
    # modes = ['pm_xy_color']
    modes = ['pm_xy']
    knn = 20
    gen_sim = 'kernnel'
    sim_lim ='minimun'
    # sim_lim ='mean'

# clus_dic = gns_cluster_finder.finder(gns1['pm_xp'], gns1['pm_yp'],
#                              # gns1['xp'], gns1['yp'], 
#                              gns1['xp'], gns1['yp'], 
#                              gns1['l'].value, gns1['b'].value,
#                              modes[0],
#                              gns1['H'],gns1['H'],
#                              knn,gen_sim,sim_lim, save_reg = pruebas1)

    clus_dic = gns_cluster_finder.finder(gns2['pm_xp'], gns2['pm_yp'],
                                 # gns1['xp'], gns1['yp'], 
                                 gns2['xp'], gns2['yp'], 
                                 gns2['l'].value, gns2['b'].value,
                                 modes[0],
                                 gns2['H'],gns2['H'],
                                 knn,gen_sim,sim_lim, save_reg = pruebas1)
#     elif destination == 1:
#         clus_dic = gns_cluster_finder.finder(gns1_mpm['pm_x'], gns1_mpm['pm_y'],
#                                      gns1_mpm['xp'], gns1_mpm['yp'], 
#                                      gns1_mpm['l'].value, gns1_mpm['b'].value,
#                                      modes[0],
#                                      gns1_mpm['H'],gns1_mpm['H'],
#                                      knn,gen_sim,sim_lim, save_reg = pruebas1)
# # %%


# =============================================================================
# zone =  'F20_f01_f06_H'
# 
# scamp_f = '/Users/amartinez/Desktop/Projects/GNS_gd/scamp/GNS0/%s/'%(zone)
# 
# try:
#     cat = Table.read(scamp_f +f'merged_{zone}_1.ocat', format = 'ascii') 
# except:
#     cat = Table.read(scamp_f +f'merged_{zone}.ocat', format = 'ascii') 
# vel_max = 50
# pm_mask = (abs(cat['PMALPHA_J2000']) <vel_max) &  (abs(cat['PMDELTA_J2000']) <vel_max) & (abs(cat['PMALPHA_J2000']) >1e-5) &  (abs(cat['PMDELTA_J2000']) >1e-5)
# cat = cat[pm_mask]
# 
# cat_c = SkyCoord(ra = cat['ALPHA_J2000'], dec = cat['DELTA_J2000'], frame = 'fk5').galactic
# gns_c = SkyCoord(l = gns1['l'], b = gns1['b'], frame = 'galactic')
# 
# 
# 
# # %
# max_seo = 50*u.mas
# idx, d2d, _ = cat_c.match_to_catalog_sky(gns_c, nthneighbor=1)
# match_mask = d2d < max_sep
# cat_m = cat[match_mask]
# gns_m = gns1[idx[match_mask]]
# 
# # %
# fig, ax = plt.subplots(1,1)
# ax.scatter(cat_c.l,cat_c.b)
# ax.scatter(cat_c.l[match_mask],cat_c.b[match_mask])
# 
# 
# dpmx = (gns_m['pm_xp'] - cat_m['PMALPHA_J2000'])
# dpmy = (gns_m['pm_yp'] - cat_m['PMDELTA_J2000'])
# 
# m_pm, lim = sig_cl(dpmx, dpmy, 3) 
# 
# dpmx_m = dpmx[m_pm]
# dpmy_m = dpmy[m_pm]
# 
# 
# fig, ax = plt.subplots(1,1)
# 
# ax.scatter(dpmx,dpmy, color = 'k', alpha = 0.3)
# ax.axvline(lims[0], color = 'r', ls = 'dashed', label = f'{sig_pm}$\sigma$')
# ax.scatter(dpmx_m,dpmy_m,edgecolor = 'k', label = '$\overline{\Delta \mu_{x}}$ = %.2f, $\sigma$ = %.2f\n''$\overline{\Delta \mu_{y}}$ = %.2f, $\sigma$ = %.2f'%(np.mean(dpmx_m),np.std(dpmx_m),np.mean(dpmy_m),np.std(dpmy_m)))
# ax.axvline(lim[1], color = 'r', ls = 'dashed')
# ax.axhline(lim[2], color = 'r', ls = 'dashed')
# ax.axhline(lim[3], color = 'r', ls = 'dashed')
# ax.legend()
# 
# =============================================================================







