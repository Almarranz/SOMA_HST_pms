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
from matplotlib.colors import PowerNorm
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
# results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
# results = ''
# zone = 'G032.03+00.05'
# zone = 'G028.20-00.05'
# zones = [ 'AFGL5180','G028.20-00.05', 'G032.03+00.05', 'G35.2-0.74N', 'G339.88-01.26','IRAS07299-1651', 'IRAS16562-3959']
zones = ['AFGL5180']
color_l = None
# zone = zones[2]
epoch = 1
for zone in zones:
    # band = '160w'
    band1 = '110w'
    band2 = '160w'
    band3 = '128n'
    band4 = '164n'
    # epoch = 2
    if band1 == '160w':
        b_name1 = 'H'
    elif band1 == '110w':
        b_name1 = 'J'
    elif band1 == '128n':
        b_name1 = 'Paβ'
    elif band1 == '164n':
        b_name1 = 'Fe II'
    if band2 == '160w':
        b_name2 = 'H'
    elif band2 == '110w':
        b_name2 = 'J'
    elif band2 == '128n':
        b_name2 = 'Paβ'
    elif band2 == '164n':
        b_name2 = 'Fe II'
    if band3 == '160w':
        b_name3 = 'H'
    elif band3 == '110w':
        b_name3 = 'J'
    elif band3 == '128n':
        b_name3 = 'Paβ'
    elif band3 == '164n':
        b_name3 = 'Fe II'
    if band4 == '160w':
        b_name4 = 'H'
    elif band4 == '110w':
        b_name4 = 'J'
    elif band4 == '128n':
        b_name4 = 'Paβ'
    elif band4 == '164n':
        b_name4 = 'Fe II'
        
    # =============================================================================
    # Alig Para
    # =============================================================================
    max_sep = 50*u.mas
    
    
    pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'
    tmp2 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band2}/epoch1/tmp/'

    results1 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band1}/epoch{epoch}/'
    results2 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band2}/epoch{epoch}/'
    results3 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band3}/epoch{epoch}/'
    results4 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band4}/epoch{epoch}/'

    # cat.write(results + f'calib_{zone}_EP{epoch}_f{band}.txt', format = 'ascii', overwrite = True )

    
    
    cat1 = Table.read(results1 +  f'calib_{zone}_EP{epoch}_f{band1}.txt', format = 'ascii')
    cat2 = Table.read(results2 +  f'calib_{zone}_EP{epoch}_f{band2}.txt', format = 'ascii')
    cat3 = Table.read(results3 +  f'calib_{zone}_EP{epoch}_f{band3}.txt', format = 'ascii')
    cat4 = Table.read(results4 +  f'calib_{zone}_EP{epoch}_f{band4}.txt', format = 'ascii')
    
    
    
    crd1 = SkyCoord(ra = cat1['ra'], dec = cat1['dec'], unit = 'degree', frame = 'fk5')
    crd2 = SkyCoord(ra = cat2['ra'], dec = cat2['dec'], unit = 'degree', frame = 'fk5')
    crd3 = SkyCoord(ra = cat3['ra'], dec = cat3['dec'], unit = 'degree', frame = 'fk5')
    crd4 = SkyCoord(ra = cat4['ra'], dec = cat4['dec'], unit = 'degree', frame = 'fk5')
    
    
    idx, d2d, _ = crd1.match_to_catalog_sky(crd2, nthneighbor=1)
    match_mask = d2d < max_sep
    cat1_m = cat1[match_mask]
    cat2_m = cat2[idx[match_mask]]
    print(40*'+')
    unicos = unique(cat1_m, keep = 'first')
    print(len(cat1_m),len(unicos))
    print(40*'+')
    
   
    
    fig, ax = plt.subplots(1,1, figsize =(4,4)) 
    ax.set_title(f'{zone} \u2605 = {len(cat1_m)}')
    h = ax.hexbin(cat1_m[b_name1] - cat2_m[b_name2], cat2_m[b_name2], norm = LogNorm(),cmap = 'viridis' )
    cbar = plt.colorbar(h, ax=ax, aspect = 40)
    ax.invert_yaxis()
    ax.set_xlabel(f'{b_name1} - {b_name2}')
    ax.set_ylabel(f'{b_name2}')
   
    
   
    if color_l is not None:
        ax.axvline(color_l, color = 'red', ls = 'dashed',label = f'J-H = {color_l}')
        ax.legend()
        mask_c = (cat1_m[b_name1] - cat2_m[b_name2]) > color_l
        
        cat1_mC = cat1_m[mask_c]
        cat1_NmC = cat1_m[np.logical_not(mask_c)]
        
        region(cat1_mC, 'ra', 'dec',
               name= f'{zone}_JH_gt{color_l}',
               save_in=pruebas,wcs='fk5', color = 'green', marker= 'cross')
        region(cat1_NmC, 'ra', 'dec',
               name= f'{zone}_H-K_st{color_l}',
               save_in=pruebas,wcs='fk5',color= 'cyan',
               marker = 'circle')
    
# %%
    # =============================================================================
    # Select diagram type
    # =============================================================================
    
    for dtype in range(1,8):
        # dtype = 7
        
        # CMDs
        # 1 -> (F110W-F160W) vs F160W
        # 2 -> (F128N-F160W) vs F160W
        # 3 -> (F164N-F160W) vs F160W
        
        # CCDs
        # 4 -> (F110W-F160W) vs (F128N-F160W)
        # 5 -> (F110W-F160W) vs (F164N-F160W)
        # 6 -> (F128N-F160W) vs (F164N-F160W)
        # 7 -> (F110W-F160W) vs (F128N-F164N)
        
        
        # =============================================================================
        # Match catalogs
        # =============================================================================
        
        idx12, d12, _ = crd2.match_to_catalog_sky(crd1)
        m12 = d12 < max_sep
        
        idx32, d32, _ = crd2.match_to_catalog_sky(crd3)
        m32 = d32 < max_sep
        
        idx42, d42, _ = crd2.match_to_catalog_sky(crd4)
        m42 = d42 < max_sep
        
        good = m12 & m32 & m42
        
        cat2m = cat2[good]          # F160W reference
        cat1m = cat1[idx12[good]]   # F110W
        cat3m = cat3[idx32[good]]   # F128N
        cat4m = cat4[idx42[good]]   # F164N
        
        
        # =============================================================================
        # Definitions
        # =============================================================================
        
        colors = {
            'J-H'      : cat1m['J'] - cat2m['H'],
            'Paβ-H'    : cat3m['Paβ'] - cat2m['H'],
            'FeII-H'   : cat4m['Fe II'] - cat2m['H'],
            'Paβ-FeII' : cat3m['Paβ'] - cat4m['Fe II'],
        }
        
        mags = {
            'J'     : cat1m['J'],
            'H'     : cat2m['H'],
            'Paβ'   : cat3m['Paβ'],
            'Fe II' : cat4m['Fe II'],
        }
        
        
        # =============================================================================
        # Choose diagram
        # =============================================================================
        
        if dtype == 1:
        
            x = colors['J-H']
            y = mags['H']
        
            xlabel = 'J - H'
            ylabel = 'H'
        
            invert_y = True
        
        
        elif dtype == 2:
        
            x = colors['Paβ-H']
            y = mags['H']
        
            xlabel = r'Pa$\beta$ - H'
            ylabel = 'H'
        
            invert_y = True
        
        
        elif dtype == 3:
        
            x = colors['FeII-H']
            y = mags['H']
        
            xlabel = 'Fe II - H'
            ylabel = 'H'
        
            invert_y = True
        
        
        elif dtype == 4:
        
            x = colors['J-H']
            y = colors['Paβ-H']
        
            xlabel = 'J - H'
            ylabel = r'Pa$\beta$ - H'
        
            invert_y = False
        
        
        elif dtype == 5:
        
            x = colors['J-H']
            y = colors['FeII-H']
        
            xlabel = 'J - H'
            ylabel = 'Fe II - H'
        
            invert_y = False
        
        
        elif dtype == 6:
        
            x = colors['Paβ-H']
            y = colors['FeII-H']
        
            xlabel = r'Pa$\beta$ - H'
            ylabel = 'Fe II - H'
        
            invert_y = False
        
        
        elif dtype == 7:
        
            x = colors['J-H']
            y = colors['Paβ-FeII']
        
            xlabel = 'J - H'
            ylabel = r'Pa$\beta$ - Fe II'
        
            invert_y = False
        
        
        else:
            raise ValueError('Unknown dtype')
        
        
        # =============================================================================
        # Plot
        # =============================================================================
        
        fig, ax = plt.subplots(figsize=(4,4))
        
        ax.set_title(f'{zone} ★ = {len(cat2m)}')
        
        h = ax.hexbin(
            x,
            y,
            gridsize=(14,10),
            norm=LogNorm(vmin=0.1),
            cmap='Blues',
            edgecolor='grey',
            lw=0.1
        )
        
        plt.colorbar(h, ax=ax, aspect=40)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        if invert_y:
            ax.invert_yaxis()
        
        
        # =============================================================================
        # Overplot cluster members
        # =============================================================================
        
        clus = Table.read(tmp2 + 'clus1_HJ.txt', format='ascii')
        
        mask = np.isin(cat2m['id'], clus['id'])
        
        ax.scatter(
            x[mask],
            y[mask],
            color='lime',
            marker='.',
            s=50,
            edgecolor = 'k',lw = 0.3,
            label='DBSCAN cluster'
        )
        
        ax.legend(loc=2)
        
        plt.tight_layout()
        plt.show()
# # =============================================================================
# #     Color-Color diagram ( Chose one)
# # =============================================================================
    

#     # max_sep = 0.1 * u.arcsec  # choose an appropriate matching radius
    
#     # Match 110W -> 160W
#     idx12, d12, _ = crd2.match_to_catalog_sky(crd1)
#     m12 = d12 < max_sep
    
#     # Match 128N -> 160W
#     idx32, d32, _ = crd2.match_to_catalog_sky(crd3)
#     m32 = d32 < max_sep
    
#     # Match 128N -> 160W
#     idx42, d42, _ = crd2.match_to_catalog_sky(crd4)
#     m42 = d42 < max_sep
    
#     # Stars detected in all three bands
#     good = m12 & m32 & m42
    
#     cat2m = cat2[good]
#     cat1m = cat1[idx12[good]]
#     cat3m = cat3[idx32[good]]
#     cat4m = cat4[idx42[good]]
    
#     # col1 = cat110[''] - cat160['mag']    # F110W - F160W
#     # col2 = cat['mag'128['mag'] - cat160['mag']   # F128N0W
    
#     fig, ax = plt.subplots(1,1, figsize =(4,4)) 
#     ax.set_title(f'{zone} \u2605 = {len(cat2m)}')
#     # h = ax.hexbin(cat1m[b_name1] - cat2m[b_name2], cat3m[b_name3], norm = LogNorm(), 
#     #              gridsize = (14,10),cmap = 'Greys_r', edgecolor = None,lw= 0.1 )
    
#     # h = ax.hexbin(cat1m[b_name1] - cat2m[b_name2], cat3m[b_name3] - cat2m[b_name2], norm = LogNorm(vmin = 0.1), 
#     #               gridsize = (14,10),cmap = 'Blues', edgecolor = 'grey',lw= 0.1 )
    
#     h = ax.hexbin(cat1m[b_name1] - cat2m[b_name2], cat3m[b_name3] - cat4m[b_name4], norm = LogNorm(vmin = 0.1), 
#                   gridsize = (14,10),cmap = 'Blues', edgecolor = 'grey',lw= 0.1 )
    
    
#     cbar = plt.colorbar(h, ax=ax, aspect = 40)
    
    
#     # h = ax.hist2d(cat1m[b_name1] - cat2m[b_name2], cat3m[b_name3] - cat2m[b_name2], norm = LogNorm(), bins = (20,20))
#     # cbar = plt.colorbar(h[3], ax=ax, aspect = 40)
    
#     # ax.axvline(color_l, color = 'red', ls = 'dashed',label = f'J-H = {color_l}')
#     # ax.invert_yaxis()
#     ax.set_xlabel(f'{b_name1} - {b_name2}')
#     ax.set_ylabel(f'{b_name3} - {b_name4}')
#     ax.legend() 
    
#     # clus = Table.read(tmp2  + 'clus1_HJ.txt', format = 'ascii')
#     clus = Table.read(tmp2  + 'clus0_HJ.txt', format = 'ascii')
    
    
#     mask = np.isin(cat2m['id'], clus['id'])
#     maskR = np.isin(clus['id'], cat2m['id'])
    
#     cat2m_c = cat2m[mask]
#     cat1m_c = cat1m[mask]
#     cat3m_c = cat3m[mask]
#     cat4m_c = cat4m[mask]
    
#     # ax.scatter(cat1m_c[b_name1] - cat2m_c[b_name2], cat3m_c[b_name3], 
#     #            color = 'red', marker = 'x',s=1)
#     ax.scatter(cat1m_c[b_name1] - cat2m_c[b_name2], cat3m_c[b_name3] - cat4m_c[b_name4], 
#                color = 'lime',marker = '.', edgecolor = 'lime',facecolor = 'none', label = 'Dbscan cluster')
#     ax.legend(loc=2)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



