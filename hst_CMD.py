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
# results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
# results = ''
# zone = 'G032.03+00.05'
# zone = 'G028.20-00.05'
zones = [ 'AFGL5180','G028.20-00.05', 'G032.03+00.05', 'G35.2-0.74N', 'G339.88-01.26','IRAS07299-1651', 'IRAS16562-3959']
# zones = [ 'IRAS07299-1651']
color_l = 1
# zone = zones[2]

for zone in zones:
    # band = '160w'
    band1 = '110w'
    band2 = '160w'
    band3 = '128n'
    # band = '128n'
    # epoch = 2
    if band1 == '160w':
        b_name1 = 'H'
    elif band1 == '110w':
        b_name1 = 'J'
    elif band1 == '128n':
        b_name1 = 'Paβ'
    elif band1 == '164N':
        b_name1 = 'Fe II'
    if band2 == '160w':
        b_name2 = 'H'
    elif band2 == '110w':
        b_name2 = 'J'
    elif band2 == '128n':
        b_name2 = 'Paβ'
    elif band2 == '164N':
        b_name2 = 'Fe II'
        
    # =============================================================================
    # Alig Para
    # =============================================================================
    max_sep = 50*u.mas
    
    
    pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'
    tmp1 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band1}/'
    tmp2 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band2}/'
    tmp3 = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band3}/'
    
    
    
    cat1 = Table.read(tmp1 +  f'pm_{zone}_{band1}_posEpo1.txt', format = 'ascii')
    cat2 = Table.read(tmp2 +  f'pm_{zone}_{band2}_posEpo1.txt', format = 'ascii')
    cat3 = Table.read(tmp3 +  f'pm_{zone}_{band3}_posEpo1.txt', format = 'ascii')

    
    
    
    crd1 = SkyCoord(ra = cat1['ra'], dec = cat1['dec'], unit = 'degree', frame = 'fk5')
    crd2 = SkyCoord(ra = cat2['ra'], dec = cat2['dec'], unit = 'degree', frame = 'fk5')
    crd3 = SkyCoord(ra = cat3['ra'], dec = cat3['dec'], unit = 'degree', frame = 'fk5')
    
    
    idx, d2d, _ = crd1.match_to_catalog_sky(crd2, nthneighbor=1)
    match_mask = d2d < max_sep
    cat1_m = cat1[match_mask]
    cat2_m = cat2[idx[match_mask]]
    print(40*'+')
    unicos = unique(cat1_m, keep = 'first')
    print(len(cat1_m),len(unicos))
    print(40*'+')
    
    clust = os.path.exists(tmp1 +'clust_0.')
    
    fig, ax = plt.subplots(1,1, figsize =(4,4)) 
    ax.set_title(f'{zone} \u2605 = {len(cat1_m)}')
    h = ax.hexbin(cat1_m[b_name1] - cat2_m[b_name2], cat2_m[b_name2], norm = LogNorm(),cmap = 'viridis' )
    cbar = plt.colorbar(h, ax=ax, aspect = 40)
    ax.axvline(color_l, color = 'red', ls = 'dashed',label = f'J-H = {color_l}')
    ax.invert_yaxis()
    ax.set_xlabel(f'{b_name1} - {b_name2}')
    ax.set_ylabel(f'{b_name2}')
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
#     Color-Color digram
# =============================================================================
    

    # max_sep = 0.1 * u.arcsec  # choose an appropriate matching radius
    
    # Match 110W -> 160W
    idx12, d12, _ = crd2.match_to_catalog_sky(crd1)
    m12 = d12 < max_sep
    
    # Match 128N -> 160W
    idx32, d32, _ = crd2.match_to_catalog_sky(crd3)
    m32 = d32 < max_sep
    
    # Stars detected in all three bands
    good = m12 & m32
    
    cat160 = cat2[good]
    cat110 = cat1[idx12[good]]
    cat128 = cat3[idx32[good]]
    
    # col1 = cat110[''] - cat160['mag']    # F110W - F160W
    # col2 = cat['mag'128['mag'] - cat160['mag']   # F128N0W
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



