
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
from matplotlib.colors import LogNorm
import cluster_finder
from hst_irZP import get_vegazp
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
# HST observation
# =============================================================================

pixSca = 0.12825 #arcsec/pixel
# pixSca = 0.00001 #arcsec/pixel#

e_pos_cat  = 0.05# in arcsec. The position errors from starfinder lists are largely overstatimated!!!

red_techn = 'Gaia'
# red_techn = 'Original'



for epoch in range(1,3):
    results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
    
    if red_techn == 'Gaia':
        cat = Table.read(results + f'calib_{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii')
        ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_f{band}_drz_sci.fits')
        wcs = WCS(ima[0].header)
        scale_pix = np.sqrt(ima[0].header['CD1_1']**2 + ima[0].header['CD2_1']**2 )*3600
        
    # if red_techn == 'Original':
    #     cat = Table.read(results + f'hst_ep{epoch}_f{band}_drz_stars{band}.txt', format = 'ascii')
    #     ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/hst_ep{epoch}_f{band}_drz.fits')
    #     wcs = WCS(ima[1].header)
    #     scale_pix = np.sqrt(ima[1].header['CD1_1']**2 + ima[1].header['CD2_1']**2 )*3600
    
    print(pixSca, scale_pix)
    
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
    
    cat['x'] = cat['x']*pixSca*-1 # These are arcsec
    cat['y'] = cat['y']*pixSca
    
    cat['sx'] = cat['sx']*pixSca # These are arcsec
    cat['sy'] = cat['sy']*pixSca
    
    cat['sxy'] = np.sqrt(cat['sx']**2 + cat['sy']**2)
    cat = cat[cat['sxy'] > 0]
    
    fig, (ax, ax2) = plt.subplots(1,2, figsize = (12,6))
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
   
    

    # cat = filter_hst_data(cat, max_e_pos = e_pos_cat)
   
# %%
  

fig, ax = plt.subplots(figsize =(10,5))
ax.set_title(f'Zone: {zone}')
num_bins = np.arange(9,22,0.5)
# Create the histogram
mag1 = cat['H']
# mag2 = cat['H']

counts1, bin_edges1 = np.histogram(mag1, bins=num_bins)
# counts2, bin_edges2 = np.histogram(mag2, bins=num_bins)
c1_mask = counts1 > 1
# c2_mask = counts2 > 10



bin_1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2
# bin_2 = (bin_edges1[:-1] + bin_edges1[1:]) / 2
# bin_2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

counts1 = counts1[c1_mask]
# counts2 = counts2[c2_mask]

bin_1 = bin_1[c1_mask]
# bin_2 = bin_2[c2_mask]
# '#1f77b4', '#ff7f0e'


    
ax.plot(bin_1, counts1, marker = '.', label = 'GNS$\,$I DR1')
# ax.plot(bin_2, counts2, marker ='.',label = 'GNS$\,$I DR2')
# ax.plot(bin_1, counts1, marker = '.', label = 'GNSI: %s stars'%(len(gns1)))
# ax.plot(bin_2, counts2, marker ='.',label = 'GNSII: %s stars'%(len(gns2)))
err1 = np.sqrt(counts1)
# err2 = np.sqrt(counts2)

ax.errorbar(bin_1, counts1, yerr=err1, fmt='.', markersize=5, capsize=3,c ='tab:blue' )
# ax.errorbar(bin_2[0:-13], counts2[0:-13], yerr=err2[0:-13], fmt='.', markersize=5, capsize=3,c ='#ff7f0e' )
ax.set_yscale('log')    
ax.set_xticks(np.arange(10,23,1))
# ax.grid()  
ax.grid(True, which="both", ls="--")    
ax.set_xlabel(f'[{band}]')  
ax.set_ylabel('#')  
       




















        


