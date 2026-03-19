
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
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams

plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
# Enable automatic plotting mode
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')





# zone = 'G032.03+00.05'
zone = 'G028.20-00.05'
band = '160w'
# band = '128n'
epoch = 1
results = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/f{band}/epoch{epoch}/'
folder = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/'

phot_folder = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/{zone}/'

pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'

# =============================================================================
# Matching
# =============================================================================
max_sep = 200*u.mas

# =============================================================================
# Pho quality
# =============================================================================
max_eH = 0.05
max_H = 10 
min_H = 15

pho = Table.read(phot_folder + f'{zone}_2MASS.tsv', format = 'ascii')
cat = Table.read(results + f'{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii')
ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_f{band}_drz_sci.fits')
wcs = WCS(ima[0].header)
cat_rd = wcs.pixel_to_world(cat['x'], cat['y'])

if epoch == 1:
    center = SkyCoord(ra = np.mean(cat_rd.ra.value), dec = np.mean(cat_rd.dec.value),
                      unit = 'degree', frame = 'icrs')

cat['ra'] = cat_rd.ra
cat['dec'] = cat_rd.dec

fig, ax = plt.subplots(1,1)
ax.scatter(pho['Hmag'], pho['e_Hmag'])
ax.set_ylim(0,0.3)
ax.set_xlabel('H')
ax.set_ylabel('err_H')
ax.axhline(max_eH, ls = 'dashed', color = 'red')
ax.axvline(max_H, ls = 'dashed', color = 'red')
ax.axvline(min_H, ls = 'dashed', color = 'red')


mask_pho = (pho['e_Hmag'] < max_eH) & (pho['Hmag'] > max_H) & (pho['Hmag'] < min_H)
pho = pho[mask_pho]

ax.scatter(pho['Hmag'], pho['e_Hmag'])


plt.show()
pho_c = SkyCoord(ra = pho['RAJ2000'], dec = pho['DEJ2000'],unit = 'degree')
cat_c = SkyCoord(ra = cat['ra'], dec = cat['dec'], frame='icrs')

idx, d2d, _ = pho_c.match_to_catalog_sky(cat_c, nthneighbor=1)
match_mask = d2d < max_sep
phot_m = pho[match_mask]
cat_m = cat[idx[match_mask]]
print(40*'+')
unicos = unique(cat_m, keep = 'first')
print(len(cat_m),len(unicos))
print(40*'+')

region(cat_m, 'ra', 'dec',
       name = f'{zone}_{band}_with2MASS',
       save_in = pruebas,
       color = 'blue',
       marker = 'cross')


hins = -2.5*np.log10(cat_m['f'])


ZPs = phot_m['Hmag'] - hins

plt.plot()
plt.hist(ZPs, bins = 'auto')












