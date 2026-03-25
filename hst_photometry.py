
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
import matplotlib.pyplot as plt
# %%plotting parametres
from matplotlib import rc
from matplotlib import rcParams

# plt.rcParams["mathtext.fontset"] = 'dejavuserif'
# rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rcParams.update({'figure.max_open_warning': 0})# a warniing for matplot lib pop up because so many plots, this turining it of
# # Enable automatic plotting mode
# # IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')





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
max_sep = 50*u.mas

# =============================================================================
# Pho quality
# =============================================================================
max_eH = 0.05
max_H = 10 
min_H = 15

pho = Table.read(phot_folder + f'{zone}_2MASS.tsv', format = 'ascii')
cat = Table.read(results + f'{zone}_EP{epoch}_f{band}_drz_sci_stars{band}.txt', format = 'ascii')
ima = fits.open(folder  + f'{zone}/gaia_alignment/Epoch{epoch}/{zone}_EP{epoch}_f{band}_drz_sci.fits')
cabeza = ima[0].header
wcs = WCS(ima[0].header)

# =============================================================================
# ZP calcultion (from: https://spacetelescope.github.io/hst_notebooks/notebooks/WFC3/zeropoints/zeropoints.html)
# =============================================================================

import tarfile
from synphot import Observation

# tar_archive = 'hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar'
# extract_to = 'hlsp_reference-atlases_hst_multi_everything_multi_v11_sed'
tar_archive = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/synphot_cdbs/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed.tar'
extract_to = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/synphot_cdbs/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed'
# %%
with tarfile.open(tar_archive, 'r') as tar:

    tar.extractall(path=extract_to, filter='data')

# os.environ['PYSYN_CDBS'] = 'hlsp_reference-atlases_hst_multi_everything_multi_v11_sed/grp/redcat/trds/'
os.environ['PYSYN_CDBS'] = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/synphot_cdbs/hlsp_reference-atlases_hst_multi_everything_multi_v11_sed/grp/redcat/trds/'

# %%

import stsynphot as stsyn
# %%

vega_url = 'https://ssb.stsci.edu/trds/calspec/alpha_lyr_stis_010.fits'
stsyn.Vega = stsyn.spectrum.SourceSpectrum.from_file(vega_url)
# %%
# detectors = ['uvis1']
# detectors = ['uvis1', 'uvis2']
detectors = ['ir']
# %%

# filtnames = ['f200lp', 'f218w', 'f225w', 'f275w', 'f280n', 'f300x', 'f336w', 'f343n', 'f350lp', 
#              'f373n', 'f390m', 'f390w', 'f395n', 'f410m', 'f438w', 'f467m', 'f469n', 'f475w', 
#              'f475x', 'f487n', 'f502n', 'f547m', 'f555w', 'f600lp', 'f606w', 'f621m', 'f625w', 
#              'f631n', 'f645n', 'f656n', 'f657n', 'f658n', 'f665n', 'f673n', 'f680n', 'f689m', 
#              'f763m', 'f775w', 'f814w', 'f845m', 'f850lp', 'f953n']
# filtnames = ['f160w']   
filtnames = [f'f{band}']   
# filtnames = ['f098m', 'f105w', 'f110w', 'f125w', 'f126n', 'f127m', 'f128n', 'f130n', 
#              'f132n', 'f139m', 'f140w', 'f153m', 'f160w', 'f164n', 'f167n']

# %%
mjd = ima[0].header['EXPSTART']   # e.g., 57610.4918468


# mjd = '60755.92021683'
# mjd = str(Time.now().mjd)
# %%

aper = '6.0'
# aper = '0.396'
# aper = '0.385'
# %%

# obsmode = f'wfc3,{detectors[0]},{filtnames[0]},mjd#{mjd},aper#{aper}'
obsmode = f'wfc3,{detectors[0]},{filtnames[0]},mjd#{mjd},aper#{aper}'
bp = stsyn.band(obsmode)
# %%
def calculate_values(detector, filt, mjd, aper):
    # parameters can be removed from obsmode as needed
    obsmode = f'wfc3, {detector}, {filt}, mjd#{mjd}, aper#{aper}'
    bp = stsyn.band(obsmode)  
    
    # STMag
    photflam = bp.unit_response(stsyn.conf.area)  # inverse sensitivity in flam
    stmag = -21.1 - 2.5 * np.log10(photflam.value)
    
    # Pivot Wavelength and bandwidth
    photplam = bp.pivot() # pivot wavelength in angstroms
    bandwidth = bp.photbw() # bandwidth in angstroms
    
    # ABMag
    abmag = stmag - 5 * np.log10(photplam.value) + 18.6921
    
    # Vegamag
    obs = Observation(stsyn.Vega, bp, binset=bp.binset)  # synthetic observation of vega in bandpass using vega spectrum
    vegamag = -1 * obs.effstim(flux_unit='obmag', area=stsyn.conf.area)
    
    return obsmode, photplam.value, bandwidth.value, photflam.value, stmag, abmag, vegamag.value
# %%
obsmode, photplam, bandwidth, photflam, stmag, abmag, ZP_Vega = calculate_values(detectors[0], filtnames[0], mjd, aper)

# print values
# print('                                      ')
print(f'Obsmode {obsmode}\nPivotWave{photplam:.1f}\nPhotflam {photflam:.4e}\nSTMAG {stmag:.3f}\nABMAG {abmag:.3f}\nVEGAMAG {ZP_Vega:.3f}')



cat['H'] =  -2.5*np.log10(cat['f']) + ZP_Vega
cat['dH'] = (2.5 / np.log(10)) * (cat['sf'] / cat['f'])




cat_rd = wcs.pixel_to_world(cat['x'], cat['y'])
if epoch == 1:
    center = SkyCoord(ra = np.mean(cat_rd.ra.value), dec = np.mean(cat_rd.dec.value),
                      unit = 'degree', frame = 'icrs')

cat['ra'] = cat_rd.ra
cat['dec'] = cat_rd.dec
# %%


fig, (ax, ax2) = plt.subplots(1,2, figsize = (10,5))
ax.set_title('2MASS')
ax.scatter(pho['Hmag'], pho['e_Hmag'], facecolor = 'none', edgecolor = 'k')
ax.set_ylim(0,0.3)
ax.set_xlabel('H')
ax.set_ylabel('err_H')
ax.axhline(max_eH, ls = 'dashed', color = 'red')
ax.axvline(max_H, ls = 'dashed', color = 'red')
ax.axvline(min_H, ls = 'dashed', color = 'red')
ax2.set_title('HST')
# ax2.scatter(cat['H'], cat['dH'])
h = ax2.hexbin(cat['H'], cat['dH'], norm = LogNorm())
plt.colorbar(h, ax = ax2                                                                                                                                                                                                                                                                                                                                                                                                                                    )
ax2.set_xlabel('H')
ax2.set_ylabel('err_H')


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


dH_cat = cat_m['H'] - phot_m['Hmag']
m_H = sigma_clip(dH_cat,3, masked=True, return_bounds=True)
dH_m = dH_cat[np.logical_not(m_H[0].mask)]
cat_mask = cat_m[np.logical_not(m_H[0].mask)]

fig, ax = plt.subplots(1,1,figsize = (7,7))
ax.axhline(np.mean(dH_m), color = 'tab:orange', label = f'$\overline{{\Delta H}}$ = {np.mean(dH_m): .2f} ')
ax.axhline(m_H[2], ls = 'dashed' ,color = 'r', label = f'$\pm\sigma$ = {np.std(dH_m): .2f} ')
ax.axhline(m_H[1], ls = 'dashed' ,color = 'r')
ax.scatter(cat_mask['H'] , dH_m)
# h = ax.hexbin(cat_m['H'] , dH_cat, norm = LogNorm())
# cbar = plt.colorbar(h, ax=ax)
ax.legend()
ax.set_ylabel('$\Delta H$', fontsize = 12)
ax.set_xlabel('[H] mag', fontsize = 12)



# %%


