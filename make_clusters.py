#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 13:10:36 2026

@author: amartinez
"""

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from spisea import synthetic, evolution, atmospheres, reddening
from spisea.imf import imf, multiplicity
from matplotlib.colors import LogNorm

# ------------------------------------------------------------------
# Step 1: Isochrone
# ------------------------------------------------------------------
# logAge = np.log10(5.01 * 10**5)
logAge = np.log10(5.01 * 10**5)
AKs = 0.84      # extinction in Ks (mag)
delta_AKs = 0.1
dist = 5700                      # distance (pc)
metallicity = 1                   # [M/H]

evo_model = evolution.MISTv1()
atm_func = atmospheres.get_merged_atmosphere
red_law = reddening.RedLawHosek18b()

filt_list = ['wfc3,ir,f110w', 'wfc3,ir,f160w' ]
# filt_list = ['wfc3,ir,f160w' ]
iso_dir = 'isochrones/'

my_iso = synthetic.IsochronePhot(logAge, AKs, dist,
                                 metallicity=metallicity,
                                 evo_model=evo_model,
                                 atm_func=atm_func,
                                 red_law=red_law,
                                 filters=filt_list,
                                 iso_dir=iso_dir)

# ------------------------------------------------------------------
# Step 2: IMF with Lu+13 unresolved multiplicity (Kroupa 2001)
# ------------------------------------------------------------------

# imf_use = 'th'
imf_use = 'Kr'


imf_multi = multiplicity.MultiplicityUnresolved()
massLimits = np.array([0.2, 0.5, 1, 120])

# Kroupa 2001
if imf_use == 'Kr':
    powers = np.array([-1.3, -2.3, -2.3])
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)


# Top-heavy: high-mass segments flattened (dN/dm ~ m^-1.3 everywhere)
if imf_use == 'toph':
    powers = np.array([-1.3, -1.3, -1.3])
    my_imf = imf.IMF_broken_powerlaw(massLimits, powers, imf_multi)
# ------------------------------------------------------------------
# Step 3: Cluster of 1000 Msun
# ------------------------------------------------------------------
mass = 50.0

# Differential extinction: 1-sigma spread of the AKs distribution (mag).
# The cluster is built by drawing each star's AKs from a normal
# distribution centred on the isochrone AKs with this dispersion.


cluster = synthetic.ResolvedClusterDiffRedden(my_iso, my_imf, mass,
                                             deltaAKs=delta_AKs)

clust = cluster.star_systems
iso = my_iso.points

print('Number of star systems: {0}'.format(len(clust)))



print('Number of star systems: {0}'.format(len(clust)))
# %%


# # ------------------------------------------------------------------
# # CMD: cluster vs. isochrone
# # ------------------------------------------------------------------
# fig, ax = plt.subplots(1, 1, figsize=(10, 10))
# ax.plot(clust['m_hst_f110w'] - clust['m_hst_f160w'],
#         clust['m_hst_f160w'],
#         'k.', ms=20, alpha=1, label='__nolegend__')
# ax.plot(iso['m_hst_f110w'] - iso['m_hst_f160w'],
#         iso['m_hst_f160w'],
#         'r-', label='Isochrone')
# ax.set_xlabel('F110W - F160W')
# ax.set_ylabel('F160W')
# ax.invert_yaxis()
# ax.legend()
# plt.savefig('cluster_cmd.png', dpi=150)
# plt.show()

# ------------------------------------------------------------------
# CMD: cluster vs. isochrone, 1 Msun marked
# ------------------------------------------------------------------
clust = cluster.star_systems
iso = my_iso.points

# Index of the isochrone point closest to 1 Msun
idx = np.where(abs(iso['mass'] - 1.0) ==
               min(abs(iso['mass'] - 1.0)))[0]

fig, ax = plt.subplots(1, 1, figsize=(10, 10))
h = ax.hexbin(clust['m_hst_f110w'] - clust['m_hst_f160w'],
        clust['m_hst_f160w'],alpha=1, label='__nolegend__', norm = LogNorm(), cmap = 'grey')

cbar = plt.colorbar(h, ax=ax,)  
# ax.plot(clust['m_hst_f110w'] - clust['m_hst_f160w'],
#         clust['m_hst_f160w'],
#         'k.', ms=10, alpha=0.2, label='__nolegend__')
ax.plot(iso['m_hst_f110w'] - iso['m_hst_f160w'],
        iso['m_hst_f160w'],
        'r-', label='Isochrone')
# ax.plot((iso['m_hst_f110w'][idx] - iso['m_hst_f160w'][idx]),
#         iso['m_hst_f160w'][idx],
#         'w*', ms=15, label='1 $M_\odot$')
ax.scatter((iso['m_hst_f110w'][idx] - iso['m_hst_f160w'][idx]),
        iso['m_hst_f160w'][idx], marker = '*',s =400,
        color = 'tab:orange', label='1 $M_\odot$', zorder = 3)
ax.set_xlim(0,4)
ax.set_xlabel('F110W - F160W')
ax.set_ylabel('F160W')
ax.invert_yaxis()
ax.legend()
plt.savefig('cluster_cmd.png', dpi=150)
plt.show()

folder_sim = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/simulated_clusters/'
clust.write(folder_sim  + f'clus_{mass:.0f}M_Age{10**logAge/10**6:.1f}Myr_imf{imf_use}_{dist}pc_AKs{AKs:.2f}_dAKs{delta_AKs}.txt', format = 'ascii')

















