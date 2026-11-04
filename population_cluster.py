# =============================================================================
# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Mon Nov  2 11:29:26 2026
# 
# @author: amartinez
# """
# import sys
# # sys.path.append("/Users/amartinez/Desktop/pythons_imports/")
# import numpy as np
# from astropy import units as u
# import sys
# import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
# from astropy.table import Table, vstack, Column
# 
# import cluster_finder
# 
# from astropy.table import Table
# 
# 
# # %%
# folder_clus = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/simulated_clusters/'
# cat_zone = 'G028.20-00.05'
# folder_cat = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{cat_zone}/f160w/'
# 
# mass = 50
# Age = 0.5
# imf = 'Kr'
# dist = 5700
# AKs = 0.84
# dAKs = 0.1
# cluster = Table.read(folder_clus + f'clus_{mass}M_Age{Age}Myr_imf{imf}_{dist}pc_AKs{AKs}_dAKs{dAKs}.txt', format = 'ascii')
# cat = Table.read(folder_cat + f'pm_{cat_zone}_160w_posEpo1.txt', format = 'ascii')
# 
# 
# 
# #!/usr/bin/env python3
# """
# Inject a synthetic SPISEA cluster into a proper-motion catalog,
# after scrambling the catalog to destroy real clusters.
# 
# Outputs:
#   cat_scrambled_with_cluster.fits  - catalog to feed the finder
#   fake_cluster_members.fits        - ground truth (rows of the fake cluster)
# """
# 
# 
# 
# for loop in range(2):
#     
#     rng = np.random.default_rng()
#     
#     # Tunable parameters -------------------------------------------------
#     SIGMA_POS_ARCSEC = 5   # spatial dispersion of cluster members (")
#     sigma_pm = 1  # mas/yr
#     # --------------------------------------------------------------------
#     
#     n_cat = len(cat)
#     n_cl = len(cluster)
#     
#     # ==================================================================
#     # STEP 1: Scramble the real catalog
#     # ==================================================================
#     # Permute the spatial block and the kinematic block independently.
#     # Each block is permuted as a unit so that (x, y, ra, dec) stay
#     # mutually consistent within a row, and likewise the PMs.
#     pos_block = ['x', 'y', 'ra', 'dec', 'xp', 'yp']
#     pm_block  = ['pm_x', 'pm_y', 'dpm_x', 'dpm_y']
#     
#     perm_pos = rng.permutation(n_cat)
#     perm_pm  = rng.permutation(n_cat)
#     
#     cat_s = cat.copy()
#     for c in pos_block:
#         cat_s[c] = cat[c][perm_pos]
#     for c in pm_block:
#         cat_s[c] = cat[c][perm_pm]
#     
#     # After this, no row retains its original (position, PM) pairing:
#     # real spatial clusters lose their kinematic coherence and vice versa.
#     
#     # ==================================================================
#     # STEP 2: Choose the cluster's position and mean PM from the catalog
#     # ==================================================================
#     anchor = rng.integers(n_cat)          # spatial anchor: a random real star
#     ra0, dec0 = cat['ra'][anchor], cat['dec'][anchor]
#     
#     # PM anchor: a random real star (keeps PMs inside the data set's range)
#     pm_anchor = rng.integers(n_cat)
#     pmx0, pmy0 = cat['pm_x'][pm_anchor], cat['pm_y'][pm_anchor]
#     
#     # Internal dispersion in mas/yr from a physical velocity dispersion:
#     
#     
#     
#     # ==================================================================
#     # STEP 3: Generate cluster positions and PMs
#     # ==================================================================
#     # Positions: isotropic Gaussian on the sky around the anchor.
#     # RA offsets are divided by cos(dec) so the angular dispersion is
#     # isotropic in both coordinates.
#     ra_cl  = ra0  + (rng.normal(0, SIGMA_POS_ARCSEC, n_cl) / 3600.0) / np.cos(np.radians(dec0))
#     dec_cl = dec0 +  rng.normal(0, SIGMA_POS_ARCSEC, n_cl) / 3600.0
#     
#     # Proper motions: Gaussian around the anchor PM
#     pmx_cl = pmx0 + rng.normal(0, sigma_pm, n_cl)
#     pmy_cl = pmy0 + rng.normal(0, sigma_pm, n_cl)
#     
#     # ==================================================================
#     # STEP 4: Convert cluster RA/Dec to detector coordinates (x, y)
#     # ==================================================================
#     # Fit an affine transform from the catalog itself, so the cluster
#     # lands on the same detector geometry as the real stars.
#     A = np.column_stack([cat['ra'], cat['dec'], np.ones(n_cat)])
#     cx, *_ = np.linalg.lstsq(A, cat['x'], rcond=None)
#     cy, *_ = np.linalg.lstsq(A, cat['y'], rcond=None)
#     
#     A_cl = np.column_stack([ra_cl, dec_cl, np.ones(n_cl)])
#     x_cl  = A_cl @ cx
#     y_cl  = A_cl @ cy
#     
#     # ==================================================================
#     # STEP 5: Fill the remaining catalog columns for cluster members
#     # ==================================================================
#     # Errors and quality flags: draw from the catalog's own distributions,
#     # by sampling random catalog rows (keeps error-error correlations).
#     donors = rng.integers(0, n_cat, n_cl)
#     
#     fake = Table()
#     fake['id']   = np.max(cat['id']) + 1 + np.arange(n_cl)
#     fake['ra']   = ra_cl*u.deg
#     fake['dec']  = dec_cl*u.deg
#     fake['x']    = x_cl
#     fake['y']    = y_cl
#     fake['pm_x'] = pmx_cl
#     fake['pm_y'] = pmy_cl
#     
#     # Magnitudes: F160W is a good proxy for H, F110W for J (approximation)
#     fake['H'] = cluster['m_hst_f160w']
#     fake['J'] = cluster['m_hst_f110w']
#     
#     # Copy photometric/astrometric error structure from random real rows
#     for c in ['f', 'sf', 'dH', 'dJ', 'sxy', 'dpm_x', 'dpm_y',
#               'xp', 'yp', 'sx', 'sy']:
#         fake[c] = cat[c][donors]
#     
#     # Reorder columns to match the catalog exactly
#     fake = fake[list(cat.colnames)]
#     
#     # ==================================================================
#     # STEP 6: Assemble and write outputs
#     # ==================================================================
#     out = vstack([cat_s, fake])
#     out['is_fake'] = np.zeros(len(out), dtype=int)
#     out['is_fake'][len(cat_s):] = 1     # injected members are flagged 1
#     
#     # out.write('cat_scrambled_with_cluster.fits', overwrite=True)
#     
#     truth = fake['id', 'ra', 'dec', 'pm_x', 'pm_y', 'H']
#     # truth = fake['id', 'ra', 'dec', 'pm_x', 'pm_y', 'H', 'mass']
#     # truth.write('fake_cluster_members.fits', overwrite=True)
#     
#     print('Injected {0} fake members among {1} total catalog entries.'.format(n_cl, len(out)))
#     print('Ground truth saved to fake_cluster_members.fits')
#     
#     # stop(201)
#     # %%
#     # import cluster_finder
#     cmap_color = 'grey'
#     fig, ax = plt.subplots(1,3, figsize = (15,5))
#     ax[1].hexbin(cat_s['ra'], cat_s['dec'], norm = LogNorm(), cmap = cmap_color)
#     ax[0].hexbin(cat_s['pm_x'], cat_s['pm_y'], norm = LogNorm(), cmap = cmap_color)
#     ax[2].hexbin(cat_s['J'] - cat_s['H'], cat_s['H'], norm = LogNorm(), cmap = cmap_color)
#     
#     ax[1].scatter(fake['ra'], fake['dec'], norm = LogNorm(), color = 'lime', edgecolor = 'k')
#     ax[0].scatter(fake['pm_x'], fake['pm_y'], norm = LogNorm(), color = 'lime', edgecolor = 'k')
#     ax[2].scatter(fake['J'] - fake['H'],fake['H'], norm = LogNorm(), color = 'lime', edgecolor = 'k')
#     ax[2].invert_yaxis()
#     ax[1].invert_xaxis()
#     
#     
#     look_for_cluster = 1
#     if look_for_cluster:
#         
#        
#         
#         modes = ['pm_xy_color']
#         # modes = ['pm_xy']
#         knn = 15
#         gen_sim = 'kernnel'
#         sim_lim ='minimun'
#         # sim_lim ='mean'
#     
#         clus_dic = cluster_finder.finder(out, 'pm_x', 'pm_y',
#                                      'ra', 'dec', 
#                                      'ra', 'dec',
#                                      modes[0],
#                                      'J','H',
#                                      knn,gen_sim,sim_lim, save_reg = None)
#     
# # %%
#     
# 
#     import numpy as np
#     from sklearn.metrics import f1_score, confusion_matrix
#     
#     import numpy as np
#     from sklearn.metrics import f1_score, confusion_matrix, recall_score, precision_score
#     
#     for cnum in range(len(clus_dic)):
#         truth = out['is_fake'] == 1
#         pred  = np.isin(out['id'], clus_dic[f'clus_{cnum}']['id'])   # your finder's member list
#         
#         # Full confusion matrix
#         tn, fp, fn, tp = confusion_matrix(truth, pred).ravel()
#         
#         print('Injected members:  {0}'.format(tp + fn))
#         print('Recovered (TP):    {0}'.format(tp))
#         print('Missed (FN):       {0}'.format(fn))
#         print('False positives:   {0}'.format(fp))
#         # Percentage of recovered stars: recall
#         recovered_pct = 100.0 * tp / (tp + fn)
#         print('Recovery rate:     {0:.1f}%'.format(recovered_pct))
#         print(30*'#')
# 
# 
# 
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:29:26 2026

@author: amartinez

Injection-recovery experiment:
  - scramble the real catalog to destroy real clusters
  - inject a synthetic SPISEA cluster with random PM dispersion
    (0.2 - 2.0 mas/yr) and random position dispersion (2 - 15 arcsec)
  - run the cluster finder and measure recovery statistics
  - accumulate the statistics over many realisations and summarise by bin

Outputs:
  injection_results.fits  - one row per realisation
"""
import sys
sys.path.append("/Users/amartinez/Desktop/pythons_imports/")
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.table import Table, vstack, Column

import cluster_finder
from sklearn.metrics import confusion_matrix

# %%
folder_clus = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/simulated_clusters/'
cat_zone = 'G028.20-00.05'
folder_cat = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{cat_zone}/f160w/'

pruebas = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/pruebas/'
mass = 50
Age = 0.5
imf = 'Kr'
dist = 5700
AKs = 0.84
dAKs = 0.1
cluster = Table.read(folder_clus + f'clus_{mass}M_Age{Age}Myr_imf{imf}_{dist}pc_AKs{AKs}_dAKs{dAKs}.txt', format='ascii')
cat = Table.read(folder_cat + f'pm_{cat_zone}_160w_posEpo1.txt', format='ascii')

n_iter = 1000
n_cat = len(cat)
n_cl = len(cluster)

plot_flag = 0        # set to 1 to render the diagnostic panels for one iteration

# %% =================================================================
# One realisation: scramble + inject
# ====================================================================
def run_injection(sigma_pm, sigma_pos, rng):
    # STEP 1: Scramble the real catalog
    perm_pos = rng.permutation(n_cat)
    perm_pm  = rng.permutation(n_cat)

    cat_s = cat.copy()
    for c in ['x', 'y', 'ra', 'dec', 'xp', 'yp']:
        cat_s[c] = cat[c][perm_pos]
    for c in ['pm_x', 'pm_y', 'dpm_x', 'dpm_y']:
        cat_s[c] = cat[c][perm_pm]

    # STEP 2: Anchors drawn from the catalog itself
    anchor = rng.integers(n_cat)
    ra0, dec0 = cat['ra'][anchor], cat['dec'][anchor]
    pm_anchor = rng.integers(n_cat)
    pmx0, pmy0 = cat['pm_x'][pm_anchor], cat['pm_y'][pm_anchor]

    # STEP 3: Positions and proper motions (direct Gaussian draws)
    ra_cl  = ra0  + (rng.normal(0, sigma_pos, n_cl) / 3600.0) / np.cos(np.radians(dec0))
    dec_cl = dec0 +  rng.normal(0, sigma_pos, n_cl) / 3600.0

    pmx_cl = pmx0 + rng.normal(0, sigma_pm, n_cl)
    pmy_cl = pmy0 + rng.normal(0, sigma_pm, n_cl)

    # STEP 4: Detector coordinates via affine transform from the catalog
    A = np.column_stack([cat['ra'], cat['dec'], np.ones(n_cat)])
    cx, *_ = np.linalg.lstsq(A, cat['x'], rcond=None)
    cy, *_ = np.linalg.lstsq(A, cat['y'], rcond=None)
    A_cl = np.column_stack([ra_cl, dec_cl, np.ones(n_cl)])
    x_cl = A_cl @ cx
    y_cl = A_cl @ cy

    # STEP 5: Fill remaining columns from random donor rows
    donors = rng.integers(0, n_cat, n_cl)

    fake = Table()
    fake['id']   = np.max(cat['id']) + 1 + np.arange(n_cl)
    fake['ra']   = ra_cl*u.deg
    fake['dec']  = dec_cl*u.deg
    fake['x']    = x_cl
    fake['y']    = y_cl
    fake['pm_x'] = pmx_cl
    fake['pm_y'] = pmy_cl
    fake['H'] = cluster['m_hst_f160w']
    fake['J'] = cluster['m_hst_f110w']
    for c in ['f', 'sf', 'dH', 'dJ', 'sxy', 'dpm_x', 'dpm_y',
              'xp', 'yp', 'sx', 'sy']:
        fake[c] = cat[c][donors]
    fake = fake[list(cat.colnames)]

    # STEP 6: Assemble with the ground-truth flag
    out = vstack([cat_s, fake])
    out['is_fake'] = np.zeros(len(out), dtype=int)
    out['is_fake'][len(cat_s):] = 1
    return out, cat_s, fake

# %% =================================================================
# Main loop
# ====================================================================
results = Table(names=('iter', 'seed', 'sigma_pm', 'sigma_pos',
                       'n_inj', 'TP', 'FN', 'FP',
                       'recovery_pct', 'precision_pct'),
                dtype=(int, int, float, float, int, int, int, int,
                       float, float))

for i in range(n_iter):
    seed = 10000 + i
    rng = np.random.default_rng(seed)

    # Random dispersions within the requested ranges
    sigma_pm  = rng.uniform(0.2, 1)      # mas/yr
    sigma_pos = rng.uniform(2.0, 12)     # arcsec

    out, cat_s, fake = run_injection(sigma_pm, sigma_pos, rng)

    print('iter {0}: injected {1} fake members among {2} entries'.format(
        i, n_cl, len(out)))

    # Optional diagnostic plot for a chosen realisation
    if plot_flag and i == 0:
        cmap_color = 'grey'
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        ax[1].hexbin(cat_s['ra'], cat_s['dec'], norm=LogNorm(), cmap=cmap_color)
        ax[0].hexbin(cat_s['pm_x'], cat_s['pm_y'], norm=LogNorm(), cmap=cmap_color)
        ax[2].hexbin(cat_s['J'] - cat_s['H'], cat_s['H'], norm=LogNorm(), cmap=cmap_color)
        ax[1].scatter(fake['ra'], fake['dec'], color='lime', edgecolor='k')
        ax[0].scatter(fake['pm_x'], fake['pm_y'], color='lime', edgecolor='k')
        ax[2].scatter(fake['J'] - fake['H'], fake['H'], color='lime', edgecolor='k')
        ax[2].invert_yaxis()
        ax[1].invert_xaxis()

    # Run the finder
    clus_dic = cluster_finder.finder(out, 'pm_x', 'pm_y',
                                     'ra', 'dec',
                                     'ra', 'dec',
                                     'pm_xy_color', 'J', 'H',
                                     15, 'kernnel', 'minimun',plotting = None,
                                     save_reg=None)

    # Score against the best-matching candidate cluster
    truth = out['is_fake'] == 1
    if clus_dic is None or len(clus_dic) == 0:
        # Finder returned nothing: complete non-detection.
        # Record it explicitly so the statistics stay complete.
        tp, fn, fp = 0, n_cl, 0
        rec  = 0.0
        prec = 0.0
        results.add_row((i, seed, sigma_pm, sigma_pos, n_cl,
                         tp, fn, fp, rec, prec))
        print('  sigma_pm={0:.2f} mas/yr, sigma_pos={1:.1f}" '
              '-> NO CLUSTER FOUND (recovery 0%)'.format(
                  sigma_pm, sigma_pos))
        continue
    best = (0, 0, 0, 0)
    if len(clus_dic) == 0:
        continue
    for cnum in range(len(clus_dic)):
        pred = np.isin(out['id'], clus_dic[f'clus_{cnum}']['id'])
        tn, fp, fn, tp = confusion_matrix(truth, pred).ravel()
        if tp > best[0]:
            best = (tp, fn, fp, tn)

    tp, fn, fp, tn = best
    rec  = 100.0 * tp / (tp + fn) if (tp + fn) > 0 else 0.0
    prec = 100.0 * tp / (tp + fp) if (tp + fp) > 0 else 0.0

    results.add_row((i, seed, sigma_pm, sigma_pos, tp + fn, tp, fn,
                     fp, rec, prec))
    print('  sigma_pm={0:.2f} mas/yr, sigma_pos={1:.1f}" '
          '-> TP={2}, FN={3}, FP={4}, recovery={5:.1f}%'.format(
              sigma_pm, sigma_pos, tp, fn, fp, rec))
    if i%100 == 0:
        print(30*'_' + f'\nIteretaion {i}')
        print(30*'_')
results.write('injection_results.fits', overwrite=True)

# %% =================================================================
# Summary by dispersion bin
# ====================================================================
pm_edges  = np.array([0.2, 0.4, 0.6, 1.0])
pos_edges = np.array([2, 4, 8, 12])

tex_name = pruebas + f'dispersion_{cat_zone}_{mass}M_Age{Age}Myr_imf{imf}_{dist}pc_AKs{AKs}_dAKs{dAKs}.tex'
col_names = ['pm_bin', 'pos_bin', 'N_ite',
             'mean_TP', 'mean_FN', 'mean_FP',
             'mean_rec', 'sig_rec']
float_cols = col_names[3:]

summary = Table(names=col_names,
                dtype=(str, str, int, float, float, float, float, float))

for i in range(len(pm_edges) - 1):
    for j in range(len(pos_edges) - 1):
        sel = ((results['sigma_pm']  >= pm_edges[i])  & (results['sigma_pm']  < pm_edges[i+1]) &
               (results['sigma_pos'] >= pos_edges[j]) & (results['sigma_pos'] < pos_edges[j+1]))
        if np.sum(sel) == 0:
            continue
        sub = results[sel]
        summary.add_row((
            '{0:.1f}-{1:.1f}'.format(pm_edges[i], pm_edges[i+1]),
            '{0:.1f}-{1:.1f}'.format(pos_edges[j], pos_edges[j+1]),
            len(sub),
            np.mean(sub['TP']), np.mean(sub['FN']), np.mean(sub['FP']),
            np.mean(sub['recovery_pct']), np.std(sub['recovery_pct'])))

# Formatted copy: all floats as strings with 2 decimals
summary_fmt = summary.copy()
for c in float_cols:
    summary_fmt[c] = ['{0:.2f}'.format(v) for v in summary[c]]

summary_fmt.write(tex_name, format='ascii.latex', overwrite=True)

# Post-process: starred environment and LaTeX-styled header
header_tex = (r'[pm] & [pos] & $N_{\mathrm{ite}}$ & '
              r'$\overline{\mathrm{TP}}$ & $\overline{\mathrm{FN}}$ & '
              r'$\overline{\mathrm{FP}}$ & $\overline{\mathrm{rec}}$ & '
              r'$\sigma_{\mathrm{rec}}$ \\')

with open(tex_name, 'r') as f:
    lines = f.read().split('\n')

for k, line in enumerate(lines):
    if 'pm_bin' in line:
        lines[k] = header_tex
        break

tex = '\n'.join(lines)
tex = tex.replace(r'\begin{table}', r'\begin{table*}')
tex = tex.replace(r'\end{table}', r'\end{table*}')

with open(tex_name, 'w') as f:
    f.write(tex)

# Terminal output with plain, readable names
print('\nSummary by dispersion bin:')
summary_fmt.pprint_all()