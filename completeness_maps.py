#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from astropy.io import fits
from scipy.ndimage import convolve


# In[2]:


def create_circular_mask(radius):
    diameter = 2 * radius + 1
    center = radius
    y, x = np.ogrid[-center:center+1, -center:center+1]
    mask = x**2 + y**2 <= radius**2
    return mask.astype(float)


# In[3]:


def sum_in_circle(arr, radius):
    mask = create_circular_mask(radius)
    return convolve(arr, mask, mode='constant', cval=0.0)


# In[4]:


chip = 'nrca1'
pointing = 6
pixscale = 0.05
radius = 50

# path = f'F212N/{chip}_p{pointing}/Completeness/'
band = '160W'
zone = 'G028.20-00.05'
path = f'/Users/amartinez/Desktop/Projects/SOMA_HST_pm/sf/results/{zone}/F{band}/epoch1/completeness/'

mag_min = 12   #minimum magnitude for which completeness is tested
mag_max = 22   #maximum magnitude for which completeness is tested
mag_step = 1 #magnitude step for which completeness is tested
mags = np.arange(mag_min,mag_max+mag_step,mag_step)
print(mags)


# In[5]:


#radius = 100
#inmap = fits.getdata(path + f'inmap{16:.4f}.fits')
#inmap = inmap[0:1000,0:1000]
#outmap = fits.getdata(path + f'outmap{16:.4f}.fits')
#outmap = outmap[0:1000,0:1000]
#in_summed = sum_in_circle(inmap, radius)
#out_summed = sum_in_circle(outmap, radius)
#in_summed = in_summed.astype('float64')
#out_summed = out_summed.astype('float64')
#comp_map = np.divide(out_summed, in_summed, out=np.zeros_like(in_summed), where=out_summed!=0)
#fits.writeto(f'{path}hhh.fits',comp_map,overwrite=True)


# In[6]:


for thismag in mags:
    inmap = fits.getdata(path + f'inmap{thismag:.0f}.fits')
    outmap = fits.getdata(path + f'outmap{thismag:.0f}.fits')
    print(f'Read in maps for magnitude {thismag:.0f}.')
    in_summed = sum_in_circle(inmap, radius)
    out_summed = sum_in_circle(outmap, radius)
    in_summed = in_summed.astype('float64')
    out_summed = out_summed.astype('float64')
    comp_map = np.divide(out_summed, in_summed, out=np.zeros_like(in_summed), where=out_summed!=0)
    fits.writeto(f'{path}completeness{thismag:.0f}.fits',comp_map,overwrite=True)
    print(f'Processed magnitude {thismag:.0f}.')


# In[7]:


completeness_mean = []

for thismag in mags:  
    inmap = fits.getdata(path + f'inmap{thismag:.0f}.fits')
    outmap = fits.getdata(path + f'outmap{thismag:.0f}.fits')
    print(f'Read in maps for magnitude {thismag:.0f}.')

    in_summed = sum_in_circle(inmap, radius).astype('float64')
    out_summed = sum_in_circle(outmap, radius).astype('float64')
    # --- completeness promedio global ---
    total_in = np.sum(in_summed)
    total_out = np.sum(out_summed)

    comp_global = total_out / total_in if total_in > 0 else 0
    completeness_mean.append(comp_global)

    print(f'Processed magnitude {thismag:.0f}. Completeness = {comp_global:.3f}')


# In[8]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

fig, ax = plt.subplots(1,1)
ax.set_title(f'{zone}')
ax.plot(mags, np.array(completeness_mean)*100, marker='o')
ax.set_xlabel(f'F{band}')
ax.set_ylabel('Completeness [%]')
# ax.set_ylim(0.4, 1.05)
ax.set_ylim(4, 110)
ax.set_yticks(np.arange(0, 110, 10))
ax.set_xticks(np.arange(11, 23, 1))
ax.grid(ls ='-.')

#plt.gca().invert_xaxis()  # importante en astronomía
#plt.grid()
# plt.show()


# In[ ]:




