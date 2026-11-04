#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 11:18:53 2026

@author: amartinez
"""

import mwdust
from astropy.coordinates import SkyCoord
import sys
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

# folder = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/MainStarsReg/'
# mains = Table.read(folder + 'main_stars.txt', format = 'ascii')

folder = '/Users/amartinez/Desktop/Projects/SOMA_HST_pm/SOMA_HST_pms_variability/MainStarsReg/'
mains = Table.read(folder + 'SOMA_HST_Reid_table.ecsv', format = 'ascii')
star = mains[-1]

print(star['Target', 'd'])
coord = SkyCoord(ra=star['RA'], dec=star['Dec'], unit=(u.hourangle, u.deg)).galactic

# Initialize for the band you want; many filters are supported
dust = mwdust.Combined19(filter='2MASS Ks')

# Cumulative extinction at (coor.l.value, coord.b.value, d):
A = dust(coord.l.value, coord.b.value, star['d']*1000)   # l, b in degrees, d in parsecs

print(A)