# HST astrophotometry pipeline

The scripts are located in `/Users/amartinez/Desktop/Projects/SOMA_HST_pm/scripts/`

## Photometry

First, we run *Starfinder* over the "drizzeled" images. (Drizzeled images have been corrected from distortions, cosmic rays etc)
> The version with the `_noAlig.pro` subfix are for the HST drizzeld imagaes that *ARE NOT ALIGN WITH GAIA*. The works in the same way but they produces star lists with different names.

1. `extractpsf_GUI_noAlig.pro` 
    Extract the psf.

 
        idl -e "extractpsf_GUI_noAlig, '<zone>', '<band>', '<epoch>'"

    > That is an apadted version of *Starfinder* widget, that uses ``PSF_EXTRACT``, that gives higher SNR PSFs than ``PSFMAKER``somehow.  If you want to use the widget type ``idl 87`` and then ``> xtstarfinder``. Find the *Starfinder widget manual* in `/Users/amartinez/Desktop/PhD/StarFinder/starfinder_manual.pdf`

2. `astrophot_noAlig.pro`
   
    Generate the stars lists.

         idl -e "astrophot_noAlig, '<zone>', '<band>', '<epoch>'"


## Astrometry

We align the stars list witn the *Gaia* stars. By default we a use a degree 2 polynomial.

4. `hst_gaia_alignment.py`

   >  `hst_photometry.py` I intruced the ZP calculation and calibration inside the alignment script.
   > For the ZP calculation a function is called ``from hst_irZP import get_vegazp``. This escripts is in */Users/amartinez/Desktop/pythons_imports/hst_irZP.py*

5. `hst_relative_alignment.py`. This calculate the relative proper motions (using one of the epochs as reference frame). Also calculates the ZP with ``from hst_irZP import get_vegazp``
   
