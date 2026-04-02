# HST astrophotometry pipeline

The scripts are located in `/Users/amartinez/Desktop/Projects/SOMA_HST_pm/scripts/`

## Photometry

First, we run *Starfinder* over the "drizzeled" images. (Drizzeled images have been corrected from distortions, cosmic rays etc)

1. `extractpsf.pro` 
    Extract the psf.

 
        idl -e \"extractpsf, '<zone>', '<band>', '<epoch>'\

    > This script generate the psf automatically, and seems to work correctly. However, there is the possiboity to use *Starfinde* widget to do so. Just open a terminal al type:

    

            idl 87
            > xtstarfinder 

    Manual for *Starfinder* widget : `/Users/amartinez/Desktop/PhD/StarFinder/starfinder_manual.pdf`

2. `astrophot.pro`
   
    Generate the stars lists.

         idl -e \"astrophot, '<zone>', '<band>', '<epoch>'\


## Astrometry

We align the stars list witn the *Gaia* stars. By default we a use a degree 2 polynomial.

4. `hst_gaia_alignment.py`

   >  `hst_photometry.py` I intruced the ZP calculation and calibration inside the alignment script.
   > For the ZP calculation a function is called ``from hst_irZP import get_vegazp``. This escripts is in */Users/amartinez/Desktop/pythons_imports/hst_irZP.py*

5. `hst_relative_alignment.py`. This calculate the relative proper motions (using one of the epochs as reference frame). Also calculates the ZP with ``from hst_irZP import get_vegazp``
   
