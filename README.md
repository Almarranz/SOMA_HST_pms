# HST astrophotometry pipeline

The scripts are located in `/Users/amartinez/Desktop/Projects/SOMA_HST_pm/scripts/`

## Photometry

First, we run *Starfinder* over the "drizzeled" images. (Drizzeled images have been corrected from distortions, cosmic rays etc)

1. `extractpsf.pro` 
    Extract the psf.

    .. code:: 
        idl -e \"extractpsf, 'G028.20-00.05', '160w', '2'\

    This script generate the psf automatically, and seems to work correctly. However, there is the possiboity to use *Starfinde* widget to do so. Just open a terminal al type:

    .. code:: 

        idl 89
        > xtstarfinder 

    Manual for *Starfinder* widget : `/Users/amartinez/Desktop/PhD/StarFinder/starfinder_manual.pdf`
