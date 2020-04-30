Habex example PROPER prescription

- Created on 22 January 2020 by John Krist, 
  Jet Propulsion Laboratory, California Institute of Technology
- Modified in February 2020 by A.J. Riggs to include MATLAB and Python
  translations and example scripts to run the PROPER model.

The script habex.pro (or .m or .py) uses calls to the PROPER propagation
library to represent the hypothetical Habex telescope and coronagraph (single
polarization channel) evaluated in the NASA study for the 2020 Decadal Survey.
It has been simplified for clarity, among other changes.

(IDL version only) The user-written function "tag_exists" (included) checks if
a parameter is provided in the optional structure "optval" using the PASSVALUE
mechanism of prop_run.  Check the file for available optional parameters.

The user needs to change the map_dir variable to point to the directory
containing the optical surface error maps. The maps are available for download
at a separate website, https://www.astro.caltech.edu/~gruane/maps.zip


Example calls in IDL would be:

  Return the complex-valued field at the image plane with the default
  sampling of 0.2 lam0/D (default lam0 is 0.5 microns):

    lam = 0.5        ;-- wavelength in microns
    gridsize = 1024  ;-- computational gridsize (must be factor of 2)
    prop_run, 'habex', field, lam, gridsize

  Return the field with the star offset in X from the focal plane mask by
  7 lam0/D (lam0 may be defined using the optional lam0 parameter):

    prop_run, 'habex', field, lam, gridsize, passvalue={xoffset:7.0}


For example calls in MATLAB and Python 3, open and run the script titled:
script_simple_habex_model_test


To simulate wavefront correction with this model in MATLAB or Python 3, run the FALCO script titled:
falco_main_Habex_VC


(c) 2020 California Institute of Technology

-- END --   
