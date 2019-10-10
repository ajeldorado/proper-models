WFIRST coronagraph Phase B PROPER model for Python 
John Krist (JPL)
7 October 2019

This directory contains routines for the Python version of the Phase B
coronagraph model for PROPER.  

-------------
Prerequisites
-------------

1) PROPER version 3.1 or higher for Python, available from:

 	proper-library.sourceforge.net

   Install the PROPER library as directed in the PROPER documentation.
   It is recommended that you install the FFTW or Intel MKL library, as
   described in the PROPER manual.

2) WFIRST CGI Phase B optical errors and masks package, available from:

	wfirst.ipac.caltech.edu

   These files should be installed on a local disk to avoid slow downs from
   repeatedly accessing them over a network.

----- 
Setup
----- 

1) After unpacking the wfirst_phaseb_proper package, enter the top
   level directory (containing setup.py) and issue the command:

	python setup.py install

   The package will be installed in your local library directory, which
   is usually something like:

     ~/.local/lib/python2.7/site-packages

2) The package needs to be told where the Phase B error maps and mask
   files are stored.  Switch to the directory you unpacked those into
   (that directory should have subdirectories such as hlc_20190210, maps, 
   pol, etc.).  Then, from within Python, enter:

	import wfirst_phaseb_proper
	wfirst_phaseb_proper.set_data_dir()

   This will set the data_dir variable in your library's __init__.py file
   to this directory (this value can be overridden during execution using 
   the data_dir parameter).

--------------
Using the code
--------------

The wfirst_phaseb.py and wfirst_phaseb_compact.py prescriptions must be run
from your local working directory - they cannot be executed from another
directory.  You can copy them from your library to the current directory from
within Python like so:

	import wfirst_phaseb_proper
	wfirst_phaseb_proper.copy_here()  

Though wfirst_phase.py and wfirst_phaseb_compact.py are copied over, you still
need to import the wfirst_phaseb_proper package whenever you use them to have
certain variables defined.

See the WFIRST Phase B prescription documentation for further details,
including the available parameters.

--------
Examples
--------

The examples are provided as demonstrations of how to use the code.  These
must be executed from the local working directory.  To copy them from your
library, get into Python and issue the commands:

	import wfirst_phaseb_proper
	wfirst_phaseb_proper.copy_examples_here()

The available examples are:

run_flatten.py
  Generate aberrated HLC images without and with wavefront flattening.
 
run_hlc_input_fields.py
  Generate aberrated HLC images with full model, then generate aberrated
  exit pupil fields with full model and use them as inputs to the compact
  model.

run_hlc.py
  Generate unaberrated Dwight's HLC using default DM wavefront patterns,
  then with DM pistons instead, then the aberrated model with DM pistons.

run_hlc_erkin.py
  Generate unaberrated Erkin's HLC using default DM wavefront patterns.

run_spc_spec.py
  Generate unaberrated and EFC-corrected SPC-SPEC models (long bandpass).

run_spc_wide.py
  Generate unaberrated and EFC-corrected SPC wide-FOV models.

-----------
DM settings
-----------

In the examples directory are DM settings for selected cases:

errors_polaxis10_dm.fits
  DM1 pistons (meters) that flatten the aberrated wavefront phase at the FPM
  (use with use_errors=1 and polaxis=10, which is the mean X+Y polarization).

hlc_dm1.fits, hlc_dm2.fits
  DM1 & DM2 pistons (meters) corresponding to an unaberrated solution for
  Dwight's HLC.

hlc_with_aberrations_dm1.fits, hlc_with_aberrations_dm2.fits
  DM1 & DM2 pistons (meters) corresponding to an aberrated solution for
  Dwight's HLC, including polarization. Use with use_errors=1 and polaxis=10.

spc-wide_with_aberrations_dm1.fits, spc-wide_with_aberrations_dm2.fits,
spc-spec_long_with_aberrations_dm1.fits, spc-spec_long_with_aberrations_dm2.fits
  DM1 & DM2 pistons (meters) corresponding to an aberration solution for
  the wide-field SPC (band 4) and SPC SPEC (band 3), including polarization.
  Use with use_errors=1 and polaxis=10.

 
