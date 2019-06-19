WFIRST coronagraph Phase B PROPER model for IDL
John Krist
28 May 2019

This directory contains IDL routines for the IDL version of the Phase B
coronagraph model for PROPER.  

-------------
Prerequisites
-------------

1) IDL Astronomy User's Library, available from:

	idlastro.gsfc.nasa.gov

2) PROPER version 3.1 or higher for IDL, available from:

 	proper-library.sourceforge.net

   Install the PROPER library as directed in the PROPER documentation.

3) Intel Math Kernel Library; a free version is available from:  

	software.intel.com/en-us/performance-libraries

4) WFIRST CGI Phase B optical errors and masks package, available from:

	wfirst.ipac.caltech.edu

   These files should be installed on a local disk to avoid slow downs from
   repeatedly accessing them over a network.

----- 
Setup
----- 

1) Switch to the directory containing the Phase B IDL code.

2) Compile the interface to the Intel MKL FFT. This requires a C compiler is
   installed. From within IDL, issue the command:

	.run compile_ffti

3) Set the "data_dir" variable in both wfirst_phaseb_compact.pro and
   wfirst_phaseb.pro to point to the directory containing the Phase B
   optical errors and masks (that directory should have subdirectories
   such as hlc_20190210, maps, pol, etc.).  You can set this variable
   when calling the prescription as well, but you might as well give it
   a valid default value.

--------------
Using the code
--------------

All of the files in the Phase B IDL code directory should be copied to your
working directory, including the compiled library for the Intel FFT interface.
You can then run the Phase B PROPER prescription using prop_run or
prop_run_multi.

See the WFIRST Phase B prescription documentation for further details,
including the available parameters.

--------
Examples
--------

The examples directory contains demonstrations of how to use the code. Each can
be executed from within idl with a .run command.  These include:

run_flatten.pro
  Generate aberrated HLC images without and with wavefront flattening.
 
run_hlc_input_fields.pro
  Generate aberrated HLC images with full model, then generate aberrated
  exit pupil fields with full model and use them as inputs to the compact
  model.

run_hlc.pro
  Generate unaberrated Dwight's HLC using default DM wavefront patterns,
  then with DM pistons instead, then the aberrated model with DM pistons.

run_hlc_erkin.pro
  Generate unaberrated Erkin's HLC using default DM wavefront patterns.

run_spc_ifs.pro
  Generate unaberrated SPC-IFS model (long bandpass).

run_spc_wide.pro
  Generate unaberrated SPC wide-FOV model.

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

