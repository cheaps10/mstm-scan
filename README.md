# README #
mstm-scan

Chad Heaps
Northwestern University
April 2017

An extension of Daniel Mackowski's MSTM code to scan wavelengths, etc. 

Unlike my dipole code, this is based on MSTM v3.0


The makefile will compile the code.
The run\_mstm.sh is just a lazy way to (possibly compile) and run the same input while debugging
You will likely need to change the compiler.  I'm on my mac where I need to use gfortran-mp-6


There are 3 example inputs
mstm.inp is basically the same as one Mackowski sends out.  I just added his refractive index stuff to a file refMSTM.dat
test1.inp  A single nanoparticle donor-acceptor calculation.  Gets it's dipole orientation from efmonomer1
test2.inp  A dimer of 60 nm gold nanoparticles 14 nm apart, roughly following the geometry from http://pubs.acs.org/doi/abs/10.1021/acsphotonics.6b00148  Gets its dipole orientation from efdimer1


Edits from Original MSTM and notes on running

MSTM uses the same trick as SCSMFO where a negative error tolerance for the Mie coefficients will enforce a VSH order of -etol.  Therefore, you can specify the basis set size with mie\_epsilon.

For some reason qabs is NaN as well as the unpolarized qabstot and qscattot.  I am not sure why and haven't looked into that yet.

I changed the code so that the scattering\_coefficient\_file saves all of the wavelengths.  For large basis sets and lots of wavelengths it gets very long.  Not memory overload long but annoyingly long.
Note: It looks like Mackowski added an option to append the amn file, I still may bypass that.

INPUT FILE

Mackowski uses a set of parameter ids to identify and read in parameters.  It makes it very easy to add parameters.  If looking at the source code modules.f90, the module spheredata is responsible for all of the i/o.  In general, parameters are read in using inputdata and stored as private variables.  Usable variables are fetched using the getspheredata, where all of the arguments are optional and it allows you to retreive anything you need from sphere data with getspheredata(privatename=publicname) in your subroutine.

I defined a second routine getdipoledata that can retreive input parameters associated with the dipole calculation

I added a set of input options to the end of the input files (after sphere positions) that include

refractive\_index\_file: This is somewhat self-explanatory.  The Johnson + Christy data is currently in the directory. The souce code is hard coded to read past 4 lines of stuff, then the number of entries then the entries.  Not very flexible...

medium\_refractive\_index: You could accomplish similar behavior using Mackowski's real\ref\_index\_scale\_factor, but that is a pain in the butt.  Caution, however, I am pretty sure Mackowski hard-coded a refractive index of 1 into some function calls that technically require n\_medium

calculation\_wavelengths: wav\_min wav\_max n\_wav  Also self-explanatory.  Creates an evenly spaced list of wavelengths in microns.  Interpolation of refractive index is done using the same DDA routine from the old code.

dipole\_coordinates:  in microns, the cartesian coordinates of the dipole.  Will be scaled by length\_scale\_factor

dipole\_calculation\_type: integer 0,1,2,3  Still not sure what I'm doing here.  0 should be no dipole anything.  3 is Wendu's donor-acceptor calculation.  1 will presumably be plane-wave calculate |E|^4 and 2 widll be dipole scattering for imaging
 
efield\_file:  The file where the code expects to find the cartesian components of thte electric field inducing the dipole.  For Wendu's calculations just the cartesian axis of the donor dipole.  Also where I may write the E-field when using calc type 1.  I need to figure that out.  The code will automatically use the last entry in efield\_file for the rest of the wavelengths if the number of entries is less than n\_wav.  So you can just specify your fixed dipole once and loop over wavelengths

acceptor\_location:  The point to evaluate the nearfield for the donor-acceptor calculation.  Also given in microns.

acceptor\_out\_file:  Where to write the electric field intensity at the acceptor location.  Writes |e|^2 right now...can modify in main.f90

nf\_calc\_wavelength: real (in microns)  After generating the evenly spaced wavelengths for calculation, the code will pick the wavelength in the list closest to the value given here to evaluate the near-field

 
TODO:

1.  Add near-field dipole evaluation
2.  Make sure the medium refractive index functions properly
3.  Figure out the phase of the dipole coefficients and scattered field.  Logan, not surprisingly, makes things pretty unclear, so the same incident coefficients give the same scattered coefficients in the two codes, but not the same scattered field.






mod:	
	rm *.mod; $(MAKE) all

