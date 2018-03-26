# README #
mstm-scan

Chad Heaps  
Northwestern University  
March 2018  

This program is an extension of Daniel Mackowski's MSTM software, which may be found at http://eng.auburn.edu/users/dmckwski/scatcodes/. If the software is used in a publication, please cite the original paper for MSTM:

Mackowski and Mishchenko, "A multiple sphere T-matrix Fortran code for use on parallel computer clusters," Journal of Quantitative Spectroscopy and Radiative Transfer, 112: 2182-2192, 2011

The code has been repeatedly modified to fit the needs of a variety of research projects in the Schatz Group.  As a result, particular functions from the original program may not function correctly and no guarantee is given to performance or accuracy from mstm-scan. 

**Significant additions**:  
1.  The option of a dipole source is offered in addition to plane wave and gaussian beam.
2.  The ability to easily scan multiple including the use of wavelength dependent refractive index files for plasmonic materials like gold and silver.
3.  An option to write the electric field at a near-field point with a wavelength scan
4.  A few extra input options including optionally writing the scattering data and selecting a particular wave length for a nearfield calculation

**Possibly altered features from original code**.
1.  The parallelization with MPI.  I'm pretty comfortable this is fine.
2.  Random orientation calculations.  I only use fixed-orientation calculations, so all of the wavelength scanning and dipole source modifications have only been done in that part of the program.
3.  Layered spheres:  I did my best to maintain the same structure as the original program in terms of refractive indices and checking if the host sphere is the medium etc.  Still, this is entirely untested.  Also, the dipole source needs to be outside of all nanoparticles.
4.  Optical activity - In most places the l/r options are filled using the same value.  The code should still run, but all of the dipole code is intended for materials with equal l/r refractive indices


**Modifications to input file**
1.  The nearfield\_output\_data option now dictates where the field intensity or components are printed for the near-field point scan
2.  sphere\_component is a list of N integers for N particles identifying the material of the particle.  This is for heterogenous clusters like one silver, one gold particle
3.  number\_sphere\_components - Just gives the number of different materials.  Added to make reading the input easier
4.  refractive\_index\_file  - A list of filenames with one on each line for each material.  Should be a total of number\_sphere\_components files listed
5.  calculation\_wavelengths - Specify the start, stop, npoints for a wavelength scan.  Endpoints included
6.  dipole\_coordinates - The cartesian coordinates, in microns, giving the location of the dipole source
7.  dipole\_calculation\_type - 0: plane wave (or Gaussian beam) 3: Use a dipole source and evaluate the near field at a single point.  Options 1 and 2 are for Raman scattering but did not make the cut in the last update.
8.  efield\_file - Used to be an actual file.  Now just specifies the orientation of the point dipole source as 'x', 'y', or 'z' for the three cartesian directions.
9.  acceptor\_location - The cartesian coordinates, in microns, where the local electric field will be evaluated.  Should be external to all spheres.
10.  acceptor\_out\_file - The output file for the local electric field at acceptor\_location.  Note: What is printed in the file is determiend by nearfield\_output\_data and the file is only written when dipole\_calculation\_type = 3.
11. nf\_calc\_wavelength - The wavelength where the scattering coefficients will be saved and the near field grid calculation will be executed if the option is selected.  It will pick the closest sampled wavelength if the precise value isn't in the list
12.  write\_scat\_data - A boolean (0,1) for writing scattering, exctinction, etc.  I reformatted the printing of the scattering data for easier plotting


**General notes**
I try to keep subroutines names consistent with the original code.  There is an additional 'getscandata' to supplement 'getrunparameters' providing an interface to the input file parameters.  spheredipolecoef is analagous to sphereplanewavecoef, etc.Whever possible, I note significant changes with my initials, CWH, and a comment on the change, if you are wondering if something has been modified.

The Makefile includes basic functionality for compilation.  Just change the compiler to match your machine.

There are a few example input files for different calculation parameters.

The scattering\_coefficient\_file should be suitable for near-field plots with the Jupyter Notebook

The dipole is added using the Green's dyadic formalism:  Relevant publications include:  

Kerker, M.; Wang, D.-S. & Chew, H. Appl. Opt., 19, 4159 (1980). OSA,  Surface enhanced Raman scattering (SERS) by molecules adsorbed at spherical particles: errata.  
Ausman, L. K. & Schatz, G. C. J. Chem. Phys., 129, 054704 (2008). Whispering-gallery mode resonators: Surface enhanced Raman scattering without plasmons.  
Heaps, C. W. & Schatz, G. C. J. Chem. Phys., 146, 224201 (2017). Modeling super-resolution SERS using a T-matrix method to elucidate molecule-nanoparticle coupling and the origins of localization errors.  

