Camelus
=======
Counts of Amplified Mass Elevations from Lensing with Ultrafast Simulation  
Chieh-An Lin (CEA Saclay)  
Release v1.31 - 2016-03-23 
<p align="center"><a href="http://species.wikimedia.org/wiki/Camelus"><img src="http://www.cosmostat.org/wp-content/uploads/2014/11/Logo_Camelus_fig_name_vertical.png" width="240px" /></a></p>

Description
-----------

Camelus is a fast weak-lensing peak-count modeling algorithm in C. It provides a prediction on peak counts from input cosmological parameters.

Here is the summary of the algorithm:
  - Sample halos from a mass function
  - Assign density profiles, randomize their positions
  - Compute the projected mass, add noise
  - Make maps and create peak catalogues

For a more detailed description, please take a look at [Lin & Kilbinger (2015a)](http://arxiv.org/abs/1410.6955).

## Requirements

The following softwares are required:
  - [cmake](http://cmake.org/cmake/resources/software.html)
  - [gcc](http://gcc.gnu.org/)
  - [gsl](http://www.gnu.org/software/gsl/)
  - [fftw](http://www.fftw.org/)
  - [nicaea v2.5](http://www.cosmostat.org/nicaea.html)

## Compilation

For Mac users, do the follows before compilation:
```
$ export CC=gcc
$ export CXX=g++
```
or use `setenv` command in tcsh.

To compile the package:
```
$ export NICAEA={PATH_OF_NICAEA}/nicaea_2.5
$ cd build
$ cmake ..
$ make
```

To get program instructions:
```
$ ./camelus
```

## Updates

Current release: Camelus v1.31

##### Forked version: implementation of shear bias as function of local density - 04/07/2017

##### New features in v1.31 - Mar 22, 2016:
  - Made installation more friendly by removing the dependency on cfitsio and mpi
  - Added the routine for computing 1-halo & 2-halo terms of the convergence profile
  - Flexible parameter space for PMC ABC
  - Remove files: FITSFunctions.c/.h

##### New features in v1.3 - Dec 09, 2015:
  - New files: constraint.c/.h
  - Allowed multiscale peaks in one data vector
  - Allowed a data matrix from several realizations
  - Used the local galaxy density as the noise level in the S/N
  - Increased the parameter dimension for PMC ABC
  - Changed the summary statistic options for PMC ABC

Unavailable features because of the exterior file dependency:
  - New files: FITSFunctions.c/.h
  - Added the mask option and nonlinear filtering

##### New features in v1.2 - Apr 06, 2015:
  - Improved the computation speed by a factor of 6~7
  - Converted the halo array structure into a binned structure, called "halo_map"
  - Converted the galaxy tree structure into a binned structure, called "gal_map"
  - New files: ABC.c/.h
  - Added the population Monte Carlo approximate Bayesian computation (PMC ABC) algorithm

##### New features in v1.1 - Jan 19, 2015:
  - Fixed the bug from calculating halo radii

##### New features in v1.0 - Oct 24, 2014:
  - Fast weak lensing peak count modeling

## References

  - [Baltz et al. (2009)](http://arxiv.org/abs/0705.0682) - JCAP, 1, 15
  - [Bartelmann & Schneider (2001)](http://arxiv.org/abs/astro-ph/9912508) - Phys. Rep., 340, 291
  - [Fan et al. (2010)](http://arxiv.org/abs/1006.5121) - ApJ, 719, 1408
  - [Hetterscheidt et al. (2005)](http://arxiv.org/abs/astro-ph/0504635) - A&A, 442, 43
  - [Lin & Kilbinger (2015a)](http://arxiv.org/abs/1410.6955) - A&A, 576, A24
  - [Lin & Kilbinger (2015b)](http://arxiv.org/abs/1506.01076) - A&A, 583, A70
  - [Lin et al. (2016)](http://arxiv.org/abs/1603.06773) - Submitted to A&A
  - [Marin et al. (2011)](http://arxiv.org/abs/1101.0955)
  - [Oguri & Hamana (2011)](http://arxiv.org/abs/1101.0650) - MNRAS, 414, 1851
  - [Takada & Jain (2003a)](http://arxiv.org/abs/astro-ph/0209167) - MNRAS, 340, 580
  - [Takada & Jain (2003b)](http://arxiv.org/abs/astro-ph/0304034) - MNRAS, 344, 857
  - [Weyant et al. (2013)](http://arxiv.org/abs/1206.2563) - ApJ, 764, 116
  - Wright & Brainerd (2000) - ApJ, 534, 34

## Contact information

Authors:
  - [Chieh-An Lin](http://linc.tw/)
  - [Martin Kilbinger](http://www.cosmostat.org/people/kilbinger/)
  - [François Lanusse](http://www.cosmostat.org/people/flanusse/)

Please feel free to send questions, feedback and bug reports to chieh-an.lin (at) cea.fr.  
Check also the package [web page](http://www.cosmostat.org/software/camelus/).

## Tutorial

Go to param and modify `.par` files to customize parameters.

`$ ./camelus 1 z`  
This computes the mass function for M from 10^9 to 10^17 [M_sol h^2/Mpc^3].

`$ ./camelus 2`  
This creates a halo catalogue from random sampling using customized parameters.

`$ ./camelus 3 z_l M z_s`  
This computes the 1-halo and 2-halo terms of the convergence profile proposed by Baltz et al. (2009), with n = 2.  
z_l for the lens redshift, M for the lens mass, z_s for the source redshift.

`$ ./camelus 4`  
This creates a lensing map using the first FFT or DC filter entry.   
Intermediate results are also given, such as the halo catalogue, the noisy/noiseless galaxy catalogues, the noisy map, and the truth map.

`$ ./camelus 5`  
This creates a peak catalogue and its histogram of S/N values using the first FFT or DC filter entry.   
Intermediate results are also given, such as the halo catalogue, the galaxy catalogue, and the maps.

`$ ./camelus 6`  
This creates a data vector of peak counts from different scales.

`$ ./camelus 7 N`  
This creates N independent realizations of peak-count vector with the same cosmology and settings from `.par` files.

`$ ./camelus 7 N Omega_m sigma_8 w0_de`  
Same as above, but the inputs Omega_m, sigma_8, and w0_de will overwrite the values from `.par` files, and creates N independent realizations.

`$ ./camelus 8`  
This gives an example of approximate Bayesian computation (ABC).  
The routine is only available for default peak pipeline settings.  
The routine requires an observation data that we have provided an example in `demo`.  
The routine gives posterior samples of parameters defined by `ABC_doParam`.  
Please check [Lin & Kilbinger (2015b)](http://arxiv.org/abs/1506.01076) for more details.

