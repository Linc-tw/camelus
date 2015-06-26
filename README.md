# Camelus
Counts of Amplified Mass Elevations from Lensing with Ultrafast Simulation  
Chieh-An Lin (CEA Saclay)  
Release v1.2 - 2015-04-06 
<p align="center"><a href="http://species.wikimedia.org/wiki/Camelus"><img src="http://www.cosmostat.org/wp-content/uploads/2014/11/Logo_Camelus_fig_name_vertical.png" width="240px" /></a></p>

## Description

Camelus is a fast weak-lensing peak-count modeling algorithm in C. It provides a prediction on peak counts from input cosmological parameters.

Here is the summary of the algorithm:
  - Sample halos from a mass function
  - assign density profiles, randomize their positions
  - Compute the projected mass, add noise
  - Make maps and create peak catalogues

For a more detailed description, please take a look at Lin & Kilbinger (2015a).

## Requirements

The following softwares are required:
  - [cmake](http://cmake.org/cmake/resources/software.html)
  - gcc
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

Current release: Camelus v1.2

##### New features in v1.2 - Apr 06, 2015:
  - Improved the computation speed by a factor of 6~7
  - Converted the halo array structure into a binned structure, called "halo_map"
  - Converted the galaxy tree structure into a binned structure, called "gal_map"
  - Added the population Monte Carlo approximate Bayesian computation (PMC ABC) algorithm

##### New features in v1.1 - Jan 19, 2015:
  - Fixed the bug from calculating halo radii

##### New features in v1.0 - Oct 24, 2014:
  - Fast weak lensing peak count modeling

## References

  - [Bartelmann & Schneider (2001)](http://arxiv.org/abs/astro-ph/9912508). Phys. Rep., 340, 291.
  - [Fan et al. (2010)](http://arxiv.org/abs/1006.5121). ApJ, 719, 1408.
  - [Lin & Kilbinger (2015a)](http://arxiv.org/abs/1410.6955). A&A, 576, A24.
  - [Lin & Kilbinger (2015b)](http://arxiv.org/abs/1506.01076). Submitted to A&A.
  - [Marin et al. (2011)](http://arxiv.org/abs/1101.0955).
  - [Takada & Jain (2003a)](http://arxiv.org/abs/astro-ph/0209167). MNRAS, 340, 580.
  - [Weyant et al. (2013)](http://arxiv.org/abs/1206.2563). ApJ, 764, 116.

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
This gives the values of the mass function from Jenkins et al. (2001) for M from 10^9 to 10^17 [M_sol h^2/Mpc^3].

`$ ./camelus 2`  
This yields a halo list from sampling using customized parameters.

`$ ./camelus 3`  
This creates a peak catalogue and its histogram of S/N values.   
A halo list, a galaxy list, a smoothed noise map, a smoothed noisy map will be made as well.

`$ ./camelus 4 N`  
This creates N independent peak lists with the same cosmology and settings from `.par` files.   
Using this command can reduce the computation time, since the cosmological computations are factorized.

`$ ./camelus 4 N Omega_m sigma_8`  
Same as above, but the inputs Omega_m and sigma_8 will overwrite the values from `.par` files, and creates N independent peak lists.

`$ ./camelus 5`  
ABC computation which requires an observation data that we have provided an example in `demo`.  
Parameters are defined in `peakParam.par`. Only several summary statistics are available.  
This gives posterior samples of Omega_m-sigma_8 constraints. Only this combination is available. 
Please check Lin & Kilbinger (2015b) for more details.

