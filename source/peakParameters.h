

  /*******************************
   **  peakParameters.h		**
   **  Chieh-An Lin		**
   **  Version 2016.03.20	**
   *******************************/


#ifndef __peakParam__
#define __peakParam__

#define _GNU_SOURCE
#define __releaseMenu__

#include <nicaea/errorlist.h>
#include <nicaea/io.h>
#include <nicaea/cosmo.h>
#include <nicaea/nofz.h>
#include <nicaea/halomodel.h>


#include "commonHeader.h"


//-- Computational constants
#define NUMBER_HALOS_MAX        5000000
#define NUMBER_GALAXIES_MAX     5000000
#define CUTOFF_FACTOR_HALO      3.0
#define CUTOFF_FACTOR_HALO_SQ   9.0
#define CUTOFF_FACTOR_GAUSSIAN  2.2
#define CUTOFF_FACTOR_M_AP_TANH 1.2
#define NORM1_STARLET           0.97873802046145453 //-- W(x, y) = psi(x, y)
#define NORM2_STARLET           0.65166118860701647
#define NORM1_M_AP_TANH         5.4323267630348244  //-- W(x) = tanh(x / 0.1) / [x (1 + exp(5 - 150 x) + exp(-47 + 50 x))]
#define NORM2_M_AP_TANH         3.8414742595565707
#define FILLING_THRESHOLD_RATIO 0.5

//-- Error type definition
#define peak_base     -200
#define peak_null     -1 + peak_base
#define peak_badValue -2 + peak_base
#define peak_overflow -3 + peak_base
#define peak_unknown  -4 + peak_base
#define peak_geometry -5 + peak_base
#define peak_mapType  -6 + peak_base
#define peak_match    -7 + peak_base
#define peak_EOF      -8 + peak_base
#define peak_setting  -9 + peak_base


//-- Peak parameters
typedef struct {
  //----------------------------------------------------------------------
  //-- Customized part
  
  //-- Halos
  double z_halo_max;                   //-- Maxmimum redshift for halos
  int N_z_halo;                        //-- Number of sampling redshift bins
  double M_min;                        //-- [M_sol/h] Minimum sampling mass
  double M_max;                        //-- [M_sol/h] Maximum sampling mass
  
  //-- Galaxies
  double z_s;                          //-- If z_s > 0, all source galaxies fixed at z_s
                                       //-- If z_s < 0, use nofz_hm.dat as redshift distribution law parameters
  int doRandGalPos;                    //-- 0 = regular, 1 = random
                                       //-- If z_s < 0, doRandGalPos is automatically set to 1 
  double n_gal;                        //-- [arcmin^-2] Galaxy number density
  double n_gal_obs;                        //-- [arcmin^-2] Galaxy number density
  double sigma_eps;                    //-- [-] Ellipticity dispersion, sigma_eps^2 = <epsilon_1^2> + <epsilon_2^2> = 2 sigma_kappa^2
  int doKappa;                         //-- 0 = gamma, 1 = kappa, 2 = g, 3 = g with linear KS
  int doMask;                          //-- 0 = without, 1 = CFHTLenS W1
  char maskPath[STRING_LENGTH_MAX];    //-- Path of the mask
  
  //-- Field & map
  field_t field;                       //-- Geometry (only 'rectangle' is available in v1.2)
  double Omega[2];                     //-- [arcmin] Field size (theta_x, theta_y)
  double theta_pix;                    //-- [arcmin] Pixel size
  
  //-- Filter
  int doSmoothing;                     //-- 0 = without, 1 = binning and FFT, 2 = direct convolution, 4 = nonlinear filtering,
                                       //-- Sum for performing simultaneously multiple techniques, e.g. 5 = FFT + nonlinear
  int FFT_nbFilters;                   //-- Number of linear filters for FFT
  filter_t *FFT_filter;                //-- Filter type for FFT
  double *FFT_scale;                   //-- [arcmin] Filter size for FFT
  int DC_nbFilters;                    //-- Number of linear filters for direct convolution
  filter_t *DC_filter;                 //-- Filter type for direct convolution
  double *DC_scale;                    //-- [arcmin] Filter size for direct convolution
  int MRLens_nbScales;                 //-- Number of scales for MRLens
  double MRLens_FDR;                   //-- False detection rate for MRLens
  
  //-- Peak historgram
  int N_nu;                            //-- [-] Number of S/N bins
  double *nu_bin;                      //-- S/N bins
  int N_kappa;                         //-- [-] Number of kappa bins
  double *kappa_bin;                   //-- kappa bins
  
  //-- ABC
  char ABC_obsPath[STRING_LENGTH_MAX]; //-- Path of the observation data
  int ABC_f;                           //-- Dimension of parameter set
  int *ABC_doParam;                    //-- Parameters to include into constraints
                                       //-- 0 = Omega_m, 1 = Omega_de, 2 = Omega_b, 3 = n_s, 4 = h_100,
                                       //-- 5 = sigma_8, 6 = w0_de,    7 = w1_de,   8 = c_0, 9 = beta_NFW
  int ABC_Q;                           //-- Number of particles
  double ABC_r_stop;                   //-- Shutoff accept rate
  char ABC_summ[STRING_LENGTH_MAX];    //-- Summary type
  
  //-- Label
  char simulName[STRING_LENGTH_MAX];   //-- Name of the parameter set
  
  //---------------------------------------------------------------------
  //-- Default part
  
  double dlogM;                        //-- [-] Bin width for mass function
  double dz_gal;                       //-- [-] Galaxy redshift binwidth
  double theta_CCD_inv[2];             //-- [arcmin^-1] Inverse of CCD mask pixel size, only used when doMask = 1
  char tempPath[STRING_LENGTH_MAX];    //-- Path of the temporary directory
  
  //----------------------------------------------------------------------
  //-- Precomputed part
  
  double area;                         //-- [arcmin^2] Field area
  double w_s;                          //-- [Mpc/h] Comoving radial distance of the source
  double D_s;                          //-- [Mpc/h] Angular diameter distance of the source
  double dz_halo;                      //-- [-] Halo redshift binwidth
  int N_M;                             //-- Number of mass bins
  int N_z_gal;                         //-- Number of galaxy redshift bins
  double sigma_half;                   //-- [-] sigma_eps / sqrt(2)
  double theta_pix_inv;                //-- [arcmin^-1] Inverse of pixel size
  double sigma_pix;                    //-- [-] Noise dispersion for each pixel
  gsl_rng *generator;                  //-- Random number generator
  
  //-- Filter-related
  int doNonlinear;                     //-- 0 = No, 1 = Yes, derived from doSmoothing
  int nbFilters;                       //-- Total number of filters
  int smootherSize;                    //-- Length of the smoother
  filter_t *filter;                    //-- Filter type
  double *FFT_sigma_noise;             //-- [-] Expected noise level for linear filters, WARNING used in paperIII, but not used in the released version
  double *DC_sigma_noise;              //-- [-] Expected noise level for linear filters, WARNING not used
  double *DC_scale_invSq;              //-- [arcmin^-2] Inverse squared/non-squared size of linear filters
  double *DC_cut_sq;                   //-- [arcmin^2] Squared/non-squared cutoff size of linear filters
  
  //-- Only used if field = rectangle
  int resol[2];                        //-- [pix] round(Omega / theta_pix)
  double *FFT_scaleInPix;              //-- [pix] Linear filter size in pixel
  double *DC_scaleInPix;               //-- [pix] Linear filter size in pixel
  int bufferSize;                      //-- [pix] 0-padding size
  int FFTSize;                         //-- [pix] FFT table size
  double FFTNormFactor;                //-- [-] FFT normalization factor
  double limit[4];                     //-- [arcmin] Field limits (theta_x_min, theta_x_max, theta_y_min, theta_y_max)
  
  //----------------------------------------------------------------------
  //-- Running part
  
  //-- Mode
  int printMode;                       //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  int doNoise;                         //-- 0 = noiseless, 1 = noisy
  
  //-- Pipeline index
  int cosmoInd;                        //-- Index for cosmology, Paper II and Paper III only
  int scaleInd;                        //-- Index for scale, for multiscale analysis
  int realizationInd;                  //-- Index for realization
  
  //-- MPI
  int MPISize;                         //-- Number of MPI processors
  int MPIInd;                          //-- Index for MPI processors
  
  //----------------------------------------------------------------------
  //-- Paper I only
  
  int cut;                             //-- aardvark only
  int nbPatches, *patchList;           //-- aardvark only
  int patchOrder, RTPatch;             //-- aardvark only
  int fastId, noiseId;                 //-- aardvark only
  
  //----------------------------------------------------------------------
} peak_param;


//-- Functions related to cosmo_hm
cosmo_hm *initialize_cosmo_hm_default(error **err);
void read_cosmo_hm(char name[], cosmo_hm **cmhm, error **err);
cosmo_hm *updateCmhm(cosmo_hm *oldCmhm, double Omega_m, double sigma_8, double w0_de, error **err);

//-- Functions related to peak_param
void free_peak_param(peak_param *peak);
void read_peak_param(char name[], peak_param *peak, error **err);
void set_peak_param(cosmo_hm *cmhm, peak_param *peak, error **err);

//-- Functions related to print
void printIntArray(int *array, int length);
void printDoubleArray(double *array, int length, int digit);
void printGalaxyInfo(peak_param *peak);
void printFilterArray(filter_t *array, int length);
void printParam(cosmo_hm *cmhm, peak_param *peak);
void printParam_complete(cosmo_hm *cmhm, peak_param *peak);


#endif

