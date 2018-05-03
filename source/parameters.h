

  /*******************************************************
   **  parameters.h					**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_PARAMETERS__
#define __CAMELUS_PARAMETERS__

#include <nicaea/errorlist.h>
#include <nicaea/io.h>
#include <nicaea/cosmo.h>
#include <nicaea/nofz.h>
#include <nicaea/halomodel.h>
#include <nicaea/hod.h>

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#ifdef __CAMELUS_USE_HEALPIX__
#include "HEALPixFunctions.h"
#endif


//-- Computational constants
#define NUMBER_HALOS_MAX        5000000
#define NUMBER_GALAXIES_MAX     5000000
#define CUTOFF_FACTOR_HALO      3.0
#define CUTOFF_FACTOR_HALO_SQ   9.0
#define CUTOFF_FACTOR_GAUSSIAN  2.2
#define CUTOFF_FACTOR_M_AP_TANH 1.2
#define CUTOFF_FACTOR_HEALPIX   1.01 //-- The maximal field size for a distance dilatation of 1% is 13.95 deg.
#define NORM1_STARLET           0.97873802046145453 //-- W(x, y) = psi(x, y)
#define NORM2_STARLET           0.65166118860701647
#define NORM1_M_AP_TANH         5.4323267630348244  //-- W(x) = tanh(x / 0.1) / [x (1 + exp(5 - 150 x) + exp(-47 + 50 x))]
#define NORM2_M_AP_TANH         3.8414742595565707
#define FILLING_THRESHOLD_RATIO 0.5
#define EPS_NUM                 1e-12 //-- Numerical tolerance
#define EPS_ADJ                 4e-12 //-- Adjust amplitude

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

#define STRING_TO_ENUM(element, selement, type, stype, j, N, err) \
  element = -1;\
  for (j=0; j<N; j++) {\
    if (strcmp(selement, stype((type)j))==0) element = (type)j;\
  }\
  testError((int)element==-1, conf_undef, "Parametrization '%s' not found in type %s, cannot assign value of type %s", *err, __LINE__)

// MKDEBUG from peakParameters.h, removed in Linc-tw
typedef enum {rectangle=0, circle=1, aardvark_hPatch04=2, aardvark_gPatch086=3} field_t;
#define NB_FIELD_T 4
#define STR_FIELD_T(i) ( \
  i==0 ? "rectangle" : \
  i==1 ? "circle" : \
  i==2 ? "hPatch04" : \
  i==3 ? "gPatch086" : \
  "")

typedef enum {gauss=0, star=1, M_ap_tanh=2, mrlens=3} filter_t;
#define NB_FILTER_T 4
#define STR_FILTER_T(i) ( \
  i==0 ? "gauss" : \
  i==1 ? "star" : \
  i==2 ? "tanh" : \
  i==3 ? "mrlens" : \
  "")



typedef enum {
  kappa_map=0, K_map=1,  kn_map=2,  KN_map=3, noise_map=4,    N_map=5, 
  gamma_map=6, R_map=7,  re_map=8,  RE_map=9, epsilon_map=10, E_map=11,
  g_map=12,    G_map=13, ge_map=14, GE_map=15, 
  peak_list=16, noise_list=17, proj_mask=18, HP_mask=19
} map_t;
#define NB_MAP_T 20
#define STR_MAP_T(i) ( \
  i==0  ? "kappa map"       : i==1 ?  "Smoothed kappa map"  : \
  i==2  ? "kappa+noise map" : i==3 ?  "Smoothed (kappa+noise) map" : \
  i==4  ? "Pixel noise map" : i==5 ?  "Smoothed pixel noise map" : \
  i==6  ? "gamma map"       : i==7 ?  "Smoothed gamma map" : \
  i==8  ? "gamma+noise map" : i==9  ? "Smoothed gamma+noise map" : \
  i==10 ? "Shape noise map" : i==11 ? "Smoothed shape noise map" : \
  i==12 ? "g map"           : i==13 ? "Smoothed g map" : \
  i==14 ? "g+noise map"     : i==15 ? "Smoothed g+noise map" : \
  i==16 ? "Peak list" : \
  i==17 ? "Noise list" : \
  i==18 ? "Projected mask" : \
  i==19 ? "HEALPix mask" : \
  "")
  
//-- Peak parameters
typedef struct {
  //----------------------------------------------------------------------
  //-- Customized part
  
  //-- Parameter files
  char pkParPath[STRING_LENGTH_MAX];		//-- [string] Path of the peak parameter file
  char hmParPath[STRING_LENGTH_MAX];		//-- [string] Path of the halo model parameter file
  char cmParPath[STRING_LENGTH_MAX];		//-- [string] Path of the cosmological parameter file
  char seed[64];				//-- [int or string] Seed for generating random numbers
						//--   If seed is int, suggested range = [0, 4294967295]
						//--   Put anything containing 'r' for a random seed
  int verbose;					//-- [int] 0 = default, 1 = all, 2 = no flush, 3 = line mode, 4 = pipeline off & realization on, 5 = MPI, 99 = silent
  
  //-- Field
  int field;					//-- [int] 0 = rectangle, 1 = projected HEALPix, 2 = HEALPix
						//--   If field = 0, ray-tracing with the Cartesian distance then filtering
						//--   If field = 1, ray-tracing with the angular distance then project then FFT
						//--   If field = 2, ray-tracing with the angular distance then direct convolution
  double Omega[2];				//-- [2 float | arcmin] Field size (theta_x, theta_y)
						//--   Ignored if field = 2
  double theta_pix;				//-- [float | arcmin] Pixel size
						//--   Forced to 1.0 if field = 2
  long nside;					//-- [int] N_side of the field
						//--   Ignored if field = 0
  long patch;					//-- [int] Ring order of the field
						//--   Ignored if field = 0
  double rotAng;				//-- [float | deg] Rotation angle of the field after projection
						//--   Forced to 0.0 if field != 1
  int HP_resol;					//-- [int] Resolution for HEALPix regular galaxies or maps
						//--   Ignored if field = 0
						//--   Has to be a power of 2
  
  //-- Halos
  char inHaloCatPath[STRING_LENGTH_MAX];	//-- [string] Path of the input halo catalogue
						//--   Leave blank if sample halos from a mass function
  int inHaloCatCol[5];				//-- [5 int] -1 = skip, 0 = first column, etc.
						//--   Successively read as pos1, pos2, redshift, comoving distance, mass
						//--   Cannot skip pos1, pos2, redshift, mass
						//--   (pos1, pos2) is (theta_x, theta_y) in [arcmin] if field = 0, (RA, DEC) in [deg] otherwise
  double z_halo_min;				//-- [float] Minimum halo redshift
  double z_halo_max;				//-- [float] Maximum halo redshift
  int N_z_halo;					//-- [int] Number of halo redshift bins
  double M_min;					//-- [float | M_sol/h] Minimum halo mass
  double M_max;					//-- [float | M_sol/h] Maximum halo mass
  double dlogM;					//-- [float] Halo mass binwidth
  
  //-- Galaxies
  char inGalCatPath[STRING_LENGTH_MAX];		//-- [string] Path of the input galaxy catalogue
						//--   Leave blank if want a regular grid or sample from a distribution
  int inGalCatCol[7];				//-- [7 int] -1 = skip, 0 = first column, etc.
						//--   Successively read as pos1, pos2, redshift, weight, kappa, e_1, e_2
						//--   Cannot skip pos1, pos2, redshift
						//--   (pos1, pos2) is (theta_x, theta_y) in [arcmin] if field = 0, (RA, DEC) in [deg] otherwise
						//--   weight should be given by shape measurement
						//--   If weight is skipped, consider uniform weighting with 2 / sigma_eps^2
						//--   (e_1, e_2) can also be (gamma_1, gamma_2) or (g_1, g_2)
  double z_s;					//-- [float] Source redshift
						//--   If z_s > 0, fix all source galaxy redshifts at z_s
						//--   If z_s <= 0, use distribution parameters defined in hmParam.par
  double dz_gal;				//-- [float] Galaxy redshift binwidth
  int doRandGalPos;				//-- [int] 0 = regular, 1 = random
						//--   Forced to 1 if z_s <= 0
  double n_gal;					//-- [float | arcmin^-2] Galaxy number density
  double sigma_eps;				//-- [float] Ellipticity dispersion, sigma_eps^2 = <epsilon_1^2> + <epsilon_2^2> = 2 sigma_kappa^2
  
  //-- Masks
  int doMask;					//-- [int] 0 = without, 1 = random mask, 2 = maskPath
  char maskPath[STRING_LENGTH_MAX];		//-- [string] Path of the mask
  int nbHoleTypes;				//-- [int] Number of hole types for random mask
						//--   Ignored if doMask != 1
  double *holeRadius;				//-- [float array | arcmin] Hole radius for random mask
						//--   Ignored if doMask != 1
  double *holeDensity;				//-- [float array | deg^-2] Hole density for random mask
						//--   Ignored if doMask != 1
  double *stripeLRatio;				//-- [float array] Ratio of stripe length to radius
						//--   Ignored if doMask != 1
  double *stripeWRatio;				//-- [float array] Ratio of stripe width to radius
						//--   Ignored if doMask != 1
  
  //-- Lensing
  int doLensing;				//-- [int] 0 = read from inGalCatPath, 1 = compute lensing
						//--   Forced to 0 if inGalCatCol does not skip any of kappa, e_1, or e_2
  int doKappa;					//-- [int] 0 = gamma, 1 = kappa, 2 = g with linear KS, 3 = g with iterative KS, 4 = g with SS
  int doSubtraction;				//-- [int] 0 = without, 1 = subtract kappa_mean, 2 = subtract the mass sheet kappa_1
						//--   Ignored if doKappa = 0
						//--   If doSubtraction = 0, the computed kappa is always positive or zero, since the underdensity is not correctly considered by mass projection
						//--   If doSubtraction = 2, a mass sheet value kappa_1 (depending on the galaxy redshift) is subtracted
  
  //-- Filters
  int doSmoothing;				//-- [int] 0 = without, 1 = binning and FFT, 2 = direct convolution
						//--   Sum for performing simultaneously multiple techniques, e.g. 3 = FFT + DC
  int FFT_nbFilters;				//-- [int] Number of linear filters for FFT
						//--   Forced to 0 if field = 2
						//--   Forced to 0 if doSmoothing without bit 1
  int *FFT_filter;				//-- [int array] 0 = Gaussian, 1 = starlet, 2 = M_ap tanh, 3 = M_ap gamma_t
						//--   2 & 3 only allowed if doKappa != 1
  double *FFT_scale;				//-- [float array | arcmin] Filter size for FFT
  int DC_nbFilters;				//-- [int] Number of linear filters for direct convolution
						//--   Forced to 0 if field = 1
						//--   Forced to 0 if doSmoothing without bit 2
  int *DC_filter;				//-- [int array] 0 = Gaussian, 1 = starlet, 2 = M_ap tanh, 3 = M_ap gamma_t
						//--   0 & 1 only allowed if doKappa = 1
						//--   2 & 3 only allowed if doKappa != 1
  double *DC_scale;				//-- [float array | arcmin] Filter size for direct convolution
  
  //-- Histograms
  int doLocalNoise;				//-- [int] 0 = uniform global noise level, 1 = local noise
  int N_nu;					//-- [int] Number of S/N bins
  double *bin_nu;				//-- [float array] Bin edge values
						//--   The number of bin_nu should be N_nu + 1
  
  //-- Outputs
  char prefix[STRING_LENGTH_MAX];		//-- [string] Prefix for output files
  int doFITS;					//-- [int] 0 = no, 1 = yes
						//--   Forced to 0 if Camelus not linked to cfitsio
  int outHaloCat;				//-- [int] 0 = no, 1 = yes
  int outGalCat;				//-- [int] 0 = no, 1 = yes
  int outMaps;					//-- [int] 0 = no, 1 = yes
  int outTruth;					//-- [int] 0 = no, 1 = yes
						//--   Output noise-free maps
  int outMask;					//-- [int] 0 = no, 1 = yes
  int outPeakList;				//-- [int] 0 = no, 1 = yes
  int outHist;					//-- [int] 0 = no, 1 = yes
  int outMultiscale;				//-- [int] 0 = no, 1 = yes
  
  //-- ABC
  int ABC_f;					//-- [int] Dimension of parameter set
  int *ABC_doParam;				//-- [int array] Parameters to include into constraints
						//--   0 = Omega_m, 1 = Omega_de, 2 = Omega_b, 3 = n_s, 4 = h_100,
						//--   5 = sigma_8, 6 = w0_de,    7 = w1_de,   8 = c_0, 9 = beta_NFW
  int ABC_Q;					//-- [int] Number of particles
  double ABC_r_stop;				//-- [float] Shutoff accept rate
  char ABC_obsPath[STRING_LENGTH_MAX];		//-- [string] Path of the observation data
  int ABC_doCorr;				//-- [int] Distance assumed to be the chi-squared in this version
						//--   0 = uncorrelated inverse covariance, 1 = full inverse covariance
  char ABC_invCovPath[STRING_LENGTH_MAX];	//-- [string] Path of the inverse covariance
  
  //----------------------------------------------------------------------
  //-- Precomputed part
  
  //-- Parameter files
  int cmParTable[128];				//-- [int array] Lookup table for cosmological parameters
  int hmParTable[128];				//-- [int array] Lookup table for halo model parameters
  int pkParTable[128];				//-- [int array] Lookup table for peak parameters
  gsl_rng *generator;				//-- [-] Random number generator
  
  //-- Field
  double area;					//-- [float | arcmin^2] Field area
						//--   If field > 0, this is patch area
  double theta_pix_inv;				//-- [float | arcmin^-1] Inverse of pixel size
  int resol[2];					//-- [int | pix]
						//--   If field = 0, round(Omega / theta_pix)
						//--   If field > 0, 6400 / nside
  
  //-- Halos
  int doInHaloCat;				//-- [int] 0 = random, 1 = inHaloCatPath
  double dz_halo;				//-- [float] Halo redshift binwidth
  int N_M;					//-- [int] Number of mass bins
  
  //-- Galaxies
  int doInGalCat;				//-- [int] 0 = random, 1 = inGalCatPath
  double w_s;					//-- [float | Mpc/h] Comoving radial distance of the source
  double D_s;					//-- [float | Mpc/h] Angular diameter distance of the source
  int N_z_gal;					//-- [int] Number of galaxy redshift bins
  int doNoise;					//-- [int] 0 = without, 1 = noise
  double sigma_half;				//-- [float] sigma_eps / sqrt(2)
  double sigma_pix;				//-- [float] Noise dispersion for each pixel
  
  //-- Filters
  int nbFilters;				//-- [int] Total number of filters
  int smootherSize;				//-- [int] Length of the smoother
  int *filter;					//-- [int array] Filter type
  
  //-- Filters - FFT
  int FFT_firstToInvert;			//-- [int] First scale index to invert gamma/g into kappa
  int FFT_firstToConserve;			//-- [int] First scale index to conserve gamma/g
  int FFT_hasGammaT;				//-- [int] 0 = no, 1 = yes
  double *FFT_scaleInPix;			//-- [int | pix] Linear filter size in pixel
  double *FFT_sigma_noise;			//-- [float array] Expected noise level for linear filters
  
  //-- Filters - DC
  int DC_hasGammaT;				//-- [int] 0 = no, 1 = yes
  double *DC_scaleInPix;			//-- [int | pix] Linear filter size in pixel
  double *DC_scale_inv;				//-- [float array | arcmin^-2] Inverse squared/non-squared size of linear filters
  double *DC_cut;				//-- [float array | arcmin^2] Squared/non-squared cutoff size of linear filters
  double *DC_sigma_noise;			//-- [float array] Expected noise level for linear filters
  
  //-- Maps & histograms
  int bufferSize;				//-- [int | pix] 0-padding size
  int FFTSize;					//-- [int | pix] FFT table size
  double FFTNormFactor;				//-- [float] FFT normalization factor
  int peakFieldResol[2];			//-- [2 int | pix] Peak field resolution
  int peakListLength;				//-- [int] Peak list length
  
  //-- Only used if field = 1 (projected HEALPix) or 2 (HEALPix)
  long nest;					//-- [int] Nest index of the HEALPix patch
  long nsidePix;				//-- [int] nside of pixels from HEALPix regular galaxies or maps
  int HP_length;				//-- [int] Number of pixels from HEALPix regular galaxies or maps
  long HP_first;				//-- [int] Nest index of the first pixel from HEALPix regular galaxies or maps
  int cap;					//-- [int] Decomposition of the HEALPix patch: cap
  int level;					//-- [int] Decomposition of the HEALPix patch: level
  int length;					//-- [int] Decomposition of the HEALPix patch: length
  int off;					//-- [int] Decomposition of the HEALPix patch: offset
  int j;					//-- [int] Decomposition of the HEALPix patch: j
  double z0;					//-- [float] Reference level of z
  double center[4];				//-- [4 float | rad] RA, DEC, sin(DEC), cos(DEC) of the center of the HEALPix patch
  double offset;				//-- [float] Offset of the center of the patch after projection
  
  //----------------------------------------------------------------------
  //-- Running part
  
  int doSurvey;					//-- [int] 0 = no multiplicative correction, 1 = yes (conflit with doKappa < 2) (move to customized?)
  int realizationInd;				//-- [int] Index for realization, in use
  int MPISize;					//-- [int] Number of MPI processors
  int MPIInd;					//-- [int] Index for MPI processors
  
  //----------------------------------------------------------------------
} peak_param;


//-- Functions related to parameter reading
void ignoreComments(char line[]);
int getKeyAndValues(char line[], char key[][256]);
int *makeIntArray(char kv[][256], int count, error **err);
double *makeDoubleArray(char kv[][256], int count, error **err);
void setPathWhichCanBeBlank(char path[], char kv[][256], int count);

//-- Functions related to cosmo_hm
cosmo_hm *initialize_cosmo_hm_default(error **err);
int findKey_cosmo(cosmo *cmPar, peak_param *pkPar, char kv[][256], int count, error **err);
void read_cosmo(char name[], cosmo *cmPar, peak_param *pkPar, error **err);
int findKey_cosmo_hm(cosmo_hm *chPar, peak_param *pkPar, char kv[][256], int count, error **err);
void read_cosmo_hm(char name[], cosmo_hm *chPar, peak_param *pkPar, error **err);
cosmo_hm *reinitialize_cosmo_hm(cosmo_hm *chPar, error **err);
void outAsciiCosmoParam(FILE *file, cosmo_hm *chPar, peak_param *pkPar);
#ifdef __CAMELUS_USE_FITS__
void outFitsCosmoParam(FITS_t *fits, cosmo_hm *chPar, peak_param *pkPar);
#endif

//-- Functions related to peak_param
peak_param *initialize_peak_param(error **err);
void free_peak_param(peak_param *pkPar);
int findKey_peak_param(peak_param *pkPar, char kv[][256], int count, error **err);
void read_peak_param(char name[], peak_param *pkPar, error **err);
int updateFromCommandLine(int argc, char *argv[], cosmo_hm *chPar, peak_param *pkPar, error **err);
void checkParam(peak_param *pkPar, error **err);
void set_peak_param(cosmo_hm *chPar, peak_param *pkPar, error **err);

//-- Functions related to printing
char *printDoKappa(int doKappa);
char *printDoSubtraction(int doSubtraction);
char *printField(int field);
char *printFilter(int filter);
void printIntArray(int *iArr, int length);
void printDoubleArray(double *lfArr, int length, double factor, int digit);
void printGalaxyInfo(peak_param *pkPar, redshift_t *sDParam);
void printFilterArray(int *iArr, int length);
void printParam(cosmo_hm *chPar, peak_param *pkPar);
void printParam_peakCustomized_1(peak_param *pkPar);
void printParam_peakCustomized_2(peak_param *pkPar);
void printParam_peakPrecomputed(peak_param *pkPar);

#endif

