

  /*******************************
   **  peakParameters.h		**
   **  Chieh-An Lin		**
   **  Version 2015.12.08	**
   *******************************/


#ifndef __peakParam__
#define __peakParam__

#define _GNU_SOURCE

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
#define CUTOFF_FACTOR_FILTER    2.2
#define CUTOFF_FACTOR_FILTER_SQ 4.84
#define WAVELET_SCALE_MAX       20.0 //-- [arcmin]
#define STARLET_NORM1           0.97873802046145453
#define STARLET_NORM2           0.65166118860701647

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


typedef enum {rectangle=0, circle=1, aardvark_hPatch04=2, aardvark_gPatch086=3} field_t;
#define NB_FIELD_T 4
#define STR_FIELD_T(i) ( \
  i==0 ? "rectangle" : \
  i==1 ? "circle" : \
  i==2 ? "hPatch04" : \
  i==3 ? "gPatch086" : \
  "")

typedef enum {gauss=0, star=1, mrlens=2} filter_t;
#define NB_FILTER_T 3
#define STR_FILTER_T(i) ( \
  i==0 ? "gauss" : \
  i==1 ? "star" : \
  i==2 ? "mrlens" : \
  "")

typedef enum {
  kappa_map=0, K_map=1,  kn_map=2,  KN_map=3, noise_map=4,    N_map=5, 
  gamma_map=6, R_map=7,  re_map=8,  RE_map=9, epsilon_map=10, E_map=11,
  g_map=12,    G_map=13, ge_map=14, GE_map=15, 
  peak_list=16, noise_list=17, mask_map=18
} mapType_t;
#define NB_MAPTYPE_T 19
#define STR_MAPTYPE_T(i) ( \
  i==0  ? "kappa map"   : i==1 ?  "K map"  : \
  i==2  ? "k+n map"     : i==3 ?  "K+N map" : \
  i==4  ? "n map"       : i==5 ?  "N map" : \
  i==6  ? "gamma map"   : i==7 ?  "R map" : \
  i==8  ? "r+e map"     : i==9  ? "R+E map" : \
  i==10 ? "epsilon map" : i==11 ? "E map" : \
  i==12 ? "g map"       : i==13 ? "G map" : \
  i==14 ? "g+e map"     : i==15 ? "G+E map" : \
  i==16 ? "Peak list" : \
  i==17 ? "Noise list" : \
  i==18 ? "Mask" : \
  "")

//-- Peak parameters
typedef struct {
  //----------------------------------------------------------------------
  //-- Customized part
  
  //-- Halos
  double z_halo_max;     //-- Maxmimum redshift for halos
  int N_z_halo;          //-- Number of sampling redshift bins
  double M_min;          //-- [M_sol/h] Minimum sampling mass
  double M_max;          //-- [M_sol/h] Maximum sampling mass
  
  //-- Galaxies
  double z_s;            //-- If z_s > 0, all source galaxies fixed at z_s
                         //-- If z_s < 0, use nofz_hm.dat as redshift distribution law parameters
  int doRandGalPos;      //-- 0 = regular, 1 = random
                         //-- If z_s < 0, doRandGalPos is automatically set to 1 
  double n_gal;          //-- [arcmin^-2] Galaxy number density
  double sigma_eps;      //-- [-] Ellipticity dispersion, sigma_eps^2 = <epsilon_1^2> + <epsilon_2^2> = 2 sigma_kappa^2
  int doKappa;           //-- 0 = gamma, 1 = kappa, 2 = g, 3 = g with linear KS
  int doMask;            //-- 0 = without, 1 = CFHTLenS W1
  
  //-- Field. map, & filter
  field_t field;         //-- Geometry (only 'rectangle' is available in v1.2)
  double Omega[2];       //-- [arcmin] Field size (theta_x, theta_y)
  double theta_pix;      //-- [arcmin] Pixel size
  int nbFilters;         //-- Number of filters
  filter_t *filter;      //-- Filter type
  double *scale;         //-- [arcmin] Filter size
  
  //-- Peak historgram
  int N_nu;              //-- [-] Number of S/N bins
  double *nu_bin;        //-- S/N bins
  int N_kappa;           //-- [-] Number of kappa bins
  double *kappa_bin;     //-- kappa bins
  
  //-- ABC
  int ABC_Q;             //-- Number of particles
  double ABC_r_stop;     //-- Shutoff accept rate
  char ABC_summ[64];     //-- Summary type
  
  //-- Label
  char *simulName;       //-- Name of the parameter set
  
  //----------------------------------------------------------------------
  //-- Default part
  
  double dlogM;          //-- [-] Bin width for mass function
  double dz_gal;         //-- [-] Galaxy redshift binwidth
  double theta_CCD_inv;  //-- [arcmin^-1] Inverse of CCD mask pixel size
  
  //----------------------------------------------------------------------
  //-- Precomputed part
  
  double area;           //-- [arcmin^2] Field area
  double w_s;            //-- [Mpc/h] Comoving radial distance of the source
  double D_s;            //-- [Mpc/h] Angular diameter distance of the source
  double dz_halo;        //-- [-] Halo redshift binwidth
  int N_M;               //-- Number of mass bins
  int N_z_gal;           //-- Number of galaxy redshift bins
  double sigma_half;     //-- [-] sigma_eps / sqrt(2)
  double theta_pix_inv;  //-- [arcmin^-1] Inverse of pixel size
  double sigma_pix;      //-- [-] Noise dispersion for each pixel
  gsl_rng *generator;    //-- Random number generator
  
  //-- Filter-related
  int nbLinFilters;      //-- Number of linear filters
  int nbNonlinFilters;   //-- Number of nonlinear filters
  int smootherSize;      //-- Length of the smoother
  filter_t *linFilter;   //-- Linear filter type
  double *linScale;      //-- [arcmin] Linear filter size
  int *nonlinScale;      //-- [-] Maximum scale order to compute of nonlinear filters
  double *sigma_noise;   //-- [-] Expected noise level for linear filters, WARNING used in paperIII, but not used in the released version
  double *scale_invSq;   //-- [arcmin^-2] Inverse squared size of linear filters
  double *cut_sq;        //-- [arcmin^2] Squared cutoff size of linear filters
  
  //-- Only used if field = rectangle
  int resol[2];          //-- [pix] round(Omega / theta_pix)
  double *linScaleInPix; //-- [pix] Linear filter size in pixel
  int bufferSize;        //-- [pix] 0-padding size
  int FFTSize;           //-- [pix] FFT table size
  double FFTNormFactor;  //-- [-] FFT normalization factor
  double limits[4];      //-- [arcmin] Field limits (theta_x_min, theta_x_max, theta_y_min, theta_y_max)
  
  //----------------------------------------------------------------------
  //-- Running part
  
  double *weight;        //-- [-] Weight for direct convolution, WARNING not initialized
  
  //-- Mode
  int printMode;         //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  int doNoise;           //-- 0 = noiseless, 1 = noisy
  int doSmoothing;       //-- 0 = without, 1 = FFT smoothing, (2 = direct convolution), 3 = nonlinear filtering, 4 = FFT + nonlinear
  
  //-- Pipeline index
  int cosmoInd;          //-- Index for cosmology, Paper II and Paper III only
  int scaleInd;          //-- Index for scale, for multiscale analysis
  int realizationInd;    //-- Index for realization
  
  //-- MPI
  int MPISize;           //-- Number of MPI processors
  int MPIInd;            //-- Index for MPI processors
  
  //----------------------------------------------------------------------
  //-- Paper I only
  
  int cut;                   //-- aardvark only
  int nbPatches, *patchList; //-- aardvark only
  int patchOrder, RTPatch;   //-- aardvark only
  int fastId, noiseId;       //-- aardvark only
  
  //----------------------------------------------------------------------
} peak_param;


//-- Functions related to cosmo_hm
cosmo_hm *initialize_cosmo_hm_default(error **err);
void read_cosmo_hm(char name[], cosmo_hm **cmhm, error **err);
void outputCosmoParam(FILE *file, cosmo_hm *cmhm, peak_param *peak);
cosmo_hm *updateCmhm(cosmo_hm *old, double Omega_m, double sigma_8, double w0_de, error **err);

//-- Functions related to peak_param
void free_peak_param(peak_param *peak);
void read_peak_param(char name[], peak_param *peak, error **err);
void set_peak_param(cosmo_hm *cmhm, peak_param *peak, error **err);

//-- Functions related to print
void printLineDoubleArray(double *array, int length, int digit);
void printGalaxyInfo(peak_param *peak);
void printParam(cosmo_hm *cmhm, peak_param *peak);
void printParam_ABC(cosmo_hm *cmhm, peak_param *peak);
void printParam_complete(cosmo_hm *cmhm, peak_param *peak);


#endif

