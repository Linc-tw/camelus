

  /*******************************
   **  peakParameters.h		**
   **  Chieh-An Lin		**
   **  Version 2015.04.05	**
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
#define NUMBER_HALOS_MAX     5000000
#define NUMBER_GALAXIES_MAX  5000000
#define CUTOFF_FACTOR_HALO   3.0
#define CUTOFF_FACTOR_FILTER 2.2

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


typedef enum {rectangle=0, circle=1, aardvark_hPatch04=2, aardvark_gPatch086=3} field_t;
#define NB_FIELD_T 4
#define STR_FIELD_T(i) ( \
  i==0 ? "rectangle" : \
  i==1 ? "circle" : \
  i==2 ? "hPatch04" : \
  i==3 ? "gPatch086" : \
  "")

typedef enum {gauss=0} filter_t;
#define NB_FILTER_T 1
#define STR_FILTER_T(i) ( \
  i==0 ? "gauss" : \
  "")

typedef enum {kappa_map=0, K_map=1, noise_map=2, N_map=3, kn_map=4, KN_map=5, noise_list=6, peak_list=7} mapType_t;
#define NB_MAPTYPE_T 8
#define STR_MAPTYPE_T(i) ( \
  i==0 ? "kappa map" : \
  i==1 ? "K map"  : \
  i==2 ? "Noise map" : \
  i==3 ? "N map" : \
  i==4 ? "k+n map" : \
  i==5 ? "K+N map" : \
  i==6 ? "Noise list" : \
  i==7 ? "Peak list" : \
  "")
#define COMMENT_MAPTYPE_T(i) ( \
  i==0 ? "unsmoothed noiseless map" : \
  i==1 ? "smoothed noiseless map"  : \
  i==2 ? "unsmoothed noise map" : \
  i==3 ? "smoothed noise map" : \
  i==4 ? "unsmoothed noisy map" : \
  i==5 ? "smoothed noisy map" : \
  i==6 ? "" : \
  i==7 ? "" : \
  "")

typedef void rand_pos_fct(gsl_rng*, double[2]);

//-- Peak parameters
typedef struct {
  //----------------------------------------------------------------------
  //-- Loaded part
  
  //-- Field
  field_t field;      //-- Geometry (only 'rectangle' is available in v1.2)
  double Omega[2];    //-- [arcmin] Field size (theta_x, theta_y)
  
  //-- Halos
  double z_halo_max;  //-- Maxmimum redshift for halos
  int N_z_halo;       //-- Number of sampling redshift bins
  double M_min;       //-- [M_sol/h] Minimum sampling mass
  double M_max;       //-- [M_sol/h] Maximum sampling mass
  
  //-- Galaxies
  int doKappa;        //-- 0 = gamma, 1 = kappa
  double z_s;         //-- If z_s > 0, all source galaxies fixed at z_s
                      //-- If z_s < 0, use nofz_hm.dat as redshift distribution law parameters
  int doRandGalPos;   //-- 0 = regular, 1 = random
                      //-- If z_s < 0, doRandGalPos is automatically set to 1 
  double n_gal;       //-- [arcmin^-2] Galaxy number density
  double sigma_eps;   //-- Ellipticity dispersion
  
  //-- Map, filter & peaks selection
  double theta_pix;   //-- [arcmin] Pixel size
  filter_t filter;    //-- Filter type
  double theta_G;     //-- [arcmin] Filter size, only used if filter = gauss
  //double fillingFactor; //-- Minimum galaxy filling factor necessary to consider using the pixel in peak counts 
  
  //-- ABC
  int ABC_d;          //-- Parameter dimension
  int ABC_p;          //-- Number of particles
  double ABC_r_stop;  //-- Shutoff success rate
  char ABC_summ[32];  //-- Summary type
  
  //-- Label
  char *simulName;
  int paperII_index;  //-- Paper II only
  
  //----------------------------------------------------------------------
  //-- Precomputed part
  
  //-- Modes
  int printMode;      //-- 0 = detailed, 1 = ABC, 2 = log file
  int doNewNoise;     //-- 0 = use information contained in nMap, 1 = generate new noise map for each realization
  
  //-- Others
  gsl_rng *generator;
  double dlogM;       //-- Bin width for mass function
  int nbMassBins;     //-- Number of mass bins
  double area;        //-- [arcmin^2] Field area
  double w_s;         //-- [Mpc/h] Comoving radial distance of the source
  double D_s;         //-- [Mpc/h] Angular diameter distance of the source
  double sigma_noise; //-- Expected noise level for a Gaussian filter, 0.024228320599590476, for aardvark
  double sigma_pix;   //-- Noise dispersion for each pixel, 0.28284271247461901 for aardvark
  double theta_G_sq;  //-- [arcmin^2]
  
  //-- Only used if field = rectangle
  double s;           //-- [pix] theta_G / sigma_pix
  int resol[2];       //-- [pix] round(Omega / theta_pix)
  int bufferSize;     //-- [pix] 0-padding size
  int FFTSize;        //-- [pix] FFT table size
  
  //-- Only used if field = aardvark_hPatch04 or aardvark_gPatch086
  rand_pos_fct *randPosFct;         //-- Function pointer to customize for an HEALPix patch
  int cut;                          //-- aardvark only
  int nbPatches, *patchList;        //-- aardvark only
  int patchOrder, RTPatch;          //-- aardvark only
  int fastNb, noiseNb, doSmoothing; //-- aardvark only
  double cutoff;      //-- [arcmin] Cutoff size for filter
  double cutoff_sq;   //-- [arcmin^2]
  //-- bufferSize, resol[2] are also used when field = aardvark_hPatch04 or aardvark_gPatch086
} peak_param;


//-- Functions related to cosmo_hm
cosmo_hm *initialize_cosmo_hm_default(error **err);
void read_cosmo_hm(char name[], cosmo_hm **cmhm, error **err);
void outputCosmoParam(FILE *file, cosmo_hm *cmhm, peak_param *peak);
cosmo_hm *updateCmhm(cosmo_hm *old, double Omega_m, double sigma_8, error **err);

//-- Functions related to peak_param
//peak_param *initialize_peak_param_default(error **err);
void free_peak_param(peak_param *peak);
void read_peak_param(char name[], peak_param *peak, error **err);
void set_peak_param(cosmo_hm *cmhm, peak_param *peak, error **err);
void printParam(cosmo_hm *cmhm, peak_param *peak);
void printParam_ABC(cosmo_hm *cmhm, peak_param *peak);
void printParam_complete(cosmo_hm *cmhm, peak_param *peak);

//-- Functions related to position randomization
void randPos_aardvark_hPatch04(gsl_rng *generator, double pos[2]);
void randPos_aardvark_gPatch086(gsl_rng *generator, double pos[2]);


#endif

