

  /*******************************
   **  constraint.h		**
   **  Chieh-An Lin		**
   **  Version 2015.12.09	**
   *******************************/


#ifndef __constraint__
#define __constraint__

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"


typedef struct {
  int d;                    //-- Dimension of data vector
  int N;                    //-- Number of realizations
  double *matrix;           //-- N * d data matrix
  gsl_vector *x_mod;        //-- Model prediction, mean of several realizations of Camelus
  gsl_vector *x_obs;        //-- Data vector from observation 
  gsl_matrix *cov;          //-- Covariance matrix for x^mod, debiased
  gsl_matrix *cov2;         //-- Just a copy of cov
  gsl_matrix *invCov;       //-- Inverse of cov, debiased
  gsl_permutation *perm;    //-- Used for matrix inversion
  KDE_arr *estArr;          //-- Array of kernel density estimators
  
  //-- Stockage
  gsl_vector *intermediate; //-- To stock invCov * (x^mod - x^obs)
  gsl_vector *B_arr;        //-- To stock B_i for copula
} likelihood_t;


//-- Functions related to data matrix
void fillMultiscale(peak_param *peak, hist_t *hist, double_mat *multiscale);
void multiscaleFromMassFct(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask,
			   FFT_arr *smoother, map_t *kMap, FFT_arr *varArr, double_arr *peakList, hist_t *hist, hist_t *hist2, double_mat *multiscale, error **err);
void fillRealization(peak_param *peak, double_mat *multiscale, double *matrix);
void outputMultiscale(char name[], double_mat *multiscale);
void outputDataMat(char name[], double_mat *dataMat);
//void outfitsCosmoParam(FITS_t *fits, cosmo_hm *cmhm, peak_param *peak);
//void outfits_double_mat(FITS_t *fits, double_mat *dataMat, int N_bin, int nbFilters);
//void outfitsDataMat(char name[], cosmo_hm *cmhm, peak_param *peak, double_mat *dataMat);

//-- Main function
void doMultiscale(cosmo_hm *cmhm, peak_param *peak, error **err);
void doDataMatrix(cosmo_hm *cmhm, peak_param *peak, int N, error **err);


/*
//-- Functions related to likelihood_t
likelihood_t *initialize_likelihood_t(int d, int N, error **err);
void free_likelihood_t(likelihood_t *LLH);
void set_likelihood_t(likelihood_t *LLH);
void setKDE_likelihood_t(likelihood_t *LLH);
void fillObservation_likelihood_t(char name[], likelihood_t *LLH, error **err);

//-- Functions related to chi^2 computation
double computeGaussianChi2(likelihood_t *LLH);
double computeCopulaChi2(likelihood_t *LLH);

void doChi2(cosmo_hm *cmhm, peak_param *peak, int N, error **err);
*/


#endif

