

  /***************************************
   **  peakSelection.h			**
   **  Chieh-An Lin, François Lanusse	**
   **  Version 2015.03.25		**
   ***************************************/


#ifndef __peakSelect__
#define __peakSelect__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"


typedef struct {
  int N;                    //-- Number of realizations
  int d;                    //-- Dimension of observable vector
  gsl_vector *X_model;      //-- Observables from model, mean of several realizations for Camelus
  gsl_vector *intermediate; //-- To stock invCov * (X_model - X_obs)
  gsl_vector *X_obs;        //-- Observables from observation 
  gsl_matrix *cov;          //-- Covariance matrix for X_model, debiased
  gsl_matrix *cov2;         //-- Just a copy of cov
  gsl_matrix *invCov;       //-- Inverse of cov, debiased
  gsl_permutation *perm;    //-- Used for matrix inversion
} chi2_t;


//-- Functions related to peak selection
int isPeak(double *kappa, int N1, int i, int j);
void selectPeaks(peak_param *peak, map_t *kMap, double_arr *peakList, error **err);
void cutSmallPeaks(double_arr *peakList, double nu_min);
void outputPeakList(char name[], peak_param *peak, double_arr *peakList);
void peakListFromMassFct(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, gal_map *gMap, 
			 map_t *kMap, map_t *nMap, FFT_t *transformer, double_arr *peakList, error **err);
void makeHist(double_arr *peakList, hist_t *hist, int silent);
void outputHist(char name[], hist_t *hist);
void inputHist(char name[], hist_t *hist, error **err);

//-- Functions related to chi2_t
chi2_t *initialize_chi2_t(int N, int d, error **err);
void free_chi2_t(chi2_t *chichi);
void update_chi2_t(chi2_t *chichi, hist_t *obsHist, double *dataMat);
double execute_chi2_t(chi2_t *chichi);

//-- Main functions
void doPeakList(char KNMap[], cosmo_hm *cmhm, peak_param *peak, error **err);
void doPeakList_repeat(cosmo_hm *cmhm, peak_param *peak, int N, error **err);
//void doPeakHistForPMC(cosmo_hm *cmhm, peak_param *peak, hist_t *hist, double *dataMat, int N, error **err);


#endif

