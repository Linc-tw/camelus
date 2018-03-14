

  /*******************************************************
   **  multiscale.h					**
   **  Version 2018.03.07				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_MULTISCALE__
#define __CAMELUS_MULTISCALE__

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_cdf.h>

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#include "parameters.h"
#include "haloSampling.h"
#include "galaxySampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"


typedef struct {
  sampler_arr *hSampArr;
  interpolator_t *k1Inter;
  halo_map *hMap;
  sampler_t *gSamp;
  gal_map *gMap;
  mask_map *mask;
  FFT_arr *FFTSmoother;
  FFT_arr *DCSmoother;
  signal_map *kMap;
  FFT_arr *variance;
  double_arr *peakList;
  hist_t *nuHist;
  double_mat *multiscale;
} pipeline_t;


//-- Functions related to data matrix
void addToMultiscale(hist_t *hist, double_mat *multiscale, int scaleInd);
void mapToMultiscale_FFT(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err);
void mapToMultiscale_DC(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err);
void mapToMultiscaleAndOutput_FFT(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err);
void mapToMultiscaleAndOutput_DC(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err);
void massFctToMultiscale(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, interpolator_t *k1Inter, halo_map *hMap, sampler_t *gSamp, gal_map *gMap, mask_map *mask, 
			 FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, 
			 double_mat *multiscale, int nbPatches, error **err);
void forwardResult(peak_param *pkPar, double_mat *multiscale, double *matrix);
void inverseCovarianceMatrix(double_mat *dataMat, gsl_matrix *invCov);

//-- Functions related to output
void outAsciiMultiscale(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *multiscale, error **err);
void outFitsMultiscale(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *multiscale);
void outputMultiscale(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *multiscale, error **err);
void outAsciiDataMat(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *dataMat, error **err);
void outFitsDataMat(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *dataMat);
void outputDataMat(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *dataMat, error **err);
void outAsciiInvCov(char name[], cosmo_hm *chPar, peak_param *pkPar, gsl_matrix *invCov, error **err);
void outFitsInvCov(char name[], cosmo_hm *chPar, peak_param *pkPar, gsl_matrix *invCov);
void outputInvCov(char name[], cosmo_hm *chPar, peak_param *pkPar, gsl_matrix *invCov, error **err);

//-- Functions related to pipeline_t
pipeline_t *initialize_pipeline_t(cosmo_hm *chPar, peak_param *pkPar, error **err);
void free_pipeline_t(pipeline_t *pipe);

//-- Main function
void doMultiscale(cosmo_hm *chPar, peak_param *pkPar, error **err);
void doDataMatrix(cosmo_hm *chPar, peak_param *pkPar, int N, error **err);

#endif

