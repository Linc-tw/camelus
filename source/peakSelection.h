

  /***************************************
   **  peakSelection.h			**
   **  Chieh-An Lin, François Lanusse	**
   **  Version 2015.12.09		**
   ***************************************/


#ifndef __peakSelect__
#define __peakSelect__

#include "commonHeader.h"
#include "FITSFunctions.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"


//-- Functions related to local variance
void fillGaussianKernelForVariance(fftw_complex *kernel, int length, int M, double scaleInPix);
void fillStarletKernelForVariance(fftw_complex *kernel, int length, int M, double scaleInPix);
void makeKernelForVariance(peak_param *peak, FFT_arr *varArr);
void computeLocalVariance(gal_map *gMap, FFT_t *var, double sigma_half_sq);
void computeLocalVariance_arr(peak_param *peak, gal_map *gMap, FFT_arr *varArr);
void kappaToSNR(peak_param *peak, gal_map *gMap, FFT_t *smoo, map_t *kMap, FFT_t *var);

//-- Functions related to peak selection
int isPeak(double *kappa, int N1, int i, int j);
int isPeak_float(float *kappa, int N1, int i, int j);
int isPeakForTable(fftw_complex *table, int M, int i, int j);
void selectPeaks(peak_param *peak, map_t *kMap, double_arr *peakList, error **err);
void selectPeaks_mrlens(char name[], peak_param *peak, gal_map *gMap, double_arr *peakList);
void cutSmallPeaks(double_arr *peakList, double nu_min);
void outputPeakList(char name[], peak_param *peak, double_arr *peakList);
void peakListFromMassFct(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask,
			 FFT_arr *smoother, map_t *kMap, FFT_arr *varArr, double_arr *peakList, error **err);

//-- Functions related to histogram
void setHist(peak_param *peak, hist_t *hist);
void setHist2(peak_param *peak, hist_t *hist);
void makeHist(double_arr *peakList, hist_t *hist, int silent);
void outputHist(char name[], hist_t *hist);

//-- Main functions
void doPeakList(char KNMap[], cosmo_hm *cmhm, peak_param *peak, error **err);
void doPeakList_repeat(cosmo_hm *cmhm, peak_param *peak, int N, error **err);


#endif

