

  /*******************************************************
   **  peakSelection.h					**
   **  Version 2018.03.07				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_PEAK_SELECTION__
#define __CAMELUS_PEAK_SELECTION__

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#include "parameters.h"
#include "haloSampling.h"
#include "galaxySampling.h"
#include "rayTracing.h"
#include "smoothing.h"


//-- Functions related to local variance
void makeKernelForVariance(peak_param *pkPar, FFT_arr *variance);
void fillPixelVariance(gal_map *gMap, FFT_t *var);
void makeLocalVariance(peak_param *pkPar, gal_map *gMap, FFT_arr *variance);
void kappaToSNR_FFT(peak_param *pkPar, gal_map *gMap, FFT_t *FFTSmoo, signal_map *kMap, FFT_t *var, int FFTScaleInd);
void kappaToSNR_DC(peak_param *pkPar, FFT_t *DCSmoo, signal_map *kMap, int DCScaleInd);

//-- Functions related to peak selection
int isPeak(double *kappa, int N1, int i, int j);
int isPeak_float(float *kappa, int N1, int i, int j);
int isPeakForTable(fftw_complex *table, int M, int i, int j);
void selectPeaks(peak_param *pkPar, signal_map *kMap, double_arr *peakList, error **err);
void cutSmallPeaks(double_arr *peakList, double nu_min);
void outAsciiPeakField(FILE *file, peak_param *pkPar);
void outAsciiPeakList(char name[], peak_param *pkPar, double_arr *peakList, int filterInd, error **err);
#ifdef __CAMELUS_USE_FITS__
void outFitsPeakField(FITS_t *fits, peak_param *pkPar);
#endif
void outFitsPeakList(char name[], peak_param *pkPar, double_arr *peakList, int filterInd);

//-- Functions related to histogram
void setHist_nu(peak_param *pkPar, hist_t *hist);
void makeHist(double_arr *peakList, hist_t *hist, int silent);
void outAsciiHistInfo(FILE *file, peak_param *pkPar);
void outAsciiHist(char name[], peak_param *pkPar, hist_t *hist, int filterInd, error **err);
#ifdef __CAMELUS_USE_FITS__
void outFitsHistInfo(FITS_t *fits, peak_param *pkPar);
#endif
void outFitsHist(char name[], peak_param *pkPar, hist_t *hist, int filterInd);

#endif

