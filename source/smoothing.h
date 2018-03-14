

  /*******************************************************
   **  smoothing.h					**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_SMOOTHING__
#define __CAMELUS_SMOOTHING__

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#include "parameters.h"
#include "haloSampling.h"
#include "galaxySampling.h"
#include "rayTracing.h"


typedef struct {
  int N1, N2;           //-- Number of pixels of each side
  int length;           //-- If type = peak_list, length = nb of peaks, otherwise length = N1 * N2
  double theta_pix;     //-- [rad] Pixel size
  double theta_pix_inv; //-- [rad^-1] Inverse of pixel size
  map_t type;           //-- Map type (defined in parameters.h)
  double *value1;       //-- Either kappa, gamma_1, or g_1
  double *value2;       //-- Either empty, gamma_2, or g_2
} signal_map;

typedef void DC_fct(peak_param*, gal_list*, FFT_t**, double*, int);


//-- Functions related to signal_map
signal_map *initialize_signal_map(int N1, int N2, double theta_pix, error **err);
void free_signal_map(signal_map *kMap);
void read_signal_map(char *name, peak_param *pkPar, signal_map *kMap, map_t type, error **err);
void getPixPos(double pos[2], double theta_pix, int i, int j);

//-- Functions related to smoothing by FFT
void fillGaussianKernel(fftw_complex *kernel, int length, int M, double scaleInPix, int doVar);
double starlet_2D(double x, double y);
void fillStarletKernel(fftw_complex *kernel, int length, int M, double scaleInPix, int doVar);
void fillMApTanhKernel(fftw_complex *kernel, int length, int M, double scaleInPix, int doVar);
void fillMApGammaTKernel(fftw_complex *kernel, int length, int M, double scaleInPix1, double scaleInPix2, int doVar);
void makeKernel(peak_param *pkPar, FFT_arr *FFTSmoother);
void smoothByFFT_arr(peak_param *pkPar, FFT_arr *FFTSmoother);

//-- Functions related to direct convolution
void DCForPair_kappa(peak_param *pkPar, gal_list *gList, FFT_t **smooArr, double pixPos[2], int index_FFT);
void DCForPair_gamma(peak_param *pkPar, gal_list *gList, FFT_t **smooArr, double pixPos[2], int index_FFT);
void DCForPixel(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, DC_fct *DCForPairFct, double pixPos[2], int index_FFT, int i_pix, int j_pix);
void DCForMap(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, DC_fct *DCForPairFct);
void smoothAndBinByDC(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother);

//-- Functions related to linear KS inversion
void gammaToKappa(fftw_complex *table, int M);
void invertByLinKS(FFT_t *smoo);
void linKSAndFFT(FFT_t *smoo);
void linKSAndFFT_arr(peak_param *pkPar, FFT_arr *FFTSmoother);

//-- Functions related to iterative KS inversion
void iterativeKS(fftw_complex *reducedShear, fftw_complex *table, fftw_plan forward, fftw_plan backward, signal_map *kMap, int M);
void zeroPadding(fftw_complex *table, int N1, int N2, int M);
void invertByIterKS(FFT_t *smoo, signal_map *kMap);

//-- Functions related to Seitz-Schneider inversion
void kappaToGamma(fftw_complex *table, int M);
void SeitzSchneider(fftw_complex *reducedShear, fftw_complex *table, fftw_plan forward, fftw_plan backward, signal_map *kMap, int M);
void invertBySS(FFT_t *smoo, signal_map *kMap);

//-- Functions related to ASCII output
void outAscii_signal_map(FILE *file, signal_map *kMap, int field);
void outAscii_fftw_complex(FILE *file, peak_param *pkPar, fftw_complex *table, map_t type);
void outAsciiNoiseInfo(FILE *file, peak_param *pkPar);
void outAsciiFilterInfo(FILE *file, peak_param *pkPar, map_t type, int filterInd);
void outAsciiAllFilterInfo(FILE *file, peak_param *pkPar);
void outAsciiMap(char name[], cosmo_hm *chPar, peak_param *pkPar, signal_map *kMap, int filterInd, error **err);
void outAsciiMapFromTable(char name[], cosmo_hm *chPar, peak_param *pkPar, fftw_complex *table, map_t type, int filterInd, error **err);
void outAsciiMask(char name[], peak_param *pkPar, gal_map *gMap, signal_map *kMap, error **err);

//-- Functions related to FITS output
#ifdef __CAMELUS_USE_FITS__
void outFits_signal_map(FITS_t *fits, signal_map *kMap, int field, double factor);
void outFits_fftw_complex(FITS_t *fits, peak_param *pkPar, fftw_complex *table, int type);
void outFitsNoiseInfo(FITS_t *fits, peak_param *pkPar);
void outFitsFilterInfo(FITS_t *fits, peak_param *pkPar, map_t type, int filterInd);
void outFitsAllFilterInfo(FITS_t *fits, peak_param *pkPar);
#endif
void outFitsMap(char name[], cosmo_hm *chPar, peak_param *pkPar, signal_map *kMap, int filterInd);
void outFitsMapFromTable(char name[], cosmo_hm *chPar, peak_param *pkPar, fftw_complex *table, map_t type, int filterInd);
void outputMapFromTable(char name[], cosmo_hm *chPar, peak_param *pkPar, fftw_complex *table, map_t type, int filterInd, error **err);

//-- Functions related to map making
void gMapToSmoother(gal_map *gMap, FFT_t *smoo);
void galaxyBinning(peak_param *pkPar, gal_map *gMap, gal_map *gMap2, FFT_t *smoo, error **err);
void invertForTrueMap(peak_param *pkPar, FFT_t *smoo, signal_map *kMap);
void pixelization(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, error **err);
void inversion(peak_param *pkPar, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap);
void makeMaps(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap, error **err);
void makeMapsAndOutput(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, gal_map *gMap2, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap, error **err);

//-- Main functions
void doKMap(cosmo_hm *chPar, peak_param *pkPar, error **err);

#endif

