

  /*******************************
   **  smoothing.h		**
   **  Chieh-An Lin		**
   **  Version 2016.03.20	**
   *******************************/


#ifndef __smooth__
#define __smooth__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"

#ifndef __releaseMenu__
  #include "FITSFunctions.h"
#endif


typedef struct {
  int N1, N2;        //-- Number of pixels of each side
  int length;        //-- If type = peak_list, length = nb of peaks, otherwise length = N1 * N2
  double theta_pix;  //-- [arcmin/deg] Pixel size
  double theta_pix_inv; //-- [arcmin^-1/deg^-1] Inverse of pixel size
  double limit[4];   //-- [arcmin/deg] (theta_x_min, theta_x_max, theta_y_min, theta_y_max) in arcmin or (RA_min, RA_max, DEC_min, DEC_max) in deg
  double center[2];  //-- [arcmin/deg] Center of the map;
  mapType_t type;    //-- Map type (defined in peakParameters.h)
  double *value1;    //-- Either kappa, gamma_1, or g_1
  double *value2;    //-- Either empty, gamma_2, or g_2
} map_t;

typedef void DC_fct(peak_param*, gal_list*, FFT_t**, double*, int);


//-- Functions related to map_t
map_t *initialize_map_t(int N1, int N2, double theta_pix, error **err);
void free_map_t(map_t *kMap);
void getPixPos_map_t(map_t *kMap, double pos[2], int i, int j);
void getPixPos(double pos[2], double limit[4], double theta_pix, int i, int j);
void read_map_t(char *name, peak_param *peak, map_t *kMap, mapType_t type, error **err);
void output_map_t(FILE *file, map_t *kMap);

//-- Functions related to smoothing by FFT
void fillGaussianKernel(fftw_complex *kernel, int length, int M, double scaleInPix);
double starlet_2D(double x, double y);
void fillStarletKernel(fftw_complex *kernel, int length, int M, double scaleInPix);
void makeKernel(peak_param *peak, FFT_arr *FFTSmoother);
void smoothByFFT_arr(peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother);

//-- Functions related to direct convolution
void DCForPair_kappa(peak_param *peak, gal_list *gList, FFT_t **DCSmooArr, double *pixPos, int index_FFT);
void DCForPair_gamma(peak_param *peak, gal_list *gList, FFT_t **DCSmooArr, double *pixPos, int index_FFT);
void DCForPixel(peak_param *peak, gal_map *gMap, FFT_arr *DCSmoother, DC_fct *DCForPairFct, int i_pix, int j_pix);
void DCForMap(peak_param *peak, gal_map *gMap, FFT_arr *DCSmoother, DC_fct *DCForPairFct);
void smoothAndBinByDC(peak_param *peak, gal_map *gMap, FFT_arr *DCSmoother);

//-- Functions related to binning and noise
void galaxyBinning(gal_map *gMap, FFT_t *FFTSmoo);
void copyToAllBefore(FFT_arr *FFTSmoother);
void addNoiseToTable(peak_param *peak, gal_map *gMap, fftw_complex *table);

//-- Functions related to linear KS inversion
void linearKS(fftw_complex *table, int M);
void invertByLinKS(peak_param *peak, FFT_t *smoo);
void invertByLinKS_arr(peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother);
void linKSAndFFT(peak_param *peak, FFT_t *FFTSmoo);
void linKSAndFFT_arr(peak_param *peak, FFT_arr *FFTSmoother);

//-- Functions related to iterative KS inversion
void iterativeKS(fftw_complex *reducedShear, fftw_complex *table, fftw_plan forward, fftw_plan backward, map_t *kMap, double FFTNormFactor, int M);
void zeroPadding(fftw_complex *table, int N1, int N2, int M);
void invertByIterKS(peak_param *peak, FFT_t *smoo, map_t *kMap);
void invertByIterKS_arr(peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, map_t *kMap);

//-- Functions related to nonlinear filtering
void MRLensFiltering(peak_param *peak, char executable[], char option[], char input[], char output[], char log[]);

//-- Functions related to output
void outputMap(char name[], cosmo_hm *cmhm, peak_param *peak, map_t *kMap);
void outputMapFromTable(char name[], cosmo_hm *cmhm, peak_param *peak, fftw_complex *table, mapType_t type);
void outfitsMapFromTable(char name[], peak_param *peak, fftw_complex *table, map_t *kMap);
void outputMask(char name[], peak_param *peak, gal_map *gMap, map_t *kMap, error **err);

//-- Functions related to map making
void pixelization(peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, error **err);
void makeMap(peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, map_t *kMap, error **err);
void makeTrueMap(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, map_t *kMap, int semiTruth, error **err);
void makeMapAndOutputAll(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, map_t *kMap, error **err);

//-- Main functions
void doKMap(char fileName[], cosmo_hm *cmhm, peak_param *peak, int doNoise, error **err);

//-- New functions
void makeMapAndOutputAll2(char fileName[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, map_t *kMap, error **err);

#endif

