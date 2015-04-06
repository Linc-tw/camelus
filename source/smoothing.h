

  /*******************************
   **  smoothing.h		**
   **  Chieh-An Lin		**
   **  Version 2015.02.23	**
   *******************************/


#ifndef __smooth__
#define __smooth__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"


typedef struct {
  int N1, N2;        //-- Number of pixels of each side
  int length;        //-- If type = peak_list, length = nb of peaks, otherwise length = N1 * N2
  double theta_pix;  //-- [arcmin/deg] Pixel size
  double center[2];  //-- [arcmin/deg] Center of the map;
  double limits[4];  //-- [arcmin/deg] (theta_x_min, theta_x_max, theta_y_min, theta_y_max) in arcmin or (RA_min, RA_max, DEC_min, DEC_max) in deg
  mapType_t type;    //-- Map type (defined in peakParameters.h)
  double kappa_mean; //-- kappa_mean = mean of kappa, excluding the buffer area
  double *kappa;
} map_t;


//-- Functions related to map_t
map_t *initialize_map_t(int N1, int N2, double theta_pix, error **err);
void free_map_t(map_t *kMap);
void getPixPos_map_t(map_t *kMap, double pos[2], int i, int j);
void subtractMean_map_t(peak_param *peak, map_t *kMap);
void read_map_t(char *name, peak_param *peak, map_t *kMap, mapType_t type, error **err);
void output_map_t(FILE *file, map_t *kMap);

//-- Functions related to smoothing and noise
void fillGaussianKernel(FFT_t *transformer, double s);
void doFFTSmoothing(peak_param *peak, map_t *kMap, FFT_t *transformer, error **err);
void fillNoise(peak_param *peak, map_t *nMap, u_int64_t seed);
void addNoise(map_t *kMap, map_t *nMap, error **err);

//-- Functions related to map making
void kappaMapFromKappa(gal_map *gMap, map_t *kMap, error **err);
void outputMap(char name[], cosmo_hm *cmhm, peak_param *peak, map_t *kMap);
void makeMap(peak_param *peak, gal_map *gMap, map_t *kMap, map_t *nMap, FFT_t *transformer, error **err);

//-- Main functions
void doKMap(char haloFileName[], cosmo_hm *cmhm, peak_param *peak, error **err);
void doNMap(cosmo_hm *cmhm, peak_param *peak, error **err);
void doKNMap(char KName[], char NName[], cosmo_hm *cmhm, peak_param *peak, error **err);


#endif

