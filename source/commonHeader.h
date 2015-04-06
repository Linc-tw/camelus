

  /*******************************
   **  commonHeader.h		**
   **  Chieh-An Lin		**
   **  Version 2015.02.25	**
   *******************************/


#ifndef __commonHdr__
#define __commonHdr__

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <fftw3.h>


//-- Constants
//-- The precision of 64 bits double is 52 bits, and 2^-52 is roughly 2e-16, so we only keep 17 digits in base 10.

//-- Mathematical constants
#define PI                 3.1415926535897932
#define HALF_PI            1.5707963267948966
#define TWO_PI             6.2831853071795865
#define FOUR_PI            12.566370614359172
#define FOUR_PI_OVER_THREE 4.1887902047863910 
#define PI_SQ              9.8696044010893586
#define PI_INV             0.31830988618379067

#define SQRT_2             1.4142135623730950
#define SQRT_2_INV         0.70710678118654752

#define EXP_1              2.7182818284590452
#define LN_10              2.3025850929940457
#define LOG10_E            0.43429448190325183

//-- Unit conversion
#define DEGREE_TO_RADIAN     0.017453292519943296
#define ARCMIN_TO_RADIAN     2.9088820866572160e-04
#define ARCSEC_TO_RADIAN     4.8481368110953599e-06
#define RADIAN_TO_DEGREE     57.295779513082321
#define RADIAN_TO_ARCMIN     3437.7467707849393
#define RADIAN_TO_ARCSEC     206264.80624709636
#define MEGA_PARSEC_TO_METER 3.0856775814671916e+22

//-- Physical constants
#define LIGHT_SPEED          299792458              //-- [m/s]
#define CRITICAL_DENSITY     2.7753619786618317e+11 //-- [M_sol h^2 / Mpc^3], mass density, 3 H^2 / 8 pi G
#define HUBBLE_DISTANCE      2997.92458             //-- [Mpc/h]; LIGHT SPEED is in [m/s], so we should write H_0 in 100000 h [m/s/Mpc]
#define FOUR_PI_G_OVER_C2    6.0135402045290702e-19 //-- [Mpc / M_sol]
#define FULL_SKY             41252.961249419271     //-- [deg^2]

//-- Computational constants
#define STRING_LENGTH_MAX 1024
#define EPS_MIN           4.4408920985006262e-16 //-- 2^-51


typedef struct {
  int length;
  int *array;
} int_arr;

typedef struct {
  int length;
  double *array;
} double_arr;

typedef struct {
  int N1, N2, length;
  double *matrix;
} double_mat;

typedef struct {
  int N1, N2, N3, length;
  double *tensor;
} double_ten3;

typedef struct {
  int length;    //-- Number of points
  double dx;     //-- Interval width
  double *x;     //-- Coordinates
  double *value; //-- Values
} interpolator_t;

typedef struct {
  int length;    //-- Number of points
  double dx;     //-- Interval width
  double *x;     //-- Coordinates
  double *pdf;   //-- Normalized pdf
  double *cdf;   //-- Normalized cdf
  double totPdf; //-- Integration over pdf, before or after normalization (see set_sampler_t)
  double x_mean; //-- Integration over x*p(x), before or after normalization (see set_sampler_t)
} sampler_t;

typedef struct {
  int length;    
  sampler_t **array;
} sampler_arr;

typedef struct {
  int N1, N2, length, length_f;
  double *before, *kernel, *after;             //-- Direct space elements
  fftw_complex *before_f, *kernel_f, *after_f; //-- Fourier space elements
  fftw_plan before_p, kernel_p, after_p;       //-- fftw_plan elements
} FFT_t;

typedef struct {
  int length;      //-- Number of bins
  double dx;       //-- Bin width
  int n_tot;       //-- Total counts
  double *x_lower; //-- Lower limit of each bin
  double x_max;    //-- Maximal limit
  int *n;          //-- Histogram
} hist_t;


//-- Functions related to int_arr
int_arr *initialize_int_arr(int length);
void free_int_arr(int_arr *iArr);
void print_int_arr(int_arr *iArr);

//-- Functions related to double_arr
double_arr *initialize_double_arr(int length);
void free_double_arr(double_arr *fArr);
void print_double_arr(double_arr *fArr);

//-- Functions related to double_mat
double_mat *initialize_double_mat(int N1, int N2);
void free_double_mat(double_mat *fMat);
void print_double_mat(double_mat *fMat);

//-- Functions related to double_ten3
double_ten3 *initialize_double_ten3(int N1, int N2, int N3);
void free_double_ten3(double_ten3 *fTen);

//-- Functions related to interpolator_t
interpolator_t *initialize_interpolator_t(int length);
void free_interpolator_t(interpolator_t *inter);
void print_interpolator_t(interpolator_t *inter);
double execute_interpolator_t(interpolator_t *inter, double x);

//-- Functions related to sampler_t
sampler_t *initialize_sampler_t(int length);
void free_sampler_t(sampler_t *samp);
void print_sampler_t(sampler_t *samp);
void set_sampler_t(sampler_t *samp, int setTotalToOne);
double execute_sampler_t(sampler_t *samp, double x);

//-- Functions related to sampler_arr
sampler_arr *initialize_sampler_arr(int length, int nbPoints);
void free_sampler_arr(sampler_arr *sampArr);

//-- Functions related to FFT_t
FFT_t *initialize_FFT_t(int N1, int N2);
void free_FFT_t(FFT_t *transformer);
void reset_FFT_t(FFT_t *transformer);
void execute_FFT_t(FFT_t *transformer);

//-- Functions related to hist_t
hist_t *initialize_hist_t(int length);
void free_hist_t(hist_t *hist);
void print_hist_t(hist_t *hist);
void set_hist_t(hist_t *hist, double x_min, double x_max);
void reset_hist_t(hist_t *hist);
void push_hist_t(hist_t *hist, double x);
void silentPush_hist_t(hist_t *hist, double x);
hist_t *deepCopy_hist_t(hist_t *oldHist);

//-- Functions related to RNG
u_int64_t renewSeed();
gsl_rng *initializeGenerator();

//-- Functions related to stopwatch
void printTime(clock_t start, clock_t stop);
void routineTime(clock_t start, clock_t stop);

//-- Math functions
int imod(int N, int i);
int imodc(int N, int i);

//-- Fast math functions
#define DIST_2D_SQ(a,b) (pow((a)[0]-(b)[0], 2) + pow((a)[1]-(b)[1], 2))
#define DIST_3D_SQ(a,b) (pow((a)[0]-(b)[0], 2) + pow((a)[1]-(b)[1], 2) + pow((a)[2]-(b)[2], 2))
#define NORM_2D_SQ(a)   (pow((a)[0], 2) + pow((a)[1], 2))
#define NORM_3D_SQ(a)   (pow((a)[0], 2) + pow((a)[1], 2) + pow((a)[2], 2))
#define SQ(a)           (pow((a), 2))
#define SUM_SQ_2(a,b)   (pow((a), 2) + pow((b), 2))
#define SUM_SQ_3(a,b,c) (pow((a), 2) + pow((b), 2) + pow((c), 2))
#define CB(a)           (pow((a), 3))

//----------------------------------------------------------------------


#endif

