

  /*******************************
   **  commonHeader.c		**
   **  Chieh-An Lin		**
   **  Version 2015.02.25	**
   *******************************/

  
#include "commonHeader.h"


//----------------------------------------------------------------------
//-- Functions related to int_arr

int_arr *initialize_int_arr(int length)
{
  int_arr *iArr = (int_arr*)malloc(sizeof(int_arr));
  iArr->length  = length;
  iArr->array   = (int*)calloc(length, sizeof(int));
  return iArr;
}

void free_int_arr(int_arr *iArr)
{
  if (iArr->array) {free(iArr->array); iArr->array = NULL;}
  free(iArr); iArr = NULL;
  return;
}

void print_int_arr(int_arr *iArr)
{
  printf("# int_arr\n");
  
  int L = (int)fmin(iArr->length, 20);
  int i;
  for (i=0; i<L; i++) printf(" %d ", iArr->array[i]);
  printf("\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to double_arr

double_arr *initialize_double_arr(int length)
{
  double_arr *fArr = (double_arr*)malloc(sizeof(double_arr));
  fArr->length     = length;
  fArr->array      = (double*)calloc(length, sizeof(double));
  return fArr;
}

void free_double_arr(double_arr *fArr)
{
  if (fArr->array) {free(fArr->array); fArr->array = NULL;}
  free(fArr); fArr = NULL;
  return;
}

void print_double_arr(double_arr *fArr)
{
  printf("# double_arr\n");
  
  int L = (int)fmin(fArr->length, 20);
  int i;
  for (i=0; i<L; i++) printf(" %.3f ", fArr->array[i]);
  printf("\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to double_mat

double_mat *initialize_double_mat(int N1, int N2)
{
  double_mat *fMat = (double_mat*)malloc(sizeof(double_mat));
  fMat->N1         = N1;
  fMat->N2         = N2;
  fMat->length     = N1 * N2;
  fMat->matrix     = (double*)calloc(fMat->length, sizeof(double));
  return fMat;
}

void free_double_mat(double_mat *fMat)
{
  if (fMat->matrix) {free(fMat->matrix); fMat->matrix = NULL;}
  free(fMat); fMat = NULL;
  return;
}

void print_double_mat(double_mat *fMat)
{
  printf("# double_mat\n");
  
  int L1 = (int)fmin(fMat->N1, 8);
  int L2 = (int)fmin(fMat->N2, 8);
  int i, j;
  for (j=0; j<L2; j++) {
    for (i=0; i<L1; i++) {
      printf(" %9g ", fMat->matrix[i+j*fMat->N1]);
    }
    printf("\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to double_ten3

double_ten3 *initialize_double_ten3(int N1, int N2, int N3)
{
  double_ten3 *fTen = (double_ten3*)malloc(sizeof(double_ten3));
  fTen->N1          = N1;
  fTen->N2          = N2;
  fTen->N3          = N3;
  fTen->length      = N1 * N2 * N3;
  fTen->tensor      = (double*)calloc(fTen->length, sizeof(double));
  return fTen;
}

void free_double_ten3(double_ten3 *fTen)
{
  if (fTen->tensor) {free(fTen->tensor); fTen->tensor = NULL;}
  free(fTen); fTen = NULL;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to interpolator_t

interpolator_t *initialize_interpolator_t(int length)
{
  //-- length = number of points
  //-- dx     = interval width
  //-- *x     = coordinate
  //-- *value = value
  interpolator_t *inter = (interpolator_t*)malloc(sizeof(interpolator_t));
  inter->length         = length;
  inter->dx             = 0.0;
  inter->x              = (double*)calloc(length, sizeof(double));
  inter->value          = (double*)calloc(length, sizeof(double));
  return inter;
}

void free_interpolator_t(interpolator_t *inter)
{
  if (inter->x)     {free(inter->x);     inter->x = NULL;}
  if (inter->value) {free(inter->value); inter->value = NULL;}
  free(inter); inter = NULL;
  return;
}

void print_interpolator_t(interpolator_t *inter)
{
  printf("# interpolator_t\n");
  printf("# Number of points = %d\n", inter->length);
  printf("# Interval = %g\n", inter->dx);
  printf("#\n");
  printf("#        x      value\n");
  
  int i;
  for (i=0; i<inter->length; i++) printf(" %9g  %9g\n", inter->x[i], inter->value[i]);
  return;
}

double execute_interpolator_t(interpolator_t *inter, double x)
{
  int length    = inter->length;
  double *coor  = inter->x;
  double *value = inter->value;
  
  if (x > coor[length-1] || x < coor[0]) {
    printf("Value out of range of interpolator\n");
    exit(1);
  }
    
  int i;
  for (i=1; i<length; i++) {
    if (coor[i] >= x) break;
  }
  double r = (x - coor[i-1]) / (coor[i] - coor[i-1]); 
  return (1-r)*value[i-1] + r*value[i];
}

//----------------------------------------------------------------------
//-- Functions related to sampler_t

sampler_t *initialize_sampler_t(int length)
{
  //-- length  = number of points
  //-- dx      = interval width
  //-- *x      = bin values
  //-- *pdf    = normalized pdf
  //-- *cdf    = normalized cdf
  //-- totPdf  = integration over pdf, before or after normalization (see set_sampler_t)
  //-- x_mean  = integration over x*p(x), before or after normalization (see set_sampler_t)
  //-- 
  //-- Usage:
  //--   - initialize it with length
  //--   - fill dx, x[0], rest of x
  //--   - fill pdf (not necessarily normalized)
  //--   - use set_sampler_t to set
  
  sampler_t *samp = (sampler_t*)malloc(sizeof(sampler_t));
  samp->length    = length;
  samp->x         = (double*)calloc(length, sizeof(double));
  samp->pdf       = (double*)calloc(length, sizeof(double));
  samp->cdf       = (double*)calloc(length, sizeof(double));
  samp->x_mean    = 0;
  return samp;
}

void free_sampler_t(sampler_t *samp)
{
  if (samp->x)   {free(samp->x);   samp->x   = NULL;}
  if (samp->pdf) {free(samp->pdf); samp->pdf = NULL;}
  if (samp->cdf) {free(samp->cdf); samp->cdf = NULL;}
  free(samp); samp = NULL;
  return;
}

void print_sampler_t(sampler_t *samp)
{
  printf("# sampler_t\n");
  int L = (int)fmin(samp->length, 20);
  int i;
  
  printf("# length = %d\n", samp->length);
  printf("# dx     = %g\n", samp->dx);
  printf("# x      =");
  for (i=0; i<L; i++) printf(" %.3f ", samp->x[i]);
  printf("\n");
  printf("# pdf    =");
  for (i=0; i<L; i++) printf(" %.3f ", samp->pdf[i]);
  printf("\n");
  printf("# cdf    =");
  for (i=0; i<L; i++) printf(" %.3f ", samp->cdf[i]);
  printf("\n");
  printf("# totPdf = %g\n", samp->totPdf);
  printf("# x_mean = %g\n", samp->x_mean);
  return;
}

void set_sampler_t(sampler_t *samp, int setTotalToOne)
{
  //-- Let p(x) = input pdf, q(x) = normalized pdf.
  //-- If setTotalToOne = 0, the input pdf is considered to contain physical information:
  //--   totPdf = integration over p(x)
  //--   x_mean = integration over x*p(x)
  //-- If setTotalToOne = 1, the input pdf is considered to be normalized, so:
  //--   totPdf = integration over q(x) = 1.0
  //--   x_mean = integration over x*q(x)
  //-- In both cases, samp->pdf always represents q(x).
  
  double *x     = samp->x;
  double *pdf   = samp->pdf;
  double *cdf   = samp->cdf;
  int length    = samp->length;
  double totPdf = 0;
  double x_mean = 0;
  int i;
  
  //-- Compute total pdf and <x>
  totPdf += 0.5 * (pdf[0]        + pdf[length-1]);
  x_mean += 0.5 * (pdf[0] * x[0] + pdf[length-1] * x[length-1]);
  for (i=1; i<=length-2; i++) {
    totPdf += pdf[i];
    x_mean += pdf[i] * x[i];
  }
  
  //-- Normalization
  for (i=0; i<length; i++) pdf[i] /= totPdf;
  samp->totPdf = totPdf * samp->dx; //-- Normalization by bin width
  samp->x_mean = x_mean * samp->dx; //-- Normalization by bin width
  if (setTotalToOne) {
    samp->totPdf  = 1.0;
    samp->x_mean /= samp->totPdf;
  }
  
  //-- Fill cdf
  cdf[0] = 0;
  for (i=1; i<length; i++) cdf[i] = cdf[i-1] + 0.5 * (pdf[i-1] + pdf[i]);
  return;
}

double execute_sampler_t(sampler_t *samp, double x)
{
  int length    = samp->length;
  double *cdf   = samp->cdf;
  double *value = samp->x;
    
  int i;
  for (i=1; i<length; i++) {
    if (cdf[i] >= x) break;
  }
  double r = (x - cdf[i-1]) / (cdf[i] - cdf[i-1]); 
  return (1-r)*value[i-1] + r*value[i];
}

//----------------------------------------------------------------------
//-- Functions related to sampler_arr

sampler_arr *initialize_sampler_arr(int length, int nbPoints)
{
  sampler_arr *sampArr  = (sampler_arr*)malloc(sizeof(sampler_arr));
  sampArr->length       = length;
  sampArr->array        = (sampler_t**)malloc(length * sizeof(sampler_t*));
  int i;
  for (i=0; i<length; i++) sampArr->array[i] = initialize_sampler_t(nbPoints);   
  return sampArr;
}

void free_sampler_arr(sampler_arr *sampArr)
{
  int i;
  if (sampArr->array) {
    for (i=0; i<sampArr->length; i++) {free_sampler_t(sampArr->array[i]); sampArr->array[i] = NULL;}
  }
  free(sampArr); sampArr = NULL;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to FFT_t

FFT_t *initialize_FFT_t(int N1, int N2)
{
  FFT_t *transformer    = (FFT_t*)malloc(sizeof(FFT_t));
  transformer->N1       = N1;
  transformer->N2       = N2;
  transformer->length   = N1 * N2;
  transformer->length_f = N1 * (N2/2 + 1);
  
  //-- Allocate direct space elements
  transformer->before   = (double*)calloc(transformer->length, sizeof(double));
  transformer->kernel   = (double*)calloc(transformer->length, sizeof(double));
  transformer->after    = (double*)calloc(2 * transformer->length_f, sizeof(double));
  
  //-- Allocate Fourier space elements
  transformer->before_f = (fftw_complex*)fftw_malloc(transformer->length_f * sizeof(fftw_complex));
  transformer->kernel_f = (fftw_complex*)fftw_malloc(transformer->length_f * sizeof(fftw_complex));
  transformer->after_f  = (fftw_complex*)fftw_malloc(transformer->length_f * sizeof(fftw_complex));
  
  //-- Allocate fftw_plan elements
  transformer->before_p = fftw_plan_dft_r2c_2d(N1, N2, transformer->before,  transformer->before_f, FFTW_ESTIMATE);
  transformer->kernel_p = fftw_plan_dft_r2c_2d(N1, N2, transformer->kernel,  transformer->kernel_f, FFTW_ESTIMATE);
  transformer->after_p  = fftw_plan_dft_c2r_2d(N1, N2, transformer->after_f, transformer->after,    FFTW_ESTIMATE);
  
  return transformer;
}

void free_FFT_t(FFT_t *transformer)
{
  if (transformer->before)   {free(transformer->before);                transformer->before = NULL;}
  if (transformer->kernel)   {free(transformer->kernel);                transformer->kernel = NULL;}
  if (transformer->after)    {free(transformer->after);                 transformer->after = NULL;}
  
  if (transformer->before_f) {fftw_free(transformer->before_f);         transformer->before_f = NULL;}
  if (transformer->kernel_f) {fftw_free(transformer->kernel_f);         transformer->kernel_f = NULL;}
  if (transformer->after_f)  {fftw_free(transformer->after_f);          transformer->after_f = NULL;}
  
  if (transformer->before_p) {fftw_destroy_plan(transformer->before_p); transformer->before_p = NULL;}
  if (transformer->kernel_p) {fftw_destroy_plan(transformer->kernel_p); transformer->kernel_p = NULL;}
  if (transformer->after_p)  {fftw_destroy_plan(transformer->after_p);  transformer->after_p = NULL;}
  
  free(transformer); transformer = NULL;
  return;
}

void reset_FFT_t(FFT_t *transformer)
{ 
  int i;
  for (i=0; i<transformer->length; i++) transformer->before[i] = 0.0;
  return;
}

void execute_FFT_t(FFT_t *transformer)
{ 
  fftw_execute(transformer->before_p);
  
  //-- Complex coefficient multiplication in Fourier space
  int i;
  for (i=0; i<transformer->length_f; i++) {
    transformer->after_f[i][0] = transformer->before_f[i][0] * transformer->kernel_f[i][0] - transformer->before_f[i][1] * transformer->kernel_f[i][1];
    transformer->after_f[i][1] = transformer->before_f[i][0] * transformer->kernel_f[i][1] + transformer->before_f[i][1] * transformer->kernel_f[i][0];
  }
  
  fftw_execute(transformer->after_p);
  
  //-- Rescale
  for (i=0; i<transformer->length; i++) transformer->after[i] /= (double)transformer->length;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to hist_t

hist_t *initialize_hist_t(int length)
{
  //-- length   = number of bins
  //-- dx       = bin width
  //-- n_tot    = total counts
  //-- *x_lower = lower limit of each bin
  //-- x_max;   = maximal limit
  //-- *n       = histogram
  hist_t *hist  = (hist_t*)malloc(sizeof(hist_t));
  hist->length  = length;
  hist->dx      = 0.0;
  hist->x_lower = (double*)calloc(length, sizeof(double));
  hist->x_max   = 0.0;
  hist->n       = (int*)calloc(length, sizeof(int));
  hist->n_tot   = 0;
  return hist;
}

void free_hist_t(hist_t *hist)
{
  if (hist->x_lower) {free(hist->x_lower); hist->x_lower = NULL;}
  if (hist->n)       {free(hist->n);       hist->n = NULL;}
  free(hist); hist = NULL;
  return;
}

void print_hist_t(hist_t *hist)
{
  printf("# hist_t\n");
  printf("# Number of bins = %d\n", hist->length);
  printf("# Bin width = %g\n", hist->dx);
  printf("# Total counts = %d\n", hist->n_tot);
  printf("#\n");
  printf("#                   bin        n\n");
  
  int i;
  for (i=0; i<hist->length; i++) printf("  [%9g, %9g]  %6d\n", hist->x_lower[i], hist->x_lower[i]+hist->dx, hist->n[i]);
  return;
}

void set_hist_t(hist_t *hist, double x_min, double x_max)
{
  int length  = hist->length;
  double dx   = (x_max - x_min) / (double)length;
  hist->dx    = dx;
  hist->n_tot = 0;
  int i;
  for (i=0; i<length; i++) {
    hist->x_lower[i] = x_min + i * dx;
    hist->n[i]       = 0;
  }
  hist->x_max = hist->x_lower[length-1] + dx;
  return;
}

void reset_hist_t(hist_t *hist)
{
  int i;
  for (i=0; i<hist->length; i++) hist->n[i] = 0;
  return;
}

void push_hist_t(hist_t *hist, double x)
{
  double *x_lower = hist->x_lower;
  if (x < x_lower[0] || x >= hist->x_max) {
    printf("Value out of range, not pushed into histogram\n");
    return;
  }
  
  int i;
  for (i=0; i<hist->length-1; i++) {
    if (x_lower[i+1] > x) break;
  }
  hist->n[i]  += 1;
  hist->n_tot += 1;
  return;
}

void silentPush_hist_t(hist_t *hist, double x)
{
  double *x_lower = hist->x_lower;
  if (x < x_lower[0])   return;
  if (x >= hist->x_max) return;
  
  int i;
  for (i=0; i<hist->length-1; i++) {
    if (x_lower[i+1] > x) break;
  }
  hist->n[i]  += 1;
  hist->n_tot += 1;
  return;
}

hist_t *deepCopy_hist_t(hist_t *oldHist)
{
  hist_t *newHist = initialize_hist_t(oldHist->length);
  newHist->dx     = oldHist->dx;
  newHist->n_tot  = oldHist->n_tot;
  newHist->x_max  = oldHist->x_max;
  int i;
  for (i=0; i<oldHist->length; i++) {
    newHist->x_lower[i] = oldHist->x_lower[i];
    newHist->n[i]       = oldHist->n[i];
  }
  return newHist;
}

//----------------------------------------------------------------------
//-- Functions related to RNG

u_int64_t renewSeed()
{
  FILE *file = fopen("/dev/urandom","r");
  u_int64_t seed;
  if (file == NULL) seed = 0;
  else {
    int state = fread(&seed, sizeof(u_int64_t), 1, file);
    fclose(file);
  }
  return seed;
}

gsl_rng *initializeGenerator()
{
  u_int64_t seed = renewSeed();
  gsl_rng *generator = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(generator, seed);
  return generator;
}

//-- Use the following functions:
//--   double gsl_ran_flat(const gsl_rng *r, double a, double b)
//--   double gsl_ran_gaussian(const gsl_rng *r, double sigma)

//----------------------------------------------------------------------
//-- Functions related to stopwatch

//-- Use the following functions:
//--   clock_t start = clock();
//--   clock_t finish = clock();
//--   double duration = (double)(finish - start) / CLOCKS_PER_SEC;

void printTime(clock_t start, clock_t stop)
{
  double duration = (double)(stop - start) / CLOCKS_PER_SEC;
  int hours       = (int)(duration / 3600);
  double remain   = fmod(duration, 3600);
  int minutes     = (int)(remain / 60);
  double seconds  = fmod(remain, 60);
  
  if (hours == 0) {
    if (minutes == 0) printf("Computation time = %.2f secs\n", duration);
    else              printf("Computation time = %d m %d s\n", minutes, (int)seconds);
  }
  else                printf("Computation time = %d h %d m %d s\n", hours, minutes, (int)seconds);
  return;
}

void routineTime(clock_t start, clock_t stop)
{
  double duration = (double)(stop - start) / CLOCKS_PER_SEC;
  int hours       = (int)(duration / 3600);
  double remain   = fmod(duration, 3600);
  int minutes     = (int)(remain / 60);
  double seconds  = fmod(remain, 60);
  
  if (hours == 0) {
    if (minutes == 0) printf("routine time = %.2f secs\n", duration);
    else              printf("routine time = %d m %d s\n", minutes, (int)seconds);
  }
  else                printf("routine time = %d h %d m %d s\n", hours, minutes, (int)seconds);
  return;
}

//----------------------------------------------------------------------
//-- Math functions

int imod(int N, int i)
{ 
  //-- mod function ranged in [0, N-1]
  if (i < 0)  return imod(N, i + N);
  if (i >= N) return imod(N, i - N);
  return i;
}

int imodc(int N, int i)
{ 
  //-- mod function centered on 0, positive side is priviledged if N is even
  if (i <= -N + N/2) return imodc(N, i + N);
  if (i > N/2)       return imodc(N, i - N);
  return i;
}

//----------------------------------------------------------------------
