

  /*******************************
   **  smoothing.c		**
   **  Chieh-An Lin		**
   **  Version 2015.12.08	**
   *******************************/


#include "smoothing.h"


//----------------------------------------------------------------------
//-- Functions related to map_t

map_t *initialize_map_t(int N1, int N2, double theta_pix, error **err)
{
  map_t *kMap     = (map_t*)malloc_err(sizeof(map_t), err);                 forwardError(*err, __LINE__,);
  kMap->N1        = N1;
  kMap->N2        = N2;
  kMap->length    = N1 * N2;
  kMap->theta_pix = theta_pix;
  kMap->theta_pix_inv = 1.0 / theta_pix;
  
  kMap->limits[0] = 0;
  kMap->limits[1] = N1 * theta_pix;
  kMap->limits[2] = 0;
  kMap->limits[3] = N2 * theta_pix;
  kMap->center[0] = 0.5 * kMap->limits[1];
  kMap->center[1] = 0.5 * kMap->limits[3];
  
  kMap->type      = kappa_map;
  kMap->value1    = (double*)malloc_err(kMap->length * sizeof(double), err); forwardError(*err, __LINE__,);
  kMap->value2    = (double*)malloc_err(kMap->length * sizeof(double), err); forwardError(*err, __LINE__,);
  return kMap;
}

void free_map_t(map_t *kMap)
{
  if (kMap->value1) {free(kMap->value1); kMap->value1 = NULL;}
  if (kMap->value2) {free(kMap->value2); kMap->value2 = NULL;}
  free(kMap); kMap = NULL;
  return;
}

void getPixPos_map_t(map_t *kMap, double pos[2], int i, int j)
{
  //-- Compute the position of the pixel indexed (i, j) and stock this information in pos.
  pos[0] = (0.5 + i) * kMap->theta_pix;
  pos[1] = (0.5 + j) * kMap->theta_pix;
  return;
}

void getPixPos(double pos[2], double limits[4], double theta_pix, int i, int j)
{
  //-- Compute the position of the pixel indexed (i, j) and stock this information in pos.
  pos[0] = (0.5 + i) * theta_pix;
  pos[1] = (0.5 + j) * theta_pix;
  return;
}

void read_map_t(char *name, peak_param *peak, map_t *kMap, mapType_t type, error **err)
{
  testErrorRet(kMap->length!=peak->resol[0]*peak->resol[1], peak_match, "Resolution match error", *err, __LINE__,);
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double *value1 = kMap->value1;
  double *value2 = kMap->value2;
  
  //-- Read
  int c = fgetc(file);
  if (type < 6) {
    while (c != EOF) {
      if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
      else {
	testErrorRet(count>=kMap->length, peak_overflow, "Too many pixels", *err, __LINE__,);
	ungetc(c, file);
	buffer2 = fscanf(file, "%*f %*f %lf\n", &value1[count]);
	count++;
      }
      c = fgetc(file);
    }
  }
  else {
    while (c != EOF) {
      if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
      else {
	testErrorRet(count>=kMap->length, peak_overflow, "Too many pixels", *err, __LINE__,);
	ungetc(c, file);
	buffer2 = fscanf(file, "%*f %*f %lf %lf\n", &value1[count], &value2[count]);
	count++;
      }
      c = fgetc(file);
    }
  }
  fclose(file);
  
  //-- Length check
  testErrorRet(count!=kMap->length, peak_match, "Array length error", *err, __LINE__,);
  kMap->type = type;
  printf("\"%s\" read\n", name);
  return;
}

void output_map_t(FILE *file, map_t *kMap)
{
  mapType_t type = kMap->type;
  
  fprintf(file, "# Number of pixels = %d\n", kMap->length);
  fprintf(file, "#\n");
  
  double pos[2];
  int jN1;
  
  int i, j;
  if (type < 6 || type == mask_map) {
    fprintf(file, "#  theta_x   theta_y    value1\n");
    fprintf(file, "# [arcmin]  [arcmin]       [-]\n");
    
    for (j=0; j<kMap->N2; j++) {
      jN1 = j * kMap->N1;
      for (i=0; i<kMap->N1; i++) {
	getPixPos_map_t(kMap, pos, i, j);
	fprintf(file, "  %8.3f  %8.3f  %8.5f\n", pos[0], pos[1], kMap->value1[i+jN1]);
      }
    }
  }
    
  else {
    fprintf(file, "#  theta_x   theta_y    value1    value2\n");
    fprintf(file, "# [arcmin]  [arcmin]       [-]       [-]\n");
    
    int index;
    for (j=0; j<kMap->N2; j++) {
      jN1 = j * kMap->N1;
      for (i=0; i<kMap->N1; i++) {
	index = i + jN1;
	getPixPos_map_t(kMap, pos, i, j);
	fprintf(file, "  %8.3f  %8.3f  %8.5f  %8.5f\n", pos[0], pos[1], kMap->value1[index], kMap->value2[index]);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to smoothing by FFT

void fillGaussianKernel(fftw_complex *kernel, int length, int M, double scaleInPix)
{
  //-- A kernel should be placed in this way
  //--  1.0  0.6   0   0   0  0.6
  //--  0.6  0.1   0   0   0  0.1
  //--    0    0   0   0   0    0
  //--    0    0   0   0   0    0
  //--    0    0   0   0   0    0
  //--  0.6  0.1   0   0   0  0.1
  //-- 
  //-- Gaussian kernel:
  //-- exp[-(x^2 + y^2)/s^2]
  //-- where s is theta_G / theta_pix
  //--
  //-- M is resolution + size of zero padding
  
  //-- Initialization
  reset_fftw_complex(kernel, length);
  
  double scaleInPix_invSq = 1.0 / (SQ(scaleInPix)); //-- [pix^2]
  double relay = CUTOFF_FACTOR_FILTER * scaleInPix; //-- CUTOFF_FACTOR_FILTER set to 2.2 in peakParameters.h, equivalent to a Gaussian trunked at 3-sigma
  int cutInPix = (int)ceil(relay);
  double cutInPix_sq = SQ(relay);
  double sum   = 0.0;
  double value;
  int i, j, mod_j_M, r_sq;
  
  for (j=-cutInPix; j<=cutInPix; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    for (i=-cutInPix; i<=cutInPix; i++) {
      r_sq = i*i + j*j;
      if (r_sq > cutInPix_sq) continue;
      value = exp(-r_sq * scaleInPix_invSq);
      sum += value;
      kernel[POS_MOD(M, i) + mod_j_M][0] = value;
    }
  }
  
  //-- Normalization
  rescaleReal_fftw_complex(kernel, length, 1.0/sum);
  return;
}

double starlet_2D(double x, double y)
{
  //-- This function computes 144 * psi(x, y).
  //--
  //-- Scaling function:
  //--   phi(t)  = 1/12 * (|t - 2|^3  - 4 |t - 1|^3   + 6 |t|^3  - 4 |t + 1|^3   + |t + 2|^3)
  //--   phi(2t) = 1/12 * (|2t - 2|^3 - 4 |2t - 1|^3  + 6 |2t|^3 - 4 |2t + 1|^3  + |2t + 2|^3)
  //--           = 8/12 * (|t - 1|^3  - 4 |t - 1/2|^3 + 6 |t|^3  - 4 |t + 1/2|^3 + |t + 1|^3)
  //--
  //-- 2D-starlet:
  //--   psi(x, y) = 4 phi(2x) phi(2y) - phi(x) phi(y)
  //--             = 1/144 [4 * 8 (|x - 1|^3 - 4 |x - 1/2|^3 + 6 |x|^3 - 4 |x + 1/2|^3 + |x + 1|^3) * 8 (|y - 1|^3 - 4 |y - 1/2|^3 + 6 |y|^3 - 4 |y + 1/2|^3 + |y + 1|^3)
  //--                      -     (|x - 2|^3 - 4 |x - 1|^3   + 6 |x|^3 - 4 |x + 1|^3   + |x + 2|^3) *   (|y - 2|^3 - 4 |y - 1|^3   + 6 |y|^3 - 4 |y + 1|^3   + |y + 2|^3)]
  //--             = 1/144 [256 (x_one - 4 x_half + 6 x_zero) (y_one - 4 y_half + 6 y_zero) - (x_two - 4 x_one  + 6 x_zero) (y_two - 4 y_one  + 6 y_zero)]
  //--   where
  //--     x_zero = |x|^3,  x_half = |x - 1/2|^3 + |x + 1/2|^3,  x_one  = |x - 1|^3 + |x + 1|^3,  x_two  = |x - 2|^3 + |x + 2|^3,
  //--     y_zero = |y|^3,  y_half = |y - 1/2|^3 + |y + 1/2|^3,  y_one  = |y - 1|^3 + |y + 1|^3,  y_two  = |y - 2|^3 + |y + 2|^3
  
  double x_zero    = fabs(pow(x,       3.0));
  double x_half    = fabs(pow(x - 0.5, 3.0)) + fabs(pow(x + 0.5, 3.0));
  double x_one     = fabs(pow(x - 1.0, 3.0)) + fabs(pow(x + 1.0, 3.0));
  double x_two     = fabs(pow(x - 2.0, 3.0)) + fabs(pow(x + 2.0, 3.0));
  double y_zero    = fabs(pow(y,       3.0));
  double y_half    = fabs(pow(y - 0.5, 3.0)) + fabs(pow(y + 0.5, 3.0));
  double y_one     = fabs(pow(y - 1.0, 3.0)) + fabs(pow(y + 1.0, 3.0));
  double y_two     = fabs(pow(y - 2.0, 3.0)) + fabs(pow(y + 2.0, 3.0));
  
  double value = 0.0;
  //-- This is not yet normalized.
  value += 256.0 * (x_one - 4*x_half + 6*x_zero) * (y_one - 4*y_half + 6*y_zero);
  value -=         (x_two - 4*x_one  + 6*x_zero) * (y_two - 4*y_one  + 6*y_zero);
  return value;
}

void fillStarletKernel(fftw_complex *kernel, int length, int M, double scaleInPix)
{
  //-- L2-norm of the 2D-starlet:
  //--   || psi ||_L2^2 = STARLET_NORM^2 = 5 I2^2 - 2 I3^2
  //--   where
  //--     I2 = 2/5 + 5/63
  //--     I3 = 1/3 + 1/5 + 1/21 + 1/48
  //--
  //-- Scaled 2D-starlet:
  //--   Let x = i/s, y = j/s, define psi_s(i, j) = psi(i/s, j/s) / (STARLET_NORM1 s^2).
  //--   Then, || psi_s ||_L1   = 1,
  //--   and   || psi_s ||_L2^2 = STARLET_NORM2^2 / (STARLET_NORM1^2 s^2)
  //--
  //--   In this case, S/N, kappa divided by sigma_noise, is given by
  //--   sigma_noise^2 = sigma_eps^2 / (2 n_g s^2) * (STARLET_NORM2 / STARLET_NORM1)^2.
  //--
  //-- Implementation:
  //--   1. Compute 144 * psi(i/s, j/s).
  //--   2. Normalize the kernel by the sum of absolute coefficients.
  //--   This cancels the factor 144 and set L1-norm = 1.
  
  //-- Initialization
  reset_fftw_complex(kernel, length);
  
  int cutInPix = MIN((int)ceil(2 * scaleInPix), M/2); //-- [pix]
  double sum   = 0.0;
  double x, y, value;
  int i, j, mod_j_M, r_sq;
  
  for (j=-cutInPix; j<=cutInPix; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    y = (double)j / scaleInPix;
    for (i=-cutInPix; i<=cutInPix; i++) {
      x = (double)i / scaleInPix;
      value = starlet_2D(x, y);
      sum  += fabs(value); //-- Compute the sum of the absolute value of coefficients
      kernel[POS_MOD(M, i) + mod_j_M][0] = value;
    }
  }
  
  //-- Normalization
  rescaleReal_fftw_complex(kernel, length, 1.0/sum);
  return;
}

void makeKernel(peak_param *peak, FFT_arr *smoother)
{
  FFT_t *smoo;
  int i;
  for (i=0; i<peak->nbLinFilters; i++) {
    smoo = smoother->array[i];
    if (peak->linFilter[i] == gauss)     fillGaussianKernel(smoo->kernel, smoo->length, smoo->N, peak->linScaleInPix[i]);
    else if (peak->linFilter[i] == star) fillStarletKernel(smoo->kernel, smoo->length, smoo->N, peak->linScaleInPix[i]);
    fftw_execute(smoo->kernel_f); //-- To Fourier space
  }
  return;
}

void smoothByFFT_arr(peak_param *peak, gal_map *gMap, FFT_arr *smoother)
{
  execute_FFT_arr(smoother);
  gMap->type += 1;
  if (peak->printMode < 2) printf("FFT smoothing done\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to smoothing and binning by direct convolution

void DCWithGaussianForPair(gal_list *gList, FFT_t **smooArr, const double *cut_sq, const double *scale_invSq, int nbLinScales, int doKappa, 
			   double *weight, double pixPos[2], int index_FFT)
{
  gal_node *gNode;
  gal_t *g;
  double *pix, dist_sq, factor;
  int i, j;
  
  if (doKappa == 1) {
    for (i=0, gNode=gList->first; i<gList->size; i++, gNode=gNode->next) {
      g = gNode->g;
      dist_sq = DIST_2D_SQ(pixPos, g->pos);
      for (j=0; j<nbLinScales; j++) {
	if (dist_sq > cut_sq[j]) continue;
	factor     = exp(-dist_sq * scale_invSq[j]);
	weight[j] += factor;
	pix        = smooArr[j]->before[index_FFT];
	pix[0]    += factor * g->kappa;
      }
    }
  }
  
  else {
    for (i=0, gNode=gList->first; i<gList->size; i++, gNode=gNode->next) {
      g = gNode->g;
      dist_sq = DIST_2D_SQ(pixPos, g->pos);
      for (j=0; j<nbLinScales; j++) {
	if (dist_sq > cut_sq[j]) continue;
	factor     = exp(-dist_sq * scale_invSq[j]);
	weight[j] += factor;
	pix        = smooArr[j]->before[index_FFT];
	pix[0]    += factor * g->gamma[0];
	pix[1]    += factor * g->gamma[1];
      }
    }
  }
  return;
}

void DCWithGaussianForPixel(gal_map *gMap, FFT_arr *smoother, const double *cut_sq, const double *scale_invSq, int nbLinScales, int doKappa, 
			    int M, int bufferSize, double theta_pix, double *weight, int i_pix, int j_pix)
{
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int index_FFT = i_pix + j_pix * M;
  gal_list **map = gMap->map;
  double pixPos[2];
  int i, j, jN1;
  
  pixPos[0] = (0.5 + i_pix) * theta_pix;
  pixPos[1] = (0.5 + j_pix) * theta_pix;
  reset_double(weight, nbLinScales);
  
  for (j=j_pix-bufferSize; j<=j_pix+bufferSize; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_pix-bufferSize; i<=i_pix+bufferSize; i++) {
      if (i < 0 || i >= N1) continue;
      DCWithGaussianForPair(map[i+jN1], smoother->array, cut_sq, scale_invSq, nbLinScales, doKappa, weight, pixPos, index_FFT);
    }
  }
  
  double *pix, factor;
  
  //-- Normalization
  for (i=0; i<nbLinScales; i++) {
    if (weight[i] == 0) continue;
    factor  = 1.0 / weight[i];
    pix     = smoother->array[i]->before[index_FFT];
    pix[0] *= factor;
    pix[1] *= factor;
  }
  return;
}

void DCWithGaussian(peak_param *peak, gal_map *gMap, FFT_arr *smoother)
{
  int N1          = gMap->N1;
  int N2          = gMap->N2;
  int M           = peak->FFTSize;
  int nbLinScales = smoother->length;
  int FFTLength   = smoother->array[0]->length;
  int i, j, k, jN1, jM, index_gMap;
  
  //-- Reset all tables
  for (k=0; k<nbLinScales; k++) reset_fftw_complex(smoother->array[k]->before, FFTLength);
  
  //-- Direct convolution
  for (j=0; j<N2; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<N1; i++) {
      index_gMap = i + jN1;
      if ((peak->printMode == 0) && (index_gMap % 50 == 0)) {
	printf("Computing direct convolution: %6.2f%% \r", 100.0 * index_gMap / (double)gMap->length);
	fflush(stdout);
      }
      DCWithGaussianForPixel(gMap, smoother, peak->cut_sq, peak->scale_invSq, nbLinScales, peak->doKappa, 
			     M, peak->bufferSize, peak->theta_pix, peak->weight, i, j);
    }
  }
  
  if (peak->printMode < 2) printf("Direct convolution done                 \n");
  return;
}

void smoothAndBinByDC(peak_param *peak, gal_map *gMap, FFT_arr *smoother)
{
  //TODO Mask needs to be update
  int Gaussian = 1;
  if (Gaussian) DCWithGaussian(peak, gMap, smoother);
  gMap->type += 1;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to binning and noise

void galaxyBinning(gal_map *gMap, FFT_t *smoo)
{
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int M  = smoo->N;
  gal_list **map            = gMap->map;
  fftw_complex *smoo_before = smoo->before;
  double *mean;
  int i, j, jN1, jM, index_FFT;
    
  for (j=0; j<M; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<M; i++) {
      index_FFT = i + jM; 
      if (i >= N1 || j >= N2) { //-- Buffer area
	smoo_before[index_FFT][0] = 0.0;
	smoo_before[index_FFT][1] = 0.0;
      }
      else {
	mean = map[i+jN1]->mean;
	smoo_before[index_FFT][0] = mean[0];
	smoo_before[index_FFT][1] = mean[1];
      }
    }
  }
  return;
}

void copyToAllBefore(FFT_arr *smoother)
{
  int i;
  for (i=1; i<smoother->length; i++) copy_fftw_complex(smoother->array[0]->before, smoother->array[i]->before, smoother->array[0]->length);
  return;
}

void addNoiseToTable(peak_param *peak, gal_map *gMap, fftw_complex *table)
{
  //-- Only used in paperIII.c
  
  int N1 = peak->resol[0];
  int N2 = peak->resol[1];
  int M  = peak->FFTSize;
  double sigma_pix   = peak->sigma_pix;
  gsl_rng *generator = peak->generator;
  
  int i, j, jM;
  for (j=0; j<M; j++) {
    jM  = j * M;
    for (i=0; i<M; i++) {
      if (i >= N1 || j >= N2) continue; //-- Buffer area
      table[i+jM][0] += gsl_ran_gaussian(generator, sigma_pix);
    }
  }
  
  gMap->type += 2;
  if (peak->printMode < 2) printf("Added noise to pixels\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to linear KS inversion

void linearKS(fftw_complex *table, int M)
{
  double k1, k2, k1_sq, k2_sq, k_norm_inv, f1, f2;
  double *pix, value;
  int jM;
  
  int i, j;
  for (j=0; j<M; j++) {
    k2 = CEN_MOD(M, j);
    jM = j * M;
    
    for (i=0; i<M; i++) {
      if (i + j == 0) continue; //-- This means if i = 0 and j = 0
      pix        = table[i+jM];
      k1         = CEN_MOD(M, i);
      k1_sq      = SQ(k1);
      k2_sq      = SQ(k2);
      k_norm_inv = 1.0 / (k1_sq + k2_sq);
      f1         = (k1_sq - k2_sq) * k_norm_inv;
      f2         = 2 * k1 * k2 * k_norm_inv;
      value      =  f1 * pix[0] + f2 * pix[1];
      pix[1]     = -f2 * pix[0] + f1 * pix[1];
      pix[0]     = value;
    }
  }
  
  table[0][0] = 0.0;
  table[0][1] = 0.0;
  return;
}

void zeroPadding(fftw_complex *table, int N1, int N2, int M)
{
  int i, j, jM, index_FFT;
  for (j=0; j<M; j++) {
    jM  = j * M;
    for (i=0; i<M; i++) {
      index_FFT = i + jM;
      if (i >= N1 || j >= N2) table[index_FFT][0] = 0.0; //-- Buffer area
      table[index_FFT][1] = 0.0;
    }
  }
  return;
}

void invertByLinKS(peak_param *peak, FFT_t *smoo)
{
  int N1 = peak->resol[0];
  int N2 = peak->resol[1];
  int M  = smoo->N;
  fftw_complex *smoo_before = smoo->before;
    
  fftw_execute(smoo->before_f);                                             //-- Go to Fourier space
  linearKS(smoo_before, M);                                                 //-- KS inversion
  fftw_execute(smoo->before_b);                                             //-- Go to direct space
  rescaleReal_fftw_complex(smoo_before, smoo->length, peak->FFTNormFactor); //-- Rescale, only real part is interesting.
  if (peak->doSmoothing != 2) zeroPadding(smoo_before, N1, N2, M);          //-- 0-padding
  return;
}

void invertByLinKS_arr(peak_param *peak, gal_map *gMap, FFT_arr *smoother)
{
  invertByLinKS(peak, smoother->array[0]);
  int i;
  if (peak->doSmoothing == 2) {for (i=1; i<smoother->length; i++) invertByLinKS(peak, smoother->array[i]);}
  else                         copyToAllBefore(smoother);
  
  gMap->type -= 6;
  if (peak->printMode < 2) printf("Linear KS done\n");
  return;
}

void linKSAndFFT(peak_param *peak, FFT_t *smoo)
{
  int FFTLength = smoo->length;
  fftw_complex *smoo_before = smoo->before;
  fftw_complex *smoo_after  = smoo->after;
  
  //-- Linear KS
  fftw_execute(smoo->before_f);   //-- Go to Fourier space
  linearKS(smoo_before, smoo->N); //-- KS inversion
  
  //-- FFT smoothing
  multiplication_fftw_complex(smoo_before, smoo->kernel, smoo_after, FFTLength); //-- Multiplication
  fftw_execute(smoo->after_b);                                                   //-- Go to direct space
  rescaleReal_fftw_complex(smoo_after, FFTLength, peak->FFTNormFactor);          //-- Rescale, only real part is interesting.
  return;
}

void linKSAndFFT_arr(peak_param *peak, FFT_arr *smoother)
{
  //-- Not used
  
  FFT_t *smoo = smoother->array[0];
  linKSAndFFT(peak, smoo);
  
  fftw_complex *smoo_before = smoo->before;
  int FFTLength = smoo->length;
  int i;
  
  for (i=1; i<smoother->length; i++) {
    smoo = smoother->array[i];
    multiplication_fftw_complex(smoo_before, smoo->kernel, smoo->after, FFTLength); //-- Multiplication
    fftw_execute(smoo->after_b);                                                    //-- Go to direct space
    rescaleReal_fftw_complex(smoo->after, FFTLength, peak->FFTNormFactor);          //-- Rescale, only real part is interesting.
  }
  if (peak->printMode < 2) printf("Linear KS and FFT smoothing done\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to iterative KS inversion

#define EPS 1e-5
void iterativeKS(fftw_complex *reducedShear, fftw_complex *table, fftw_plan forward, fftw_plan backward, map_t *kMap, double FFTNormFactor, int M)
{
  //-- At output, kappa is stocked in kMap->value1 and table.
  
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  int FFTLength = M * M;
  double *copy = kMap->value1;
  
  reset_double(copy, kMap->length);
  //reset_fftw_complex(table, FFTLength);
  
  double factor, diff;
  int i, j, jN1, jM, index_kMap, index_FFT;
  
  //-- Iterations
  do {
    //-- Update gamma to table
    for (j=0; j<M; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<M; i++) {
	index_FFT = i + jM;
	if (i >= N1 || j >= N2) {
	  table[index_FFT][0] = 0.0;
	  table[index_FFT][1] = 0.0;
	  continue;
	}
	factor = 1.0 - MIN(0.99, copy[i+jN1]); //-- 1 - kappa, with regularization
	table[index_FFT][0] = factor * reducedShear[index_FFT][0]; //-- gamma_1
	table[index_FFT][1] = factor * reducedShear[index_FFT][1]; //-- gamma_2
      }
    }
    fftw_execute(forward);                                     //-- Go to Fourier space
    linearKS(table, M);                                        //-- KS inversion
    fftw_execute(backward);                                    //-- Go to direct space
    rescaleReal_fftw_complex(table, FFTLength, FFTNormFactor); //-- Rescale, only Re(kappa) is interesting.
    
    diff = 0.0;
    for (j=0; j<N2; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<N1; i++) {
	index_FFT = i + jM;
	index_kMap = i + jN1;
	diff = MAX(diff, fabs(table[index_FFT][0] - copy[index_kMap])); //-- Determine norm_infty
	copy[index_kMap] = table[index_FFT][0];                         //-- Update kappa
      }
    }
    
  } while (diff > EPS);
  
  return;
}
#undef EPS

void invertByIterKS(peak_param *peak, FFT_t *smoo, map_t *kMap)
{
  int M = smoo->N;
  fftw_complex *smoo_before = smoo->before;
  fftw_complex *smoo_after  = smoo->after;
  
  copy_fftw_complex(smoo_before, smoo_after, smoo->length);                                           //-- Move reduced shear to smoo->after
  iterativeKS(smoo_after, smoo_before, smoo->before_f, smoo->before_b, kMap, peak->FFTNormFactor, M); //-- Iterative KS
  if (peak->doSmoothing != 2) zeroPadding(smoo_before, peak->resol[0], peak->resol[1], M);            //-- 0-padding
  return;
}

void invertByIterKS_arr(peak_param *peak, gal_map *gMap, FFT_arr *smoother, map_t *kMap)
{
  invertByIterKS(peak, smoother->array[0], kMap);
  int i;
  if (peak->doSmoothing == 2) {for (i=1; i<smoother->length; i++) invertByIterKS(peak, smoother->array[i], kMap);}
  else                         copyToAllBefore(smoother);
  
  gMap->type -= 12;
  if (peak->printMode < 2) printf("Iterative KS done\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to nonlinear filtering

void MRLensFiltering(peak_param *peak, char executable[], char option[], char input[], char output[], char log[])
{
  //-- This function calls an external executable wl_t2_filter, takes input, and returns the result as output.
  //-- Example:
  //--   executable = "../../mrlens/wl_t2_filter"
  //--   option     = "-k -C2 -c2 -s0.05 -F1 -n5"
  //--   input      = "kappaMap_unsmoothed.fits"
  //--   output     = "kappaMap_mrlens.fits"
  
  char command[STRING_LENGTH_MAX];
  sprintf(command, "%s %s %s %s %s", executable, option, input, output, log);
  int buffer = system(command);
  if      (peak->printMode < 2)  printf("\"%s\" made\n", output);
  else if (peak->printMode == 2) printf("MRLens done, ");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to output

void outputMap(char name[], cosmo_hm *cmhm, peak_param *peak, map_t *kMap)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# %s\n", STR_MAPTYPE_T(kMap->type));
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file, "# Filter = %s, scale = %g [arcmin] = %g [pix]\n", STR_FILTER_T(peak->filter[0]), peak->scale[0], peak->linScaleInPix[0]);
  fprintf(file, "# sigma_eps = %g, sigma_pix = %g\n", peak->sigma_eps, peak->sigma_pix);
  fprintf(file, "#\n");
  
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  
  output_map_t(file, kMap);
  
  fclose(file);
  if (peak->printMode < 2)  printf("\"%s\" made\n", name);
  return;
}

void outputMapFromTable(char name[], cosmo_hm *cmhm, peak_param *peak, fftw_complex *table, mapType_t type)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# %s\n", STR_MAPTYPE_T(type));
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file, "# Filter = %s, scale = %g [arcmin] = %g [pix]\n", STR_FILTER_T(peak->filter[0]), peak->scale[0], peak->linScaleInPix[0]);
  fprintf(file, "# sigma_eps = %g, sigma_pix = %g\n", peak->sigma_eps, peak->sigma_pix);
  fprintf(file, "#\n");
  
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  
  int N1 = peak->resol[0];
  int N2 = peak->resol[1];
  int M  = peak->FFTSize;
  double *pix, pos[2];
  int i, j, jM;
  
  fprintf(file, "# Number of pixels = %d\n", N1 * N2);
  fprintf(file, "#\n");
  
  if (type < 6) {
    fprintf(file, "#  theta_x   theta_y    value1\n");
    fprintf(file, "# [arcmin]  [arcmin]       [-]\n");
    
    for (j=0; j<N2; j++) {
      jM = j * M;
      for (i=0; i<N1; i++) {
	getPixPos(pos, peak->limits, peak->theta_pix, i, j);
	fprintf(file, "  %8.3f  %8.3f  %8.5f\n", pos[0], pos[1], table[i+jM][0]);
      }
    }
  }
    
  else {
    fprintf(file, "#  theta_x   theta_y    value1    value2\n");
    fprintf(file, "# [arcmin]  [arcmin]       [-]       [-]\n");
    
    for (j=0; j<N2; j++) {
      jM = j * M;
      for (i=0; i<N1; i++) {
	pix = table[i+jM];
	getPixPos(pos, peak->limits, peak->theta_pix, i, j);
	fprintf(file, "  %8.3f  %8.3f  %8.5f  %8.5f\n", pos[0], pos[1], pix[0], pix[1]);
      }
    }
  }
  
  fclose(file);
  if (peak->printMode < 2)  printf("\"%s\" made\n", name);
  return;
}

void outfitsMapFromTable(char name[], peak_param *peak, fftw_complex *table, map_t *kMap)
{
  //-- Copy values to kMap->value2
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  int M = peak->FFTSize;
  double *kappa = kMap->value2;
  int i, j, jN1, jM;
  
  for (j=0; j<N2; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<N1; i++) kappa[i+jN1] = table[i+jM][0]; //-- kappa
  }
  
  FITS_t *fits = initializeImageWriter_FITS_t(name);
  write2DImage_double(fits, TDOUBLE, FLOAT_IMG, N1, N2, (void*)kappa);
  free_FITS_t(fits);
  
  if (peak->printMode < 2)  printf("\"%s\" made\n", name);
  return;
}

void outputMask(char name[], peak_param *peak, gal_map *gMap, map_t *kMap, error **err)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# %s\n", STR_MAPTYPE_T(kMap->type));
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file, "#\n");
  
  setFillingThreshold(peak, gMap, err); forwardError(*err, __LINE__,);
  double threshold = gMap->fillingThreshold;
  gal_list **map = gMap->map;
  int i;
  for (i=0; i<gMap->length; i++) {
    if (threshold > (double)map[i]->size) kMap->value1[i] = 1.0; //-- Masked
    else                                  kMap->value1[i] = 0.0;
  }
  
  kMap->type = mask_map;
  output_map_t(file, kMap);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to map making

void pixelization(peak_param *peak, gal_map *gMap, FFT_arr *smoother, error **err)
{
  int i;
  
  if (peak->doSmoothing == 2) {
    //-- Direct convolution
    smoothAndBinByDC(peak, gMap, smoother); //-- smoother reset in this function
  }
  
  else {
    //-- Bin and stock information in smoother->array[0]->before
    signalMean_gal_map(peak, gMap, err); forwardError(*err, __LINE__,);
    galaxyBinning(gMap, smoother->array[0]); //-- Only in the first one
    if (peak->printMode < 2) printf("Binned galaxies\n");
  }
  return;
}

void makeMap(peak_param *peak, gal_map *gMap, FFT_arr *smoother, map_t *kMap, error **err)
{
  //-- Map making main function: kappa/gamma/g, noiseless/noisy, unsmoothed/smoothed
  
  char option[STRING_LENGTH_MAX];
  char input[STRING_LENGTH_MAX];
  char output[STRING_LENGTH_MAX];
  char log[STRING_LENGTH_MAX];
  
  pixelization(peak, gMap, smoother, err); forwardError(*err, __LINE__,);
  
  //-- Inversion
  if      (peak->doKappa == 2) invertByIterKS_arr(peak, gMap, smoother, kMap);
  else if (peak->doKappa == 1) copyToAllBefore(smoother);
  else                         invertByLinKS_arr(peak, gMap, smoother);
  
  if (peak->doSmoothing == 3 || peak->doSmoothing == 4) {
    //-- MRLens
    sprintf(option, "-k -C2 -c2 -s0.05 -F1 -n%d", peak->nonlinScale[0]);
    sprintf(input, "kappaMap_unsmoothed_MPI%d.fits", peak->MPIInd);
    sprintf(output, "kappaMap_mrlens_MPI%d.fits", peak->MPIInd);
    sprintf(log, "> log_mrlens_MPI%d", peak->MPIInd);
    outfitsMapFromTable(input, peak, smoother->array[0]->before, kMap);
    MRLensFiltering(peak, "../../mrlens/wl_t2_filter", option, input, output, log);
  }
  
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) {
    //-- FFT smoothing
    smoothByFFT_arr(peak, gMap, smoother);
  }
  return;
}

void makeTrueMap(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, FFT_arr *smoother, map_t *kMap, int semiTruth, error **err)
{
  int doKappa = peak->doKappa;
  if (doKappa == 1 && semiTruth == 1) return;
  if (doKappa == 0 && semiTruth == 0) return;
  
  reset_gal_map(gMap);
  read_gal_map(name, cmhm, peak, gMap, err); forwardError(*err, __LINE__,);
  
  //-- Semi-true kappa
  if (semiTruth == 1) {
    pixelization(peak, gMap, smoother, err);     forwardError(*err, __LINE__,);
    if (doKappa == 2) invertByIterKS_arr(peak, gMap, smoother, kMap);
    else              invertByLinKS_arr(peak, gMap, smoother);
  }
  
  //-- True kappa
  else {
    peak->doKappa = 1;
    pixelization(peak, gMap, smoother, err);     forwardError(*err, __LINE__,);
    peak->doKappa = doKappa;
  }
  return;
}
  
void makeMapAndOutputAll(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, FFT_arr *smoother, map_t *kMap, error **err)
{
  //-- Map making main function: kappa/gamma/g/g-linear, noiseless/noisy, unsmoothed/smoothed
  //--
  //-- Bin kappa, FFT (1, 1)                           | DC (1, 2)
  //--   Catalogue:     kappa          in gMap         |   Catalogue:          kappa          in gMap
  //--   Binning:       binned kappa   in smoo->before |   Direct convolution: smoothed kappa in smoo->before
  //--   FFT smoothing: smoothed kappa in smoo->after  |
  //--                                                 |
  //-- Bin gamma, linear KS, FFT (0, 1)                | DC, linear KS (0, 2)
  //--   Catalogue:     gamma          in gMap         |   Catalogue:          gamma          in gMap
  //--   Binning:       binned gamma   in smoo->before |   Direct convolution: smoothed gamma in smoo->before
  //--   Linear KS:     binned kappa   in smoo->before |   Linear KS:          smoothed kappa in smoo->before
  //--   FFT smoothing: smoothed kappa in smoo->after  |
  //--                                                 |
  //-- Bin g, iterative KS, FFT (2, 1)                 | DC, iterative KS (2, 2)
  //--   Catalogue:     g              in gMap         |   Catalogue:          g              in gMap
  //--   Binning:       binned g       in smoo->before |   Direct convolution: smoothed g     in smoo->before
  //--   Iterative KS:  kappa          in smoo->before |   Iterative KS:       smoothed kappa in smoo->before
  //--   FFT smoothing: smoothed kappa in smoo->after  |
  //--                                                 |
  //-- Bin g, linear KS, FFT (3, 1)                    | DC, iterative KS (3, 2)
  //--   Catalogue:     g              in gMap         |   Catalogue:          g              in gMap
  //--   Binning:       binned g       in smoo->before |   Direct convolution: smoothed g     in smoo->before
  //--   Linear KS:     kappa          in smoo->before |   Linear KS:          smoothed kappa in smoo->before
  //--   FFT smoothing: smoothed kappa in smoo->after  |
  
  char option[STRING_LENGTH_MAX];
  int doKappa       = peak->doKappa;
  int doSmoothing   = peak->doSmoothing;
  FFT_t *firstSmoo  = smoother->array[0];
  FFT_t *outputSmoo = smoother->array[0]; //-- This index is changeable.
  
  //-- Mask
  outputMask("mask", peak, gMap, kMap, err); forwardError(*err, __LINE__,);
  
  //-- Pixelization
  pixelization(peak, gMap, smoother, err);   forwardError(*err, __LINE__,);
  if (doSmoothing == 2) {
    if (doKappa == 1) outputMapFromTable("kappaMap",    cmhm, peak, outputSmoo->before, gMap->type);
    else              outputMapFromTable("gammaOrGMap", cmhm, peak, outputSmoo->before, gMap->type);
  }
  else {
    if (doKappa == 1) outputMapFromTable("kappaMap_unsmoothed",    cmhm, peak, firstSmoo->before, gMap->type);
    else              outputMapFromTable("gammaOrGMap_unsmoothed", cmhm, peak, firstSmoo->before, gMap->type);
  }
  
  //-- Inversion
  if (doKappa == 2) {
    invertByIterKS_arr(peak, gMap, smoother, kMap);
    if (doSmoothing == 2) outputMapFromTable("kappaMap",            cmhm, peak, outputSmoo->before, gMap->type);
    else                  outputMapFromTable("kappaMap_unsmoothed", cmhm, peak, firstSmoo->before,  gMap->type);
  }
  else if (doKappa == 0 || doKappa == 3) {
    invertByLinKS_arr(peak, gMap, smoother);
    if (doSmoothing == 2) outputMapFromTable("kappaMap",            cmhm, peak, outputSmoo->before, gMap->type);
    else                  outputMapFromTable("kappaMap_unsmoothed", cmhm, peak, firstSmoo->before,  gMap->type);
  }
  
  if (doSmoothing == 3 || doSmoothing == 4) {
    //-- MRLens
    sprintf(option, "-k -C2 -c2 -s0.05 -F1 -n%d", peak->nonlinScale[0]);
    outfitsMapFromTable("kappaMap_unsmoothed.fits", peak, firstSmoo->before, kMap);
    MRLensFiltering(peak, "../../mrlens/wl_t2_filter", option, "kappaMap_unsmoothed.fits", "kappaMap_mrlens.fits", "> log_mrlens");
  }
  
  if (doSmoothing == 1 || doSmoothing == 4) {
    //-- FFT smoothing
    smoothByFFT_arr(peak, gMap, smoother);
    outputMapFromTable("kappaMap", cmhm, peak, outputSmoo->after, gMap->type);
  }
  
  if (doKappa != 1) {
    //-- Semi-truth
    makeTrueMap("galCat", cmhm, peak, gMap, smoother, kMap, 1, err); forwardError(*err, __LINE__,);
    outputMapFromTable("kappaMap_semiTruth", cmhm, peak, firstSmoo->before, gMap->type);
  }
  
  if (doKappa != 0) {
    //-- Truth
    makeTrueMap("galCat", cmhm, peak, gMap, smoother, kMap, 0, err); forwardError(*err, __LINE__,);
    outputMapFromTable("kappaMap_truth", cmhm, peak, firstSmoo->before, gMap->type);
  }
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doKMap(char fileName[], cosmo_hm *cmhm, peak_param *peak, int doNoise, error **err)
{
  int N1_mask = (int)(peak->Omega[0] * peak->theta_CCD_inv);
  int N2_mask = (int)(peak->Omega[1] * peak->theta_CCD_inv);
  peak->doNoise     = doNoise;
  
  halo_map *hMap     = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                     forwardError(*err, __LINE__,);
  gal_map *gMap      = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask = initialize_short_mat(N1_mask, N2_mask);
  if (peak->doMask == 1) fillMask_CFHTLenS_W1(peak, CCDMask);
  FFT_arr *smoother  = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernel(peak, smoother);
  map_t *kMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  
  if (fileName == NULL) {
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
    setMassSamplers(cmhm, peak, sampArr, err);                    forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err);                forwardError(*err, __LINE__,);
    outputFastSimul("haloCat", cmhm, peak, hMap);
    free_sampler_arr(sampArr);
  }
  else {
    read_halo_map(fileName, cmhm, hMap, err);                     forwardError(*err, __LINE__,);
  }
  
  cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err); forwardError(*err, __LINE__,);
  lensingCatalogueAndOutputAll(cmhm, peak, hMap, gMap, err);      forwardError(*err, __LINE__,);
  makeMapAndOutputAll(cmhm, peak, gMap, smoother, kMap, err);     forwardError(*err, __LINE__,);

  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(smoother);
  free_map_t(kMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

