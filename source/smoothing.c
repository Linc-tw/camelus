

  /*******************************************************
   **  smoothing.c					**
   **  Version 2018.03.15				**
   **							**
   **  References:					**
   **  - Hetterscheidt et al. (2005) - A&A, 442, 43	**
   **  - Leonard et al. (2012) - MNRAS, 423, 34		**
   **  - Seitz & Schneider (1995) - A&A, 297, 287	**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "smoothing.h"


//----------------------------------------------------------------------
//-- Functions related to signal_map

signal_map *initialize_signal_map(int N1, int N2, double theta_pix, error **err)
{
  signal_map *kMap    = (signal_map*)malloc_err(sizeof(signal_map), err);    forwardError(*err, __LINE__, NULL);
  kMap->N1            = N1;
  kMap->N2            = N2;
  kMap->length        = N1 * N2;
  kMap->theta_pix     = theta_pix;
  kMap->theta_pix_inv = 1.0 / theta_pix;
  
  kMap->type   = kappa_map;
  kMap->value1 = (double*)malloc_err(kMap->length * sizeof(double), err); forwardError(*err, __LINE__, NULL);
  kMap->value2 = (double*)malloc_err(kMap->length * sizeof(double), err); forwardError(*err, __LINE__, NULL);
  return kMap;
}

void free_signal_map(signal_map *kMap)
{
  if (kMap) {
    if (kMap->value1) {free(kMap->value1); kMap->value1 = NULL;}
    if (kMap->value2) {free(kMap->value2); kMap->value2 = NULL;}
    free(kMap); kMap = NULL;
  }
  return;
}

void read_signal_map(char *name, peak_param *pkPar, signal_map *kMap, map_t type, error **err)
{
  //-- Open
  testErrorRet(kMap->length!=pkPar->resol[0]*pkPar->resol[1], peak_match, "Resolution match error", *err, __LINE__,);
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  
  double *value1 = kMap->value1;
  double *value2 = kMap->value2;
  int count = 0;
  char buffer[STRING_LENGTH_MAX];
  
  //-- Read
  int c = fgetc(file);
  if (type < 6) {
    while (c != EOF) {
      if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
      else {
	testErrorRet(count>=kMap->length, peak_overflow, "Too many pixels", *err, __LINE__,);
	ungetc(c, file);
	fscanf(file, "%*f %*f %lf\n", &value1[count]);
	count++;
      }
      c = fgetc(file);
    }
  }
  else {
    while (c != EOF) {
      if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
      else {
	testErrorRet(count>=kMap->length, peak_overflow, "Too many pixels", *err, __LINE__,);
	ungetc(c, file);
	fscanf(file, "%*f %*f %lf %lf\n", &value1[count], &value2[count]);
	count++;
      }
      c = fgetc(file);
    }
  }
  fclose(file);
  
  //-- Check
  testErrorRet(count!=kMap->length, peak_match, "Array length error", *err, __LINE__,);
  kMap->type = type;
  printf("Read \"%s\"\n", name);
  return;
}

void getPixPos(double pos[2], double theta_pix, int i, int j)
{
  //-- Compute the position of the pixel indexed (i, j) and stock this information in pos.
  pos[0] = (0.5 + i) * theta_pix;
  pos[1] = (0.5 + j) * theta_pix;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to smoothing by FFT

void fillGaussianKernel(fftw_complex *kernel, int length, int M, double scaleInPix, int doVar)
{
  //-- M is resolution + size of zero padding
  //-- 
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
  //-- For variance:
  //-- If K = sum(w_i kappa_i) / sum(w_i), then variance = sum(w_i^2 sigma_i^2) / sum(w_i)^2.
  
  //-- Reset
  reset_fftw_complex(kernel, length);
  
  double scaleInPix_invSq = 1.0 / (SQ(scaleInPix)); //-- [pix^-2]
  double cutInPix    = CUTOFF_FACTOR_GAUSSIAN * scaleInPix; //-- CUTOFF_FACTOR_GAUSSIAN set to 2.2 in peakParameters.h, equivalent to a Gaussian trunked at 3-sigma
  double cutInPix_sq = SQ(cutInPix); //-- [pix^2]
  int cutSize = (int)ceil(cutInPix); //-- [pix]
  double sum  = 0.0;
  
  double value;
  int i, j, mod_j_M, r_sq;
  
  for (j=-cutSize; j<=cutSize; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    for (i=-cutSize; i<=cutSize; i++) {
      r_sq = i*i + j*j;
      if (r_sq > cutInPix_sq) continue;
      value = exp(-r_sq * scaleInPix_invSq);
      sum  += value;
      kernel[POS_MOD(M, i) + mod_j_M][0] = (doVar == 0) ? value : SQ(value);
    }
  }
  
  //-- Normalization
  if (doVar == 0) rescaleReal_fftw_complex(kernel, length, 1.0 / sum);
  else            rescaleReal_fftw_complex(kernel, length, 1.0 / SQ(sum));
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
  
  double x_zero = fabs(pow(x,       3.0));
  double x_half = fabs(pow(x - 0.5, 3.0)) + fabs(pow(x + 0.5, 3.0));
  double x_one  = fabs(pow(x - 1.0, 3.0)) + fabs(pow(x + 1.0, 3.0));
  double x_two  = fabs(pow(x - 2.0, 3.0)) + fabs(pow(x + 2.0, 3.0));
  double y_zero = fabs(pow(y,       3.0));
  double y_half = fabs(pow(y - 0.5, 3.0)) + fabs(pow(y + 0.5, 3.0));
  double y_one  = fabs(pow(y - 1.0, 3.0)) + fabs(pow(y + 1.0, 3.0));
  double y_two  = fabs(pow(y - 2.0, 3.0)) + fabs(pow(y + 2.0, 3.0));
  
  double value = 0.0;
  //-- This is not yet normalized.
  value += 256.0 * (x_one - 4*x_half + 6*x_zero) * (y_one - 4*y_half + 6*y_zero);
  value -=         (x_two - 4*x_one  + 6*x_zero) * (y_two - 4*y_one  + 6*y_zero);
  return value;
}

void fillStarletKernel(fftw_complex *kernel, int length, int M, double scaleInPix, int doVar)
{
  //-- L2-norm of the 2D-starlet:
  //--   || psi ||_L2^2 = NORM_STARLET^2 = 5 I2^2 - 2 I3^2
  //--   where
  //--     I2 = 2/5 + 5/63
  //--     I3 = 1/3 + 1/5 + 1/21 + 1/48
  //--
  //-- Scaled 2D-starlet:
  //--   Let x = i/s, y = j/s, define psi_s(i, j) = psi(i/s, j/s) / (NORM1_STARLET s^2).
  //--   Then, || psi_s ||_L1   = 1,
  //--   and   || psi_s ||_L2^2 = NORM2_STARLET^2 / (NORM1_STARLET^2 s^2)
  //--
  //--   In this case, S/N, kappa divided by sigma_noise, is given by
  //--   sigma_noise^2 = sigma_eps^2 / (2 n_g s^2) * (NORM2_STARLET / NORM1_STARLET)^2.
  //--
  //-- Implementation:
  //--   1. Compute 144 * psi(i/s, j/s).
  //--   2. Normalize the kernel by the sum of absolute coefficients.
  //--   This cancels the factor 144 and set L1-norm = 1.
  
  //-- Reset
  reset_fftw_complex(kernel, length);
  
  int cutSize = MIN((int)ceil(2 * scaleInPix), M/2); //-- [pix]
  double sum  = 0.0;
  double x, y, value;
  int i, j, mod_j_M;
  
  for (j=-cutSize; j<=cutSize; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    y = (double)j / scaleInPix;
    for (i=-cutSize; i<=cutSize; i++) {
      x = (double)i / scaleInPix;
      value = starlet_2D(x, y);
      sum  += fabs(value); //-- Compute the sum of the absolute value of coefficients
      kernel[POS_MOD(M, i) + mod_j_M][0] = (doVar == 0) ? value : SQ(value);
    }
  }
  
  //-- Normalization
  if (doVar == 0) rescaleReal_fftw_complex(kernel, length, 1.0 / sum);
  else            rescaleReal_fftw_complex(kernel, length, 1.0 / SQ(sum));
  return;
}

#define TANH_A        6.0
#define TANH_B      150.0
#define TANH_C       47.0
#define TANH_D       50.0
#define TANH_X_C      0.1
#define TANH_X_C_INV 10.0
#define EPS          1e-9
void fillMApTanhKernel(fftw_complex *kernel, int length, int M, double scaleInPix, int doVar)
{
  //-- See Hetterscheidt et al. (2005), Eq. (9)
  //-- Dimensionless Q function:
  //--   Q(x) = tanh(x / x_c) / [pi (x / x_c) (1 + exp(a - bx) + exp(-c + dx))]
  //--   with a = 6, b = 150, c = 47, d = 50, x_c = 0.1
  //--
  //-- Implementation:
  //--   Q(x) = tanh(x / x_c) / [x (1 + exp(a - bx) + exp(-c + dx))]
  
  //-- Reset
  reset_fftw_complex(kernel, length);
  
  double scaleInPix_inv = 1.0 / scaleInPix; //-- [pix^-1]
  double cutInPix    = CUTOFF_FACTOR_M_AP_TANH * scaleInPix; //-- CUTOFF_FACTOR_M_AP_TANH set to 1.2 in peakParameters.h
  double cutInPix_sq = SQ(cutInPix); //-- [pix^2]
  int cutSize = (int)ceil(cutInPix); //-- [pix]
  double sum  = 0.0;
  
  double x, phase, factor;
  int i, j, mod_j_M, r_sq, index;
  
  for (j=-cutSize; j<=cutSize; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    for (i=-cutSize; i<=cutSize; i++) {
      r_sq = i*i + j*j;
      if (r_sq > cutInPix_sq) continue;
      
      index  = POS_MOD(M, i) + mod_j_M;
      phase  = PI + 2.0 * atan2(j, i); //-- phase = the angle in the gamma_1-gamma_2 space
      x      = sqrt(r_sq) * scaleInPix_inv;
      x      = MAX(x, EPS);
      factor = tanh(x * TANH_X_C_INV) / (x * (1.0 + exp(TANH_A - TANH_B * x) + exp(-TANH_C + TANH_D * x)));
      sum   += fabs(factor);
      if (doVar == 0) {
	kernel[index][0] =  factor * cos(phase);
	kernel[index][1] = -factor * sin(phase); //-- The sign minus is due to i.
      }
      else kernel[index][0] = SQ(factor);
    }
  }
  
  //-- Normalization
  if (doVar == 0) rescale_fftw_complex(kernel, length, 1.0 / sum);
  else            rescaleReal_fftw_complex(kernel, length, 1.0 / SQ(sum));
  return;
}
#undef TANH_A
#undef TANH_B
#undef TANH_C
#undef TANH_D
#undef TANH_X_C
#undef TANH_X_C_INV
#undef EPS

void fillMApGammaTKernel(fftw_complex *kernel, int length, int M, double scaleInPix1, double scaleInPix2, int doVar)
{
  //-- Reset
  reset_fftw_complex(kernel, length);
  
  double cutInPix1_sq = SQ(scaleInPix1); //-- [pix^2]
  double cutInPix2_sq = SQ(scaleInPix2); //-- [pix^2]
  int cutSize = (int)ceil(scaleInPix2);  //-- [pix]
  double sum  = 0.0;
  
  double phase;
  int i, j, mod_j_M, r_sq, index;
  
  for (j=-cutSize; j<=cutSize; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    for (i=-cutSize; i<=cutSize; i++) {
      r_sq = i*i + j*j;
      if (r_sq <= cutInPix1_sq) continue;
      if (r_sq > cutInPix2_sq) continue;
      
      index = POS_MOD(M, i) + mod_j_M;
      phase = PI + 2.0 * atan2(j, i); //-- phase = the angle in the gamma_1-gamma_2 space
      sum  += 1.0;
      if (doVar == 0) {
	kernel[index][0] =  cos(phase);
	kernel[index][1] = -sin(phase); //-- The sign minus is due to i.
      }
      else kernel[index][0] = 1.0;
    }
  }
  
  //-- Normalization
  if (doVar == 0) rescale_fftw_complex(kernel, length, 1.0 / sum);
  else            rescaleReal_fftw_complex(kernel, length, 1.0 / SQ(sum));
  return;
}

void makeKernel(peak_param *pkPar, FFT_arr *FFTSmoother)
{
  FFT_t *smoo;
  int i;
  for (i=0; i<pkPar->FFT_nbFilters; i++) {
    smoo = FFTSmoother->array[i];
    if      (pkPar->FFT_filter[i] == 0) fillGaussianKernel(smoo->kernel, smoo->length, smoo->N, pkPar->FFT_scaleInPix[i], 0); //-- doVar = 0
    else if (pkPar->FFT_filter[i] == 1) fillStarletKernel(smoo->kernel, smoo->length, smoo->N, pkPar->FFT_scaleInPix[i], 0);  //-- doVar = 0
    else if (pkPar->FFT_filter[i] == 2) fillMApTanhKernel(smoo->kernel, smoo->length, smoo->N, pkPar->FFT_scaleInPix[i], 0);  //-- doVar = 0
    else if (pkPar->FFT_filter[i] == 3) fillMApGammaTKernel(smoo->kernel, smoo->length, smoo->N, pkPar->FFT_scaleInPix[i], pkPar->FFT_scaleInPix[i+1], 0); //-- doVar = 0
    fftw_execute(smoo->kernel_f); //-- To Fourier space
  }
  return;
}

void smoothByFFT_arr(peak_param *pkPar, FFT_arr *FFTSmoother)
{
  execute_FFT_arr(FFTSmoother);
  if (pkPar->verbose < 3) printf("Smoothed by FFT\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to smoothing and binning by direct convolution

void DCForPair_kappa(peak_param *pkPar, gal_list *gList, FFT_t **smooArr, double pixPos[2], int index_FFT)
{
  //-- For pixel-pixel pair
  
  int DC_nbFilters     = pkPar->DC_nbFilters;
  double *DC_scale_inv = pkPar->DC_scale_inv;
  double *DC_cut       = pkPar->DC_cut;
  double sin_ctrDEC    = sin(pixPos[1]);
  double cos_ctrDEC    = cos(pixPos[1]);
  
  gal_node *gNode;
  gal_t *g;
  double *pix; 
  double dtheta, dtheta_sq, dtheta_x, dtheta_y, dRA, angle, factor;
  int i, j, dthetaComputed, dthetaSqComputed;
  
  for (i=0, gNode=gList->first; i<gList->size; i++, gNode=gNode->next) {
    g = gNode->g;
    dthetaComputed   = 0;
    dthetaSqComputed = 0;
    
    for (j=0; j<DC_nbFilters; j++) {
      //-- Gaussian
      if (pkPar->DC_filter[j] == 0) {
	//-- Compute dtheta_sq, field = 1 already excluded
	if (dthetaSqComputed == 0) {
	  if (pkPar->field == 0) dtheta_sq = DIST_2D_SQ(pixPos, g->pos);
	  else                   dtheta_sq = SQ(SPHE_DIST(pixPos, g->pos));
	  dthetaSqComputed = 1;
	}
	
	if (dtheta_sq > DC_cut[j]) continue; //-- DC_cut & DC_scale_inv squared for the Gaussian
	factor  = exp(-dtheta_sq * DC_scale_inv[j]) * g->weight; //-- factor is a weight from the kernel time a weight from the shape measurement.
	pix     = smooArr[j]->before[index_FFT];
	pix[0] += factor * g->kappa;
	pix     = smooArr[j]->kernel[index_FFT];
	pix[0] += SQ(factor);   //-- Compute the smoothed variance
	pix[1] += fabs(factor); //-- Compute the smoothed number count (sum of weights)
      }
      
      //-- Starlet
      else if (pkPar->DC_filter[j] == 1) {
	//-- Compute dtheta, field = 1 already excluded
	if (dthetaComputed == 0) {
	  if (pkPar->field == 0) {
	    dtheta_x = g->pos[0] - pixPos[0];
	    dtheta_y = g->pos[1] - pixPos[1];
	  }
	  else { //-- field = 1 already excluded
	    dtheta   = SPHE_DIST(pixPos, g->pos);
	    dRA      = g->pos[0] - pixPos[0];
	    dtheta_x = g->sinCosDEC[1] * sin(dRA);
	    dtheta_y = (cos_ctrDEC * g->sinCosDEC[0] - sin_ctrDEC * g->sinCosDEC[1] * cos(dRA));
	    angle    = atan2(dtheta_y, dtheta_x) + pkPar->rotAng; //-- For field = 0 or 2, rotAng forced to 0.0; angle = the angle of the relative position
	    dtheta_x = dtheta * cos(angle);
	    dtheta_y = dtheta * sin(angle);
	  }
	  dthetaComputed = 1;
	}
	
	if (fabs(dtheta_x) > DC_cut[j]) continue;
	if (fabs(dtheta_y) > DC_cut[j]) continue;
	factor  = starlet_2D(dtheta_x * DC_scale_inv[j], dtheta_y * DC_scale_inv[j]) * g->weight;
	pix     = smooArr[j]->before[index_FFT];
	pix[0] += factor * g->kappa;
	pix     = smooArr[j]->kernel[index_FFT];
	pix[0] += SQ(factor);   //-- Compute the smoothed variance
	pix[1] += fabs(factor); //-- Compute the smoothed number count (sum of weights)
      }
      
      //-- Do nothing for aperture mass
      else if (pkPar->DC_filter[j] == 2) ;
      else if (pkPar->DC_filter[j] == 3) ;
    }
  }
  return;
}

#define TANH_A        6.0
#define TANH_B      150.0
#define TANH_C       47.0
#define TANH_D       50.0
#define TANH_X_C      0.1
#define TANH_X_C_INV 10.0
#define EPS          1e-9
void DCForPair_gamma(peak_param *pkPar, gal_list *gList, FFT_t **smooArr, double pixPos[2], int index_FFT)
{
  //-- For pixel-pixel pair
  
  int DC_nbFilters     = pkPar->DC_nbFilters;
  double *DC_scale_inv = pkPar->DC_scale_inv;
  double *DC_cut       = pkPar->DC_cut;
  double sin_ctrDEC    = sin(pixPos[1]);
  double cos_ctrDEC    = cos(pixPos[1]);
  
  gal_node *gNode;
  gal_t *g;
  double *pix;
  double dtheta, dtheta_x, dtheta_y, dRA, phase, gamma_t;
  double x, factor;
  int i, j, phaseComputed;
  
  for (i=0, gNode=gList->first; i<gList->size; i++, gNode=gNode->next) {
    g = gNode->g;
    
    //-- Compute dtheta, field = 1 already excluded
    if (pkPar->field == 0) dtheta = sqrt(DIST_2D_SQ(pixPos, g->pos));
    else                   dtheta = SPHE_DIST(pixPos, g->pos);
    phaseComputed = 0;
    
    for (j=0; j<DC_nbFilters; j++) {
      //-- Do nothing for Gaussian & starlet
      if (pkPar->DC_filter[j] == 0) ;
      else if (pkPar->DC_filter[j] == 1) ;
      
      //-- M_ap tanh
      else if (pkPar->DC_filter[j] == 2) {
	if (dtheta > DC_cut[j]) continue;
	
	//-- Compute phase & gamma_t, field = 1 already excluded
	if (phaseComputed == 0) {
	  if (pkPar->field == 0) {
	    dtheta_x = g->pos[0] - pixPos[0];
	    dtheta_y = g->pos[1] - pixPos[1];
	  }
	  else {
	    dRA      = g->pos[0] - pixPos[0];
	    dtheta_x = g->sinCosDEC[1] * sin(dRA);
	    dtheta_y = (cos_ctrDEC * g->sinCosDEC[0] - sin_ctrDEC * g->sinCosDEC[1] * cos(dRA));
	  }
	  phase   = PI + 2.0 * (atan2(dtheta_y, dtheta_x) + pkPar->rotAng); //-- phase = the angle in the gamma_1-gamma_2 space, see rayTracing.c
	  gamma_t = g->gamma[0] * cos(phase) + g->gamma[1] * sin(phase);
	  phaseComputed = 1;
	}
	
	x       = dtheta * DC_scale_inv[j];
	x       = MAX(x, EPS);
	factor  = tanh(x * TANH_X_C_INV) / (x * (1.0 + exp(TANH_A - TANH_B * x) + exp(-TANH_C + TANH_D * x))) * g->weight; //-- factor is a weight from the kernel time a weight from the shape measurement.
	pix     = smooArr[j]->before[index_FFT];
	pix[0] += factor * gamma_t;
	pix     = smooArr[j]->kernel[index_FFT];
	pix[0] += SQ(factor);   //-- Compute the smoothed variance
	pix[1] += fabs(factor); //-- Compute the smoothed number count (sum of weights)
      }
      
      //-- M_ap gamma_t
      else if (pkPar->DC_filter[j] == 3) {
	if (dtheta <= DC_cut[j]) continue;
	if (dtheta > DC_cut[j+1]) continue;
	
	//-- Compute phase & gamma_t, field = 1 already excluded
	if (phaseComputed == 0) {
	  if (pkPar->field == 0) {
	    dtheta_x = g->pos[0] - pixPos[0];
	    dtheta_y = g->pos[1] - pixPos[1];
	  }
	  else {
	    dRA      = g->pos[0] - pixPos[0];
	    dtheta_x = g->sinCosDEC[1] * sin(dRA);
	    dtheta_y = (cos_ctrDEC * g->sinCosDEC[0] - sin_ctrDEC * g->sinCosDEC[1] * cos(dRA));
	  }
	  phase   = PI + 2.0 * (atan2(dtheta_y, dtheta_x) + pkPar->rotAng); //-- phase = the angle in the gamma_1-gamma_2 space, see rayTracing.c
	  gamma_t = g->gamma[0] * cos(phase) + g->gamma[1] * sin(phase);
	  phaseComputed = 1;
	}
	
	pix     = smooArr[j]->before[index_FFT];
	pix[0] += g->weight * gamma_t;
	pix     = smooArr[j]->kernel[index_FFT];
	pix[0] += SQ(g->weight);   //-- Compute the smoothed variance
	pix[1] += fabs(g->weight); //-- Compute the smoothed number count (sum of weights)
      }
    }
  }
  return;
}
#undef TANH_A
#undef TANH_B
#undef TANH_C
#undef TANH_D
#undef TANH_X_C
#undef TANH_X_C_INV
#undef EPS

void DCForPixel(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, DC_fct *DCForPairFct, double pixPos[2], int index_FFT, int i_pix, int j_pix)
{
  //-- For a map pixel and the whole galaxy catalogue
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int bufferSize  = pkPar->bufferSize;
  gal_list **map  = gMap->map;
  FFT_t **smooArr = DCSmoother->array;
  
  int i, j, jN1;
  
  for (j=j_pix-bufferSize; j<=j_pix+bufferSize; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_pix-bufferSize; i<=i_pix+bufferSize; i++) {
      if (i < 0 || i >= N1) continue;
      DCForPairFct(pkPar, map[i+jN1], smooArr, pixPos, index_FFT);
    }
  }
  
  double *pix;
  double totWeight;
  
  //-- Normalization
  for (i=0; i<DCSmoother->length; i++) {
    pix       = smooArr[i]->kernel[index_FFT];
    totWeight = pix[1];
    if (totWeight == 0.0) continue;
    totWeight = 1.0 / totWeight;
    pix[0]   *= SQ(totWeight);
    pix       = smooArr[i]->before[index_FFT];
    pix[0]   *= totWeight;
    pix[1]   *= totWeight;
  }
  return;
}

void DCForMap(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, DC_fct *DCForPairFct)
{
  int M  = DCSmoother->array[0]->N;
  int N1 = (pkPar->field < 2) ? gMap->N1 : M;
  int N2 = (pkPar->field < 2) ? gMap->N2 : M;
  int verbose      = pkPar->verbose;
  int DC_nbFilters = DCSmoother->length;
  int FFTLength    = DCSmoother->array[0]->length;
  long nsidePix    = pkPar->nside * pkPar->HP_resol;
  long first       = pkPar->nest * pkPar->HP_resol * pkPar->HP_resol;
  double length        = (pkPar->field < 2) ? (double)gMap->length : (double)FFTLength;
  double theta_pix     = gMap->theta_pix;
  double theta_pix_inv = gMap->theta_pix_inv;
  double offset        = pkPar->offset;
  
  double pixPos[2], projPixPos[2];
  int localNest, i_pix, j_pix;
  int i, j, k, jN1, index_gMap;
  
  //-- Reset all tables
  for (k=0; k<DC_nbFilters; k++) reset_fftw_complex(DCSmoother->array[k]->before, FFTLength);
  for (k=0; k<DC_nbFilters; k++) reset_fftw_complex(DCSmoother->array[k]->kernel, FFTLength);
  
  //-- Direct convolution
  for (j=0; j<N2; j++) {
    jN1 = j * N1;
    for (i=0; i<N1; i++) {
      index_gMap = i + jN1;
      
      if (verbose < 2 && index_gMap % 50 == 0) {
	printf("Computing direct convolution: %6.2f%%\r", 100.0 * index_gMap / length);
	fflush(stdout);
      }
      
      if (pkPar->field == 0) {
	pixPos[0] = (0.5 + i) * theta_pix;
	pixPos[1] = (0.5 + j) * theta_pix;
	i_pix = i;
	j_pix = j;
      }
      else {
#ifdef __CAMELUS_USE_HEALPIX__
	ijPixToLocalNest(pkPar->HP_resol, i, j, &localNest);
	pix2ang_nest(nsidePix, first+localNest, &pixPos[1], &pixPos[0]); //-- theta & phi
	pixPos[1] = HALF_PI - pixPos[1]; //-- DEC in [rad]
	
	project(pixPos, pkPar->center, projPixPos);
	rotate(projPixPos, pkPar->rotAng, projPixPos);
	i_pix = (int)((projPixPos[0] + offset) * theta_pix_inv);
	j_pix = (int)((projPixPos[1] + offset) * theta_pix_inv);
#endif
      }
      
      DCForPixel(pkPar, gMap, DCSmoother, DCForPairFct, pixPos, i+j*M, i_pix, j_pix);
    }
  }
  
  //-- Rescale variance
  if (pkPar->doNoise) {
    for (k=0; k<DC_nbFilters; k++) rescaleReal_fftw_complex(DCSmoother->array[k]->kernel, FFTLength, SQ(pkPar->sigma_half));
  }
  if (pkPar->verbose < 3) printf("Computed direct convolution          \n");
  return;
}

void smoothAndBinByDC(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother)
{
  //-- Stock signal in DCSmoother->array[i]->before
  //-- Stock variance in DCSmoother->array[i]->kernel[ind][0]
  //-- Stock filling factor in DCSmoother->array[i]->kernel[ind][1]
  
  DC_fct *DCForPairFct;
  if (pkPar->doKappa == 1) DCForPairFct = DCForPair_kappa;
  else                     DCForPairFct = DCForPair_gamma;
  DCForMap(pkPar, gMap, DCSmoother, DCForPairFct);
  gMap->type += 1;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to linear KS inversion

void gammaToKappa(fftw_complex *table, int M)
{
  //-- \hat{kappa} = [Re(\hat{gamma}) + i Im(\hat{gamma})] * [(l1^2 - l2^2 - i * 2 l1 l2) / (l1^2 + l2^2)]
  
  double l1, l2, l1_sq, l2_sq, l_norm_inv, f1, f2;
  double *pix, value;
  int jM;
  
  int i, j;
  for (j=0; j<M; j++) {
    l2 = CEN_MOD(M, j);
    jM = j * M;
    
    for (i=0; i<M; i++) {
      if (i + j == 0) continue; //-- This means if i = 0 and j = 0
      pix        = table[i+jM];
      l1         = CEN_MOD(M, i);
      l1_sq      = SQ(l1);
      l2_sq      = SQ(l2);
      l_norm_inv = 1.0 / (l1_sq + l2_sq);
      f1         = (l1_sq - l2_sq) * l_norm_inv;
      f2         = 2 * l1 * l2 * l_norm_inv;
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

void invertByLinKS(FFT_t *smoo)
{
  fftw_complex *smoo_before = smoo->before;
    
  fftw_execute(smoo->before_f);                                          //-- Go to Fourier space
  gammaToKappa(smoo_before, smoo->N);                                    //-- KS inversion
  fftw_execute(smoo->before_b);                                          //-- Go to direct space
  rescaleReal_fftw_complex(smoo_before, smoo->length, smoo->normFactor); //-- Rescale, only real part is interesting.
  return;
}

void linKSAndFFT(FFT_t *smoo)
{
  int FFTLength = smoo->length;
  fftw_complex *smoo_before = smoo->before;
  fftw_complex *smoo_after  = smoo->after;
  
  //-- Linear KS
  fftw_execute(smoo->before_f);       //-- Go to Fourier space
  gammaToKappa(smoo_before, smoo->N); //-- KS inversion
  
  //-- FFT smoothing
  multiplication_fftw_complex(smoo_before, smoo->kernel, smoo_after, FFTLength); //-- Multiplication
  fftw_execute(smoo->after_b);                                                   //-- Go to direct space
  rescaleReal_fftw_complex(smoo_after, FFTLength, smoo->normFactor);             //-- Rescale, only real part is interesting.
  return;
}

void linKSAndFFT_arr(peak_param *pkPar, FFT_arr *FFTSmoother)
{
  //-- Not used
  
  FFT_t *smoo = FFTSmoother->array[0];
  linKSAndFFT(smoo);
  
  fftw_complex *smoo_before = smoo->before;
  int FFTLength = smoo->length;
  int i;
  
  for (i=1; i<FFTSmoother->length; i++) {
    smoo = FFTSmoother->array[i];
    multiplication_fftw_complex(smoo_before, smoo->kernel, smoo->after, FFTLength); //-- Multiplication
    fftw_execute(smoo->after_b);                                                    //-- Go to direct space
    rescaleReal_fftw_complex(smoo->after, FFTLength, smoo->normFactor);             //-- Rescale, only real part is interesting.
  }
  if (pkPar->verbose < 3) printf("Inverted by linear KS and smoothed by FFT\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to iterative KS inversion

#define EPS 1e-5
void iterativeKS(fftw_complex *reducedShear, fftw_complex *table, fftw_plan forward, fftw_plan backward, signal_map *kMap, int M)
{
  //-- At output, kappa is stocked in kMap->value1 and table.
  
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  int FFTLength = M * M;
  double FFTNormFactor= 1.0 / (double)FFTLength;
  double *copy = kMap->value1;
  
  reset_double(copy, kMap->length);
  
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
    gammaToKappa(table, M);                                    //-- KS inversion
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

void invertByIterKS(FFT_t *smoo, signal_map *kMap)
{
  fftw_complex *smoo_before = smoo->before;
  fftw_complex *smoo_after  = smoo->after;
  copy_fftw_complex(smoo_before, smoo_after, smoo->length);                            //-- Move reduced shear to smoo->after
  iterativeKS(smoo_after, smoo_before, smoo->before_f, smoo->before_b, kMap, smoo->N); //-- Iterative KS
  return;
}

//----------------------------------------------------------------------
//-- Functions related to Seitz-Schneider inversion


void kappaToGamma(fftw_complex *table, int M)
{
  //-- \hat{gamma} = [Re(\hat{kappa}) + i Im(\hat{kappa})] * [(l1^2 - l2^2 + i * 2 l1 l2) / (l1^2 + l2^2)]
  //-- Almost identical to gammaToKappa
  
  double l1, l2, l1_sq, l2_sq, l_norm_inv, f1, f2;
  double *pix, value;
  int jM;
  
  int i, j;
  for (j=0; j<M; j++) {
    l2 = CEN_MOD(M, j);
    jM = j * M;
    
    for (i=0; i<M; i++) {
      if (i + j == 0) continue; //-- This means if i = 0 and j = 0
      pix        = table[i+jM];
      l1         = CEN_MOD(M, i);
      l1_sq      = SQ(l1);
      l2_sq      = SQ(l2);
      l_norm_inv = 1.0 / (l1_sq + l2_sq);
      f1         = (l1_sq - l2_sq) * l_norm_inv;
      f2         = 2 * l1 * l2 * l_norm_inv;
      value      = f1 * pix[0] - f2 * pix[1];
      pix[1]     = f2 * pix[0] + f1 * pix[1];
      pix[0]     = value;
    }
  }
  
  table[0][0] = 0.0;
  table[0][1] = 0.0;
  return;
}

#define EPS 1e-5
void SeitzSchneider(fftw_complex *reducedShear, fftw_complex *table, fftw_plan forward, fftw_plan backward, signal_map *kMap, int M)
{
  //-- See Seitz & Schneider (1995)
  //-- At output, kappa is stocked in kMap->value1 and table.
  
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  int FFTLength = M * M;
  double FFTNormFactor= 1.0 / (double)FFTLength;
  double *copy    = kMap->value1;
  double *norm_sq = kMap->value2;
  
  reset_double(copy, kMap->length);
  reset_double(norm_sq, kMap->length);
  
  double factor;
  int i, j, jN1, jM, index_kMap, index_FFT;
  
  //-- Transform g into delta
  for (j=0; j<M; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<M; i++) {
      index_kMap = i + jN1;
      index_FFT  = i + jM;
      if (i >= N1 || j >= N2) {
	reducedShear[index_FFT][0] = 0.0;
	reducedShear[index_FFT][1] = 0.0;
	continue;
      }
      
      norm_sq[index_kMap] = SQ(reducedShear[index_FFT][0]) + SQ(reducedShear[index_FFT][1]); //-- |g|^2
      factor = 2.0 / (1.0 + norm_sq[index_kMap]); //-- 2 / (1 + |g|^2)
      reducedShear[index_FFT][0] *= factor;       //-- delta_1
      reducedShear[index_FFT][1] *= factor;       //-- delta_2
      norm_sq[index_kMap] *= SQ(factor);          //-- |delta|^2
    }
  }
  
  double diff, sign;
  
  //-- Iterations
  do {
    //-- Compute old gamma:
    //-- \hat{gamma} = [Re(\hat{kappa}) + i Im(\hat{kappa})] * [(l1^2 - l2^2 + i * 2 l1 l2) / (l1^2 + l2^2)]
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
	table[index_FFT][0] = copy[i+jN1];
	table[index_FFT][1] = 0.0;
      }
    }
    
    fftw_execute(forward);                                 //-- Go to Fourier space
    kappaToGamma(table, M);                                //-- Comput \hat{gamma}
    fftw_execute(backward);                                //-- Go to direct space
    rescale_fftw_complex(table, FFTLength, FFTNormFactor); //-- Rescale complex
    
    //-- Compute sign
    for (j=0; j<N2; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<N1; i++) {
	index_FFT = i + jM;
	sign = SQ(1.0 - copy[i+jN1]) - SQ(table[index_FFT][0]) - SQ(table[index_FFT][1]); //-- (1 - kappa)^2 -|gamma|^2
	table[index_FFT][0] = (sign > 0.0) ? 1.0 : -1.0; //-- Stock sign
      }
    }
    
    //-- Compute new gamma:
    //-- gamma = (1 - sgn * sqrt(1 - |delta|^2)) * (delta / |delta|^2) * (1 - kappa)
    for (j=0; j<M; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<M; i++) {
	index_kMap = i + jN1;
	index_FFT  = i + jM;
	if (i >= N1 || j >= N2) {
	  table[index_FFT][0] = 0.0;
	  table[index_FFT][1] = 0.0;
	  continue;
	}
	
	factor  = 1.0 - table[index_FFT][0] * sqrt(1.0 - norm_sq[index_kMap]); //-- 1 - sgn * sqrt(1 - |delta|^2)
	factor /= MAX(EPS, norm_sq[index_kMap]);                               //-- Devide by |delta|^2
	factor *= 1.0 - copy[index_kMap];                                      //-- Multiply by 1 - kappa
	table[index_FFT][0] = reducedShear[index_FFT][0] * factor;             //-- gamma_1
	table[index_FFT][1] = reducedShear[index_FFT][1] * factor;             //-- gamma_2
      }
    }
    
    fftw_execute(forward);                                     //-- Go to Fourier space
    gammaToKappa(table, M);                                    //-- KS inversion
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

void invertBySS(FFT_t *smoo, signal_map *kMap)
{
  fftw_complex *smoo_before = smoo->before;
  fftw_complex *smoo_after  = smoo->after;
  copy_fftw_complex(smoo_before, smoo_after, smoo->length);                               //-- Move reduced shear to smoo->after
  SeitzSchneider(smoo_after, smoo_before, smoo->before_f, smoo->before_b, kMap, smoo->N); //-- Iterative KS
  return;
}

//----------------------------------------------------------------------
//-- Functions related to ASCII output

void outAscii_signal_map(FILE *file, signal_map *kMap, int field)
{
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  map_t type = kMap->type;
  double *value1 = kMap->value1;
  double *value2 = kMap->value2;
  double theta_pix = kMap->theta_pix;
  
  if (type == proj_mask) fprintf(file, "# 0 = activated, 1 = masked\n");
  fprintf(file, "# Number of pixels = %d\n", kMap->length);
  fprintf(file, "#\n");
  
  double pos[2];
  int i, j, jN1, index;
  
  if (field < 2) {
    if (type == proj_mask) {
      fprintf(file, "#  theta_x   theta_y  mask\n");
      fprintf(file, "# [arcmin]  [arcmin]   [-]\n");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fprintf(file, "  %8.3f  %8.3f  %4.0f\n", pos[0]*RADIAN_TO_ARCMIN, pos[1]*RADIAN_TO_ARCMIN, value1[i+jN1]);
	}
      }
    }
    
    else if (type < 6) {
      fprintf(file, "#  theta_x   theta_y        value\n");
      fprintf(file, "# [arcmin]  [arcmin]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fprintf(file, "  %8.3f  %8.3f  %11.4e\n", pos[0]*RADIAN_TO_ARCMIN, pos[1]*RADIAN_TO_ARCMIN, value1[i+jN1]);
	}
      }
    }
      
    else {
      fprintf(file, "#  theta_x   theta_y       value1       value2\n");
      fprintf(file, "# [arcmin]  [arcmin]          [-]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  index = i + jN1;
	  getPixPos(pos, theta_pix, i, j);
	  fprintf(file, "  %8.3f  %8.3f  %11.4e  %11.4e\n", pos[0]*RADIAN_TO_ARCMIN, pos[1]*RADIAN_TO_ARCMIN, value1[index], value2[index]);
	}
      }
    }
  }
  
  else {
    if (type < 6) {
      fprintf(file, "#   i      j          value\n");
      fprintf(file, "# [pix]  [pix]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) fprintf(file, "  %5d  %5d  %11.4e\n", i, j, value1[i+jN1]);
      }
    }
      
    else {
      fprintf(file, "#   i      j         value1       value2\n");
      fprintf(file, "# [pix]  [pix]          [-]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  index = i + jN1;
	  fprintf(file, "  %5d  %5d  %11.4e  %11.4e\n", i, j, value1[index], value2[index]);
	}
      }
    }
  }
  
  return;
}

void outAscii_fftw_complex(FILE *file, peak_param *pkPar, fftw_complex *table, map_t type)
{
  int N1 = (pkPar->field < 2) ? pkPar->resol[0] : pkPar->HP_resol;
  int N2 = (pkPar->field < 2) ? pkPar->resol[1] : pkPar->HP_resol;
  int M  = pkPar->FFTSize;
  int field = pkPar->field;
  double theta_pix = pkPar->theta_pix;
  
  double *pix, pos[2];
  int i, j, jM;
  
  if (type == proj_mask) fprintf(file, "# 0 = activated, 1 = masked\n");
  fprintf(file, "# Number of pixels = %d\n", N1 * N2);
  fprintf(file, "#\n");
  
  if (pkPar->field < 2) {
    if (type == proj_mask) {
      fprintf(file, "#  theta_x   theta_y  mask\n");
      fprintf(file, "# [arcmin]  [arcmin]   [-]\n");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fprintf(file, "  %8.3f  %8.3f  %4.0f\n", pos[0]*RADIAN_TO_ARCMIN, pos[1]*RADIAN_TO_ARCMIN, table[i+jM][0]);
	}
      }
    }
    
    else if (type < 6) {
      fprintf(file, "#  theta_x   theta_y        value\n");
      fprintf(file, "# [arcmin]  [arcmin]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fprintf(file, "  %8.3f  %8.3f  %11.4e\n", pos[0]*RADIAN_TO_ARCMIN, pos[1]*RADIAN_TO_ARCMIN, table[i+jM][0]);
	}
      }
    }
    
    else {
      fprintf(file, "#  theta_x   theta_y       value1       value2\n");
      fprintf(file, "# [arcmin]  [arcmin]          [-]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  pix = table[i+jM];
	  getPixPos(pos, theta_pix, i, j);
	  fprintf(file, "  %8.3f  %8.3f  %11.4e  %11.4e\n", pos[0]*RADIAN_TO_ARCMIN, pos[1]*RADIAN_TO_ARCMIN, pix[0], pix[1]);
	}
      }
    }
  }
  
  else {
    if (type < 6) {
      fprintf(file, "#    i      j         value\n");
      fprintf(file, "# [pix]  [pix]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  fprintf(file, "  %5d  %5d  %11.4e\n", i, j, table[i+jM][0]);
	}
      }
    }
    
    else {
      fprintf(file, "#    i      j        value1       value2\n");
      fprintf(file, "# [pix]  [pix]          [-]          [-]\n");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  pix = table[i+jM];
	  fprintf(file, "  %5d  %5d  %11.4e  %11.4e\n", i, j, pix[0], pix[1]);
	}
      }
    }
  }
  
  return;
}

void outAsciiNoiseInfo(FILE *file, peak_param *pkPar)
{
  if (pkPar->doNoise) fprintf(file, "# n_gal = %g [arcmin^-2], sigma_eps = %g, sigma_pix = %g\n", pkPar->n_gal/RADIAN_SQ_TO_ARCMIN_SQ, pkPar->sigma_eps, pkPar->sigma_pix);
  return;
}

void outAsciiFilterInfo(FILE *file, peak_param *pkPar, map_t type, int filterInd)
{
  if (type % 2 == 0) return;
  
  fprintf(file, "# doSmoothing = %d\n", pkPar->doSmoothing);
  
  int W = pkPar->FFT_nbFilters;
  int j = filterInd;
  
  if (j < W) {
    if (pkPar->filter[j] != 3) 
      fprintf(file, "# FFT filter = %s, FFT scale = %g [arcmin] = %g [pix]\n", 
	      printFilter(pkPar->FFT_filter[j]), pkPar->FFT_scale[j]*RADIAN_TO_ARCMIN, pkPar->FFT_scaleInPix[j]);
    else 
      fprintf(file, "# FFT filter = %s, FFT scale = %g, %g [arcmin] = %g, %g [pix]\n", 
	      printFilter(pkPar->FFT_filter[j]), pkPar->FFT_scale[j]*RADIAN_TO_ARCMIN, pkPar->FFT_scale[j+1]*RADIAN_TO_ARCMIN, pkPar->FFT_scaleInPix[j], pkPar->FFT_scaleInPix[j+1]);
  }
  else {
    if (pkPar->filter[j] != 3) 
      fprintf(file, "# DC filter = %s, DC scale = %g [arcmin] = %g [pix]\n", 
	      printFilter(pkPar->DC_filter[j-W]), pkPar->DC_scale[j-W]*RADIAN_TO_ARCMIN, pkPar->DC_scaleInPix[j-W]);
    else 
      fprintf(file, "# DC filter = %s, DC scale = %g, %g [arcmin] = %g, %g [pix]\n", 
	      printFilter(pkPar->DC_filter[j-W]), pkPar->DC_scale[j-W]*RADIAN_TO_ARCMIN, pkPar->DC_scale[j+1-W]*RADIAN_TO_ARCMIN, pkPar->DC_scaleInPix[j-W], pkPar->DC_scaleInPix[j+1-W]);
  }
  return;
}

void outAsciiAllFilterInfo(FILE *file, peak_param *pkPar)
{
  int FFT_nbFilters = pkPar->FFT_nbFilters;
  int *FFT_filter   = pkPar->FFT_filter;
  double *FFT_scale = pkPar->FFT_scale;
  int DC_nbFilters  = pkPar->DC_nbFilters;
  int *DC_filter    = pkPar->DC_filter;
  double *DC_scale  = pkPar->DC_scale;
  int i;
  
  fprintf(file, "Number of filters = %d\n", pkPar->nbFilters);
  
  if (FFT_nbFilters > 0) {
    fprintf(file, "FFT filters = %s", printFilter(FFT_filter[0]));
    for (i=1; i<FFT_nbFilters; i++) fprintf(file, ", %s", printFilter(FFT_filter[i]));
    fprintf(file, "\n");
    
    fprintf(file, "FFT scales = %.3f", FFT_scale[0]*RADIAN_TO_ARCMIN);
    for (i=1; i<FFT_nbFilters+pkPar->FFT_hasGammaT; i++) fprintf(file, ", %.3f", FFT_scale[i]*RADIAN_TO_ARCMIN);
    fprintf(file, " [arcmin]\n");
  }
  
  if (DC_nbFilters > 0) {
    fprintf(file, "DC filters = %s", printFilter(DC_filter[0]));
    for (i=1; i<DC_nbFilters; i++) fprintf(file, ", %s", printFilter(DC_filter[i]));
    fprintf(file, "\n");
  
    fprintf(file, "FFT scales = %.3f", DC_scale[0]*RADIAN_TO_ARCMIN);
    for (i=1; i<DC_nbFilters+pkPar->DC_hasGammaT; i++) fprintf(file, ", %.3f", DC_scale[i]*RADIAN_TO_ARCMIN);
    fprintf(file, " [arcmin]\n");
  }
  return;
}

void outAsciiMap(char name[], cosmo_hm *chPar, peak_param *pkPar, signal_map *kMap, int filterInd, error **err)
{
  if (pkPar->outMaps == 0 && pkPar->outTruth == 0) return;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# %s\n", STR_MAP_T(kMap->type));
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiLensingInfo(file, pkPar);
  outAsciiNoiseInfo(file, pkPar);
  outAsciiFilterInfo(file, pkPar, kMap->type, filterInd);
  fprintf(file, "#\n");
  
  outAscii_signal_map(file, kMap, pkPar->field);
  
  fclose(file);
  if (pkPar->verbose < 3)  printf("Outputed \"%s\"\n", name);
  return;
}

void outAsciiMapFromTable(char name[], cosmo_hm *chPar, peak_param *pkPar, fftw_complex *table, map_t type, int filterInd, error **err)
{
  if (pkPar->outMaps == 0 && pkPar->outTruth == 0) return;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# %s\n", STR_MAP_T(type));
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiLensingInfo(file, pkPar);
  outAsciiNoiseInfo(file, pkPar);
  outAsciiFilterInfo(file, pkPar, type, filterInd);
  fprintf(file, "#\n");
  
  outAscii_fftw_complex(file, pkPar, table, type);
  
  fclose(file);
  if (pkPar->verbose < 3)  printf("Outputed \"%s\"\n", name);
  return;
}

void outAsciiMask(char name[], peak_param *pkPar, gal_map *gMap, signal_map *kMap, error **err)
{
  if (pkPar->outMask == 0) return;
  
  kMap->type = proj_mask;
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# %s\n", STR_MAP_T(kMap->type));
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  fprintf(file, "#\n");
  
  setFillingThreshold(pkPar, gMap, err); forwardError(*err, __LINE__,);
  
  double threshold = gMap->fillingThreshold;
  gal_list **map = gMap->map;
  int i;
  
  for (i=0; i<gMap->length; i++) {
    if (map[i]->totWeight < threshold) kMap->value1[i] = 1.0; //-- Masked
    else                               kMap->value1[i] = 0.0;
  }
  
  outAscii_signal_map(file, kMap, 0); //-- field = 0
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to FITS output

#ifdef __CAMELUS_USE_FITS__
void outFits_signal_map(FITS_t *fits, signal_map *kMap, int field, double factor)
{
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  map_t type = kMap->type;
  double *value1 = kMap->value1;
  double *value2 = kMap->value2;
  double theta_pix = kMap->theta_pix;
  
  double pos[2];
  float fBuff;
  short hBuff;
  int i, j, jN1, index;
  
  if (field < 2) {
    if (type == proj_mask) {
      addColumn(fits, "THX",  TFLOAT, "arcmin");
      addColumn(fits, "THY",  TFLOAT, "arcmin");
      addColumn(fits, "MASK", TSHORT, "-       ");
      updateComment(fits, "TTYPE3", "0 = activated, 1 = masked");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fBuff = (float)(pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
	  fBuff = (float)(pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
	  hBuff = (short)value1[i+jN1];     writeTableColumn(fits, 2, 1, &hBuff);
	  nextRow(fits);
	}
      }
    }
    
    else if (type < 6) {
      addColumn(fits, "THX",   TFLOAT, "arcmin");
      addColumn(fits, "THY",   TFLOAT, "arcmin");
      addColumn(fits, "VALUE", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fBuff = (float)(pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
	  fBuff = (float)(pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
	  fBuff = (float)value1[i+jN1];     writeTableColumn(fits, 2, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
    
    else {
      addColumn(fits, "THX",    TFLOAT, "arcmin");
      addColumn(fits, "THY",    TFLOAT, "arcmin");
      addColumn(fits, "VALUE1", TFLOAT, "-       ");
      addColumn(fits, "VALUE2", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  index = i + jN1;
	  getPixPos(pos, theta_pix, i, j);
	  fBuff = (float)(pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
	  fBuff = (float)(pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
	  fBuff = (float)value1[index];     writeTableColumn(fits, 2, 1, &fBuff);
	  fBuff = (float)value2[index];     writeTableColumn(fits, 3, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
  }
  
  else {
    if (type < 6) {
      addColumn(fits, "IPIX",  TINT, "pix");
      addColumn(fits, "JPIX",  TINT, "pix");
      addColumn(fits, "VALUE", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	                                     writeTableColumn(fits, 0, 1, &i);
	                                     writeTableColumn(fits, 1, 1, &j);
	  fBuff = (float)value1[i+jN1]; writeTableColumn(fits, 2, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
      
    else {
      addColumn(fits, "IPIX",   TINT, "pix");
      addColumn(fits, "JPIX",   TINT, "pix");
      addColumn(fits, "VALUE1", TFLOAT, "-       ");
      addColumn(fits, "VALUE2", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jN1 = j * N1;
	for (i=0; i<N1; i++) {
	  index = i + jN1;
	                                     writeTableColumn(fits, 0, 1, &i);
	                                     writeTableColumn(fits, 1, 1, &j);
	  fBuff = (float)value1[index]; writeTableColumn(fits, 2, 1, &fBuff);
	  fBuff = (float)value2[index]; writeTableColumn(fits, 3, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
  }
  
  addLineSpread(fits);
  addKeyword(fits, TINT, "NBPIX", &kMap->length, "[-] Number of pixels");
  return;
}

void outFits_fftw_complex(FITS_t *fits, peak_param *pkPar, fftw_complex *table, int type)
{
  int N1 = (pkPar->field < 2) ? pkPar->resol[0] : pkPar->HP_resol;
  int N2 = (pkPar->field < 2) ? pkPar->resol[1] : pkPar->HP_resol;
  int M  = pkPar->FFTSize;
  int field = pkPar->field;
  double theta_pix = pkPar->theta_pix;
  double factor = (field == 0) ? RADIAN_TO_ARCMIN : RADIAN_TO_DEGREE;
  
  double *pix, pos[2];
  float fBuff;
  int i, j, jM;
  
  if (field < 2) {
    if (type < 6) {
      addColumn(fits, "THX",   TFLOAT, "arcmin");
      addColumn(fits, "THY",   TFLOAT, "arcmin");
      addColumn(fits, "VALUE", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  getPixPos(pos, theta_pix, i, j);
	  fBuff = (float)(pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
	  fBuff = (float)(pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
	  fBuff = (float)table[i+jM][0];    writeTableColumn(fits, 2, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
    
    else {
      addColumn(fits, "THX",    TFLOAT, "arcmin");
      addColumn(fits, "THY",    TFLOAT, "arcmin");
      addColumn(fits, "VALUE1", TFLOAT, "-       ");
      addColumn(fits, "VALUE2", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  pix = table[i+jM];
	  getPixPos(pos, theta_pix, i, j);
	  fBuff = (float)(pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
	  fBuff = (float)(pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
	  fBuff = (float)pix[0];            writeTableColumn(fits, 2, 1, &fBuff);
	  fBuff = (float)pix[1];            writeTableColumn(fits, 3, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
  }
  
  else {
    if (type < 6) {
      addColumn(fits, "IPIX",  TINT, "pix");
      addColumn(fits, "JPIX",  TINT, "pix");
      addColumn(fits, "VALUE", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
					 writeTableColumn(fits, 0, 1, &i);
					 writeTableColumn(fits, 1, 1, &j);
	  fBuff = (float)table[i+jM][0]; writeTableColumn(fits, 2, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
    
    else {
      addColumn(fits, "IPIX",   TINT, "pix");
      addColumn(fits, "JPIX",   TINT, "pix");
      addColumn(fits, "VALUE1", TFLOAT, "-       ");
      addColumn(fits, "VALUE2", TFLOAT, "-       ");
      
      for (j=0; j<N2; j++) {
	jM = j * M;
	for (i=0; i<N1; i++) {
	  pix = table[i+jM];
				 writeTableColumn(fits, 0, 1, &i);
				 writeTableColumn(fits, 1, 1, &j);
	  fBuff = (float)pix[0]; writeTableColumn(fits, 2, 1, &fBuff);
	  fBuff = (float)pix[1]; writeTableColumn(fits, 3, 1, &fBuff);
	  nextRow(fits);
	}
      }
    }
  }
  
  int nbPix = N1 * N2;
  addLineSpread(fits);
  addKeyword(fits, TINT, "NBPIX", &nbPix, "[-] Number of pixels");
  return;
}

void outFitsNoiseInfo(FITS_t *fits, peak_param *pkPar)
{
  addLineSpread(fits);
  
  double lfBuff;
  
  lfBuff = pkPar->n_gal / RADIAN_SQ_TO_ARCMIN_SQ;
  addKeyword(fits, TDOUBLE, "NGAL",     &lfBuff,           "[arcmin^-2] Galaxy number density");
  addKeyword(fits, TDOUBLE, "SIGMAEPS", &pkPar->sigma_eps, "[-] Ellipticity dispersion");
  lfBuff = pkPar->theta_pix * RADIAN_TO_ARCMIN;
  addKeyword(fits, TDOUBLE, "THPIX",    &lfBuff,           "[arcmin] Pixel size");
  return;
}

void outFitsFilterInfo(FITS_t *fits, peak_param *pkPar, map_t type, int filterInd)
{
  if (type % 2 == 0) return;
  
  addLineSpread(fits);
  addKeyword(fits, TINT,    "DOSMOOTH", &pkPar->doSmoothing,          "[-] 0 = without, 1 = binning and FFT, 2 = direct convolution");
  addComment(fits,                                                    "Sum for performing simultaneously multiple techniques, e.g. 3 = FFT + DC");
  
  int W = pkPar->FFT_nbFilters;
  int j = filterInd;
  double lfBuff;
  
  if (j < W) {
    if (pkPar->filter[j] != 3) {
      addKeyword(fits, TINT,    "FFTFILT",  &pkPar->FFT_filter[j],        "[-] 0 = Gaussian, 1 = starlet, 2 = M_ap tanh, 3 = M_ap gamma_t");
      lfBuff = pkPar->FFT_scale[j] * RADIAN_TO_ARCMIN;
      addKeyword(fits, TDOUBLE, "FFTSCALE", &lfBuff,                      "[arcmin] Filter size for FFT");
      addKeyword(fits, TDOUBLE, "FFTPIX",   &pkPar->FFT_scaleInPix[j],    "[pix] Filter size for FFT");
    }
    else {
      addKeyword(fits, TINT,    "FFTFILT",  &pkPar->FFT_filter[j],        "[-] 0 = Gaussian, 1 = starlet, 2 = M_ap tanh, 3 = M_ap gamma_t");
      lfBuff = pkPar->FFT_scale[j] * RADIAN_TO_ARCMIN;
      addKeyword(fits, TDOUBLE, "FFTSCALE", &lfBuff,                      "[arcmin] Filter size for FFT");
      lfBuff = pkPar->FFT_scale[j+1] * RADIAN_TO_ARCMIN;
      addKeyword(fits, TDOUBLE, "FFTSCAL2", &lfBuff,                      "[arcmin] Filter size for FFT");
      addKeyword(fits, TDOUBLE, "FFTPIX",   &pkPar->FFT_scaleInPix[j],    "[pix] Filter size for FFT");
      addKeyword(fits, TDOUBLE, "FFTPIX2",  &pkPar->FFT_scaleInPix[j+1],  "[pix] Filter size for FFT");
    }
  }
  else {
    if (pkPar->filter[j] != 3) {
      addKeyword(fits, TINT,    "FFTFILT",  &pkPar->DC_filter[j-W],       "[-] 0 = Gaussian, 1 = starlet, 2 = M_ap tanh, 3 = M_ap gamma_t");
      lfBuff = pkPar->DC_scale[j-W] * RADIAN_TO_ARCMIN;
      addKeyword(fits, TDOUBLE, "FFTSCALE", &lfBuff,                      "[arcmin] Filter size for FFT");
      addKeyword(fits, TDOUBLE, "FFTPIX",   &pkPar->DC_scaleInPix[j-W],   "[pix] Filter size for FFT");
    }
    else {
      addKeyword(fits, TINT,    "FFTFILT",  &pkPar->DC_filter[j-W],       "[-] 0 = Gaussian, 1 = starlet, 2 = M_ap tanh, 3 = M_ap gamma_t");
      lfBuff = pkPar->DC_scale[j-W] * RADIAN_TO_ARCMIN;
      addKeyword(fits, TDOUBLE, "FFTSCALE", &lfBuff,                      "[arcmin] Filter size for FFT");
      lfBuff = pkPar->DC_scale[j+1-W] * RADIAN_TO_ARCMIN;
      addKeyword(fits, TDOUBLE, "FFTSCAL2", &lfBuff,                      "[arcmin] Filter size for FFT");
      addKeyword(fits, TDOUBLE, "FFTPIX",   &pkPar->DC_scaleInPix[j-W],   "[pix] Filter size for FFT");
      addKeyword(fits, TDOUBLE, "FFTPIX2",  &pkPar->DC_scaleInPix[j+1-W], "[pix] Filter size for FFT");
    }
  }
  return;
}

void outFitsAllFilterInfo(FITS_t *fits, peak_param *pkPar)
{
  int FFT_nbFilters = pkPar->FFT_nbFilters;
  int *FFT_filter   = pkPar->FFT_filter;
  double *FFT_scale = pkPar->FFT_scale;
  int DC_nbFilters  = pkPar->DC_nbFilters;
  int *DC_filter    = pkPar->DC_filter;
  double *DC_scale  = pkPar->DC_scale;
  
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  int i;
  
  addKeyword(fits, TINT,   "NBFILT", &pkPar->nbFilters, "[-] Number of filters");
  
  if (FFT_nbFilters > 0) {
    for (i=0; i<FFT_nbFilters; i++) {
      sprintf(name, "FFT%d", i);
      if (FFT_filter[i] < 3) {
	sprintf(name2, "%d, %.3f", FFT_filter[i], FFT_scale[i]);
	addKeyword(fits, TSTRING, name,    name2,             "[-, arcmin] Type, scale");
      }
      else {
	sprintf(name2, "%d, %.3f, %.3f", FFT_filter[i], FFT_scale[i], FFT_scale[i+1]);
	addKeyword(fits, TSTRING, name,    name2,             "[-, arcmin] Type, scale");
      }
    }
  }
  
  if (DC_nbFilters > 0) {
    for (i=0; i<FFT_nbFilters; i++) {
      sprintf(name, "DC%d", i);
      if (DC_filter[i] < 3) {
	sprintf(name2, "%d, %.3f", DC_filter[i], DC_scale[i]);
	addKeyword(fits, TSTRING, name,    name2,             "[-, arcmin] Type, scale");
      }
      else {
	sprintf(name2, "%d, %.3f, %.3f", DC_filter[i], DC_scale[i], DC_scale[i+1]);
	addKeyword(fits, TSTRING, name,    name2,             "[-, arcmin] Type, scale");
      }
    }
  }
  return;
}
#endif

void outFitsMap(char name[], cosmo_hm *chPar, peak_param *pkPar, signal_map *kMap, int filterInd)
{
  if (pkPar->outMaps == 0 && pkPar->outTruth == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  double factor = (pkPar->field == 0) ? RADIAN_TO_ARCMIN : RADIAN_TO_DEGREE;
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  outFits_signal_map(fits, kMap, pkPar->field, factor);
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsLensingInfo(fits, pkPar);
  outFitsNoiseInfo(fits, pkPar);
  outFitsFilterInfo(fits, pkPar, kMap->type, filterInd);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

void outFitsMapFromTable(char name[], cosmo_hm *chPar, peak_param *pkPar, fftw_complex *table, map_t type, int filterInd)
{
  if (pkPar->outMaps == 0 && pkPar->outTruth == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  outFits_fftw_complex(fits, pkPar, table, type);
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsLensingInfo(fits, pkPar);
  outFitsNoiseInfo(fits, pkPar);
  outFitsFilterInfo(fits, pkPar, type, filterInd);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3)  printf("Outputed \"%s\"\n", name);
#endif
  return;
}

void outputMapFromTable(char name[], cosmo_hm *chPar, peak_param *pkPar, fftw_complex *table, map_t type, int filterInd, error **err)
{
  char name2[STRING_LENGTH_MAX];
  if (pkPar->doFITS == 0) {
    outAsciiMapFromTable(name, chPar, pkPar, table, type, filterInd, err);
    forwardError(*err, __LINE__,);
  }
  else {
    sprintf(name2, "%s.fits", name);
    outFitsMapFromTable(name2, chPar, pkPar, table, type, filterInd);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to map making

void gMapToSmoother(gal_map *gMap, FFT_t *smoo)
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

void galaxyBinning(peak_param *pkPar, gal_map *gMap, gal_map *gMap2, FFT_t *smoo, error **err)
{
  if (pkPar->field < 2) {
    signalMean_gal_map(pkPar, gMap, err);  forwardError(*err, __LINE__,);
    gMapToSmoother(gMap, smoo);
  }
  else {
    rebinGalaxiesForHEALPix(pkPar, gMap, gMap2, err);
    signalMean_gal_map(pkPar, gMap2, err); forwardError(*err, __LINE__,);
    gMapToSmoother(gMap2, smoo);
  }
  if (pkPar->verbose < 3) printf("Binned galaxies\n");
  return;
}

void invertForTrueMap(peak_param *pkPar, FFT_t *smoo, signal_map *kMap)
{
  if (pkPar->doKappa == 3) {
    invertByIterKS(smoo, kMap);
    zeroPadding(smoo->before, kMap->N1, kMap->N2, smoo->N); //-- 0-padding
    if (pkPar->verbose < 3) printf("Inverted by iterative KS\n");
  }
  else if (pkPar->doKappa == 4) {
    invertBySS(smoo, kMap);
    zeroPadding(smoo->before, kMap->N1, kMap->N2, smoo->N); //-- 0-padding
    if (pkPar->verbose < 3) printf("Inverted by SS\n");
  }
  else {
    invertByLinKS(smoo);
    zeroPadding(smoo->before, kMap->N1, kMap->N2, smoo->N); //-- 0-padding
    if (pkPar->verbose < 3) printf("Inverted by linear KS\n");
  }
  return;
}

void pixelization(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, error **err)
{
  if (pkPar->DC_nbFilters) smoothAndBinByDC(pkPar, gMap, DCSmoother); //-- smoother reset in this function
  
  if (pkPar->FFT_nbFilters) {
    //-- Bin and stock information in FFTSmoother->array[0]->before
    //-- Field = 2 already excluded
    galaxyBinning(pkPar, gMap, NULL, FFTSmoother->array[0], err); forwardError(*err, __LINE__,); //-- gMap2 = NULL
  }
  return;
}

void inversion(peak_param *pkPar, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap)
{
  int k = pkPar->FFT_firstToInvert;
  int M = FFTSmoother->array[0]->N;
  int print = 0;
  int FFTLength = FFTSmoother->array[0]->length;
  int FFT_nbFilters = pkPar->FFT_nbFilters;
  fftw_complex *smoo_before = FFTSmoother->array[0]->before;
  FFT_t *smoo;
  int i;
  
  //-- In practice, this part is never executed, because DC_filter = 1 or 2 is forbidden with doKappa != 1.
  for (i=0; i<pkPar->DC_nbFilters; i++) {
    if (pkPar->DC_filter[i] < 2) {
      if (pkPar->doKappa == 1) continue;
      if      (pkPar->doKappa == 3) invertByIterKS(DCSmoother->array[i], kMap);
      else if (pkPar->doKappa == 4) invertBySS(DCSmoother->array[i], kMap);
      else                          invertByLinKS(DCSmoother->array[i]);
      print = 1;
    }
  }
  
  if (FFT_nbFilters) {
    //-- If doKappa = 1, copy kappa to all tables
    if (pkPar->doKappa == 1) {
      for (i=1; i<FFT_nbFilters; i++) copy_fftw_complex(smoo_before, FFTSmoother->array[i]->before, FFTLength);
      return; //-- No inversion to do
    }
    
    //-- Copy gamma to M_ap tables
    for (i=1; i<FFT_nbFilters; i++) {
      if (pkPar->FFT_filter[i] >= 2) copy_fftw_complex(smoo_before, FFTSmoother->array[i]->before, FFTLength);
    }
    
    if (k == FFT_nbFilters) return; //-- No inversion to do
    if (k > 0) copy_fftw_complex(smoo_before, FFTSmoother->array[k]->before, FFTLength); //-- Copy to first to invert
    
    //-- Invert
    smoo = FFTSmoother->array[k];
    if      (pkPar->doKappa == 3) invertByIterKS(smoo, kMap);
    else if (pkPar->doKappa == 4) invertBySS(smoo, kMap);
    else                          invertByLinKS(smoo);
    zeroPadding(smoo->before, kMap->N1, kMap->N2, M); //-- 0-padding
    print = 1;
    
    //-- Copy kappa to the rest
    for (i=k; i<FFT_nbFilters; i++) {
      if (pkPar->FFT_filter[i] < 2) copy_fftw_complex(smoo->before, FFTSmoother->array[i]->before, FFTLength);
    }
  }
  
  if (pkPar->verbose < 3 && print) {
    if      (pkPar->doKappa == 3) printf("Inverted by iterative KS\n");
    else if (pkPar->doKappa == 4) printf("Inverted by SS\n");
    else                          printf("Inverted by linear KS\n");
  }
  return;
}

void makeMaps(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap, error **err)
{
  //-- Map making main function: kappa/gamma/g, noiseless/noisy, unsmoothed/smoothed
    
  //-- Binning + DC
  pixelization(pkPar, gMap, FFTSmoother, DCSmoother, err); forwardError(*err, __LINE__,);
  
  //-- Inversion
  inversion(pkPar, FFTSmoother, DCSmoother, kMap);
  
  //-- FFT smoothing
  if (pkPar->FFT_nbFilters) smoothByFFT_arr(pkPar, FFTSmoother);
  return;
}

void makeMapsAndOutput(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, gal_map *gMap2, FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap, error **err)
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
  //-- Bin g, linear KS, FFT (3, 1)                    | DC, iterative KS (3, 2)
  //--   Catalogue:     g              in gMap         |   Catalogue:          g              in gMap
  //--   Binning:       binned g       in smoo->before |   Direct convolution: smoothed g     in smoo->before
  //--   Linear KS:     kappa          in smoo->before |   Linear KS:          smoothed kappa in smoo->before
  //--   FFT smoothing: smoothed kappa in smoo->after  |
  //--                                                 |
  //-- Bin g, iterative KS, FFT (2, 1)                 | DC, iterative KS (2, 2)
  //--   Catalogue:     g              in gMap         |   Catalogue:          g              in gMap
  //--   Binning:       binned g       in smoo->before |   Direct convolution: smoothed g     in smoo->before
  //--   Iterative KS:  kappa          in smoo->before |   Iterative KS:       smoothed kappa in smoo->before
  //--   FFT smoothing: smoothed kappa in smoo->after  |
  
  int doNoise       = pkPar->doNoise;
  int doKappa       = pkPar->doKappa;
  int FFT_nbFilters = pkPar->FFT_nbFilters;
  int DC_nbFilters  = pkPar->DC_nbFilters;
  int typeForNoise  = doNoise * 2;
  int typeForG      = (int)(doKappa >= 2) * 6;
  FFT_t *smoo       = FFTSmoother->array[0];
  
  reset_fftw_complex(smoo->before, smoo->length);
  char name[STRING_LENGTH_MAX];
  int j;
  
  if (pkPar->outTruth == 1) {
    //-- True gamma or g
    if (doKappa != 1) {
      galaxyBinning(pkPar, gMap, gMap2, smoo, err); forwardError(*err, __LINE__,);
      sprintf(name, "%sgammaOrGMap_truth", pkPar->prefix);
      outputMapFromTable(name, chPar, pkPar, smoo->before, gamma_map+typeForG, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
      
      if (pkPar->field < 2) {
	invertForTrueMap(pkPar, smoo, kMap);
	sprintf(name, "%skappaMap_semiTruth", pkPar->prefix);
	outputMapFromTable(name, chPar, pkPar, smoo->before, kappa_map, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
      }
    }
    
    //-- True kappa
    if (doKappa > 0) {
      pkPar->doKappa = 1;
      galaxyBinning(pkPar, gMap, gMap2, smoo, err); forwardError(*err, __LINE__,);
      pkPar->doKappa = doKappa;
      sprintf(name, "%skappaMap_truth", pkPar->prefix);
      outputMapFromTable(name, chPar, pkPar, smoo->before, kappa_map, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
    }
  }
  
  if (doNoise && pkPar->doLensing) {
    //-- Add noise
    addNoiseToGalaxies(pkPar, gMap);
    if (pkPar->doFITS == 0) {
      sprintf(name, "%sgalCat_noisy", pkPar->prefix);
      outAsciiGalCat(name, chPar, pkPar, gMap, err); forwardError(*err, __LINE__,);
    }
    else {
      sprintf(name, "%sgalCat_noisy.fits", pkPar->prefix);
      outFitsGalCat(name, chPar, pkPar, gMap);
    }
    
    //-- Unsmoothed map for DC-only case
    if (FFT_nbFilters == 0 && DC_nbFilters > 0) {
      galaxyBinning(pkPar, gMap, gMap2, smoo, err); forwardError(*err, __LINE__,);
      if (doKappa == 1) {
	sprintf(name, "%skappaMap_unsmoothed", pkPar->prefix);
	outputMapFromTable(name, chPar, pkPar, smoo->before, kappa_map+typeForNoise, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
      }
      else {
	sprintf(name, "%sgammaOrGMap_unsmoothed", pkPar->prefix);
	outputMapFromTable(name, chPar, pkPar, smoo->before, gamma_map+typeForG+typeForNoise, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
      }
    }
  }
  
  //-- Mask
  sprintf(name, "%smask", pkPar->prefix);
  outAsciiMask(name, pkPar, gMap, kMap, err); forwardError(*err, __LINE__,);
  
  //-- Pixelization
  pixelization(pkPar, gMap, FFTSmoother, DCSmoother, err); forwardError(*err, __LINE__,);
  if (doKappa != 1) {
    if (doNoise && FFT_nbFilters) {
      sprintf(name, "%sgammaOrGMap_unsmoothed", pkPar->prefix);
      outputMapFromTable(name, chPar, pkPar, smoo->before, gamma_map+typeForG+typeForNoise, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
    }
    for (j=0; j<DC_nbFilters; j++) {
      if (pkPar->DC_filter[j] < 2) {
	sprintf(name, "%sgammaOrGMap_DC%d", pkPar->prefix, j);
	outputMapFromTable(name, chPar, pkPar, DCSmoother->array[j]->before, R_map+typeForG+typeForNoise, j+FFT_nbFilters, err); forwardError(*err, __LINE__,);
      }
    }
  }
  
  //-- Inversion
  inversion(pkPar, FFTSmoother, DCSmoother, kMap);
  if (doNoise && FFT_nbFilters && pkPar->FFT_firstToInvert < pkPar->FFT_nbFilters) {
    sprintf(name, "%skappaMap_unsmoothed", pkPar->prefix);
    outputMapFromTable(name, chPar, pkPar, smoo->before, kappa_map+typeForNoise, 0, err); forwardError(*err, __LINE__,); //-- filterInd = 0
  }
  for (j=0; j<DC_nbFilters; j++) {
    sprintf(name, "%skappaMap_DC%d", pkPar->prefix, j);
    outputMapFromTable(name, chPar, pkPar, DCSmoother->array[j]->before, K_map+typeForNoise, j+FFT_nbFilters, err); forwardError(*err, __LINE__,);
  }
  
  //-- FFT smoothing
  if (FFT_nbFilters) smoothByFFT_arr(pkPar, FFTSmoother);
  for (j=0; j<FFT_nbFilters; j++) {
    sprintf(name, "%skappaMap_FFT%d", pkPar->prefix, j);
    outputMapFromTable(name, chPar, pkPar, FFTSmoother->array[j]->after, K_map+typeForNoise, j, err); forwardError(*err, __LINE__,);
  }
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doKMap(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  sampler_arr *hSampArr   = initialize_sampler_arr(pkPar->N_z_halo, pkPar->N_M+1);
  interpolator_t *k1Inter = initialize_interpolator_t(pkPar->N_z_gal+1);
  setMassSampAndK1Inter(chPar, pkPar, hSampArr, NULL, k1Inter, err);                                        forwardError(*err, __LINE__,); //-- lambda = NULL
  halo_map *hMap          = initialize_halo_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);   forwardError(*err, __LINE__,);
  sampler_t *gSamp        = initialize_sampler_t(pkPar->N_z_gal+1);
  setGalaxySampler(chPar, pkPar, gSamp, err);                                                               forwardError(*err, __LINE__,);
  gal_map *gMap           = initialize_gal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);    forwardError(*err, __LINE__,);
  gal_map *gMap2          = initialize_gal_map(pkPar->HP_resol, pkPar->HP_resol, pkPar->theta_pix, err);    forwardError(*err, __LINE__,);
  mask_map *mask          = initializeMask(pkPar, err);                                                     forwardError(*err, __LINE__,);
  FFT_arr *FFTSmoother    = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  makeKernel(pkPar, FFTSmoother);
  FFT_arr *DCSmoother     = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  signal_map *kMap        = initialize_signal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err); forwardError(*err, __LINE__,);
  
  readCatOrMakeSimulAndOutput(chPar, pkPar, hSampArr, hMap, err);                   forwardError(*err, __LINE__,);
  cleanOrMakeOrResample(chPar, pkPar, gSamp, gMap, mask, err);                      forwardError(*err, __LINE__,);
  lensingCatalogueAndOutput(chPar, pkPar, hMap, gMap, k1Inter, err);                forwardError(*err, __LINE__,);
  makeMapsAndOutput(chPar, pkPar, gMap, gMap2, FFTSmoother, DCSmoother, kMap, err); forwardError(*err, __LINE__,);
  
  free_sampler_arr(hSampArr);
  free_interpolator_t(k1Inter);
  free_halo_map(hMap);
  free_sampler_t(gSamp);
  free_gal_map(gMap);
  free_gal_map(gMap2);
  free_mask_map(mask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_signal_map(kMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

