

  /*******************************************************
   **  peakSelection.c					**
   **  Version 2018.03.11				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "peakSelection.h"


//----------------------------------------------------------------------
//-- Functions related to local variance

void makeKernelForVariance(peak_param *pkPar, FFT_arr *variance)
{
  FFT_t *var;
  int i;
  for (i=0; i<pkPar->FFT_nbFilters; i++) {
    var = variance->array[i];
    if      (pkPar->FFT_filter[i] == 0) fillGaussianKernel(var->kernel, var->length, var->N, pkPar->FFT_scaleInPix[i], 1);  //-- doVar = 1
    else if (pkPar->FFT_filter[i] == 1) fillStarletKernel(var->kernel, var->length, var->N, pkPar->FFT_scaleInPix[i], 1);   //-- doVar = 1
    else if (pkPar->FFT_filter[i] == 2) fillMApTanhKernel(var->kernel, var->length, var->N, pkPar->FFT_scaleInPix[i], 1);   //-- doVar = 1
    else if (pkPar->FFT_filter[i] == 3) fillMApGammaTKernel(var->kernel, var->length, var->N, pkPar->FFT_scaleInPix[i], pkPar->FFT_scaleInPix[i+1], 1); //-- doVar = 1
    fftw_execute(var->kernel_f); //-- To Fourier space
  }
  return;
}

void fillPixelVariance(gal_map *gMap, FFT_t *var)
{
  //-- Fill the shape noise variance for each pixel
  //-- The shape noise variance is the inverse of the total weight: sigma_noise^2 = 1 / totWeight.
  //-- With uniform weighting, the weight of each galaxy is 1 / sigma_half^2 or 2 / sigma_eps^2.
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int M  = var->N;
  double threshold         = gMap->fillingThreshold;
  fftw_complex *var_before = var->before;
  
  double totWeight;
  int i, j, jN1, jM, index_FFT;
  
  //-- Fill the total variance in var->before
  for (j=0; j<M; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<M; i++) {
      index_FFT = i + jM;
      if (i >= N1 || j >= N2) { //-- Buffer area
	var_before[index_FFT][0] = 0.0;
	var_before[index_FFT][1] = 0.0;
      }
      else {
	totWeight = gMap->map[i+jN1]->totWeight;
	var_before[index_FFT][0] = (totWeight < threshold) ? 0.0 :  (1.0 / totWeight); //-- The total variance is the inversed total weight.
	var_before[index_FFT][1] = 0.0;
      }
    }
  }
  
  execute_FFT_t(var);
  return;
}

void makeLocalVariance(peak_param *pkPar, gal_map *gMap, FFT_arr *variance)
{
  //-- Take into account the average by filtering
  
  if (pkPar->FFT_nbFilters == 0) return;
  
  FFT_t *var = variance->array[0];
  fillPixelVariance(gMap, var);
  
  fftw_complex *var_before = var->before;
  int FFTLength = var->length;
  int i;
  
  for (i=1; i<variance->length; i++) {
    var = variance->array[i];
    multiplication_fftw_complex(var_before, var->kernel, var->after, FFTLength); //-- Multiplication
    fftw_execute(var->after_b);                                                  //-- Go to direct space
    rescaleReal_fftw_complex(var->after, FFTLength, pkPar->FFTNormFactor);       //-- Rescale, only real part is interesting.
  }
  return;
}

void kappaToSNR_FFT(peak_param *pkPar, gal_map *gMap, FFT_t *FFTSmoo, signal_map *kMap, FFT_t *var, int FFTScaleInd)
{
  //-- If local noise is used:
  //--   Retrieve kappa from FFTSmoo->after
  //--   Divide it by local std, which is stocked in var->after (SNR = kappa / std_local)
  //--   Stock SNR in kMap->value1
  //--
  //-- If global noise is used:
  //--   Retrieve kappa from FFTSmoo->after
  //--   Divide it by a filter-dependent global value (SNR = kappa / sigma_global)
  //--   Stock SNR in kMap->value1
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int M  = var->N;
  double sigma_noise_inv      = 1.0 / pkPar->FFT_sigma_noise[FFTScaleInd]; //-- Globle noise level
  double threshold            = gMap->fillingThreshold;
  fftw_complex *FFTSmoo_table = FFTSmoo->after;
  fftw_complex *var_after     = var->after; //-- Contain the variance
  double *SNR                 = kMap->value1;
  
  double totWeight;
  int i, j, jN1, jM, index_kMap;
  
  if (pkPar->doLocalNoise) {
    for (j=0; j<N2; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<N1; i++) {
	index_kMap      = i + jN1;
	totWeight       = gMap->map[index_kMap]->totWeight;
	SNR[index_kMap] = (totWeight < threshold) ? -DBL_MAX : (FFTSmoo_table[i+jM][0] / sqrt(var_after[i+jM][0]));
      }
    }
  }
  else {
    for (j=0; j<N2; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<N1; i++) {
	index_kMap      = i + jN1;
	totWeight       = gMap->map[index_kMap]->totWeight;
	SNR[index_kMap] = (totWeight < threshold) ? -DBL_MAX : (FFTSmoo_table[i+jM][0] * sigma_noise_inv);
      }
    }
  }
  return;
}

void kappaToSNR_DC(peak_param *pkPar, FFT_t *DCSmoo, signal_map *kMap, int DCScaleInd)
{
  //-- If local noise is used:
  //--   Compute the threshold of filling factor
  //--   Retrieve kappa from DCSmoo->before
  //--   Divide it by local std, which is stocked in DCSmoo->kernel (SNR = kappa / std_local)
  //--   Stock SNR in kMap->value1
  //--
  //-- If global noise is used:
  //--   Retrieve kappa from DCSmoo->before
  //--   Divide it by a filter-dependent global value (SNR = kappa / sigma_global)
  //--   Stock SNR in kMap->value1
  
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  int M  = DCSmoo->N;
  int bufferSize   = pkPar->bufferSize;
  double threshold = 0.0;
  fftw_complex *DCSmoo_kernel = DCSmoo->kernel; //-- Contain the variance and the filling factor
  
  int i, j, jM;
  
  //-- Compute the threshold of filling factor
  if (pkPar->doMask == 0) threshold = EPS_NUM; //-- Should not set to zero or negative values, otherwise empty pixels yield +infty in the variance calculation
  else {
    //-- Compute the filling threshold
    for (j=bufferSize; j<N2-bufferSize; j++) { 
      jM = j * M;
      for (i=bufferSize; i<N1-bufferSize; i++) threshold += DCSmoo_kernel[i+jM][1];
    }
    threshold /= (double)((N1 - 2*bufferSize) * (N2 - 2*bufferSize)); //-- Average filling factor
    threshold *= FILLING_THRESHOLD_RATIO;                             //-- Set threshold to the half of average
  }
  
  double sigma_noise_inv = 1.0 / pkPar->FFT_sigma_noise[DCScaleInd]; //-- Globle noise level
  fftw_complex *DCSmoo_before = DCSmoo->before;
  double *SNR = kMap->value1;
  double filling;
  int jN1, index_kMap, index_FFT;
  
  if (pkPar->doLocalNoise) {
    for (j=0; j<N2; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<N1; i++) {
	index_kMap      = i + jN1;
	index_FFT       = i + jM;
	filling         = DCSmoo_kernel[index_FFT][1];
	SNR[index_kMap] = (filling < threshold) ? -DBL_MAX : (DCSmoo_before[index_FFT][0] / sqrt(DCSmoo_kernel[index_FFT][0]));
      }
    }
  }
  else {
    for (j=0; j<N2; j++) {
      jN1 = j * N1;
      jM  = j * M;
      for (i=0; i<N1; i++) {
	index_kMap      = i + jN1;
	index_FFT       = i + jM;
	filling         = DCSmoo_kernel[index_FFT][1];
	SNR[index_kMap] = (filling < threshold) ? -DBL_MAX : (DCSmoo_before[index_FFT][0] * sigma_noise_inv);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to peak selection

int isPeak(double *SNR, int N1, int i, int j)
{
  double centerValue = SNR[i+j*N1];
  int m, n, jn_N1;
  
  for (n=-1; n<=1; n++) {
    jn_N1 = (j + n) * N1;
    for (m=-1; m<=1; m++) {
      if (m==0 && n==0) continue;
      if (SNR[(i+m)+jn_N1] >= centerValue) return 0;
    }
  }
  return 1;
}

int isPeak_float(float *SNR, int N1, int i, int j)
{
  float centerValue = SNR[i+j*N1];
  int m, n, jn_N1;
  for (n=-1; n<=1; n++) {
    jn_N1 = (j + n) * N1;
    for (m=-1; m<=1; m++) {
      if (m==0 && n==0) continue;
      if (SNR[(i+m)+jn_N1] >= centerValue) return 0;
    }
  }
  return 1;
}

int isPeakForTable(fftw_complex *table, int M, int i, int j)
{
  //-- table should have been turned into S/N map.
  
  double centerValue = table[i+j*M][0];
  int m, n, jn_M;
  for (n=-1; n<=1; n++) {
    jn_M = (j + n) * M;
    for (m=-1; m<=1; m++) {
      if (m==0 && n==0) continue;
      if (table[(i+m)+jn_M][0] >= centerValue) return 0;
    }
  }
  return 1;
}

void selectPeaks(peak_param *pkPar, signal_map *kMap, double_arr *peakList, error **err)
{
  //-- kMap should have been turned into S/N map.
  //-- Mask has been taken into account in kappaToSNR.
  
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  int count = 0;
  int bufferSize = (pkPar->field < 2) ? pkPar->bufferSize : 1;
  double *SNR    = kMap->value1;
  double *lfArr  = peakList->array;
  int i, j, jN1;
  
  for (j=bufferSize; j<N2-bufferSize; j++) {
    jN1 = j * N1;
    for (i=bufferSize; i<N1-bufferSize; i++) {
      if (isPeak(SNR, N1, i, j)) {
	lfArr[count] = SNR[i+jN1];
	count++;
      }
    }
  }
  
  peakList->length = count;
  return;
}

void cutSmallPeaks(double_arr *peakList, double nu_min)
{
  int count = 0;
  int i;
  for (i=0; i<peakList->length; i++) {
    if (peakList->array[i] >= nu_min) {
      peakList->array[count] = peakList->array[i];
      count++;
    }
  }
  peakList->length = count;
  return;
}

void outAsciiPeakField(FILE *file, peak_param *pkPar)
{
  fprintf(file, "# Buffer size = %d [pix], peak field = (%d, %d) [pix]\n", pkPar->bufferSize, pkPar->peakFieldResol[0], pkPar->peakFieldResol[1]);
  return;
}

void outAsciiPeakList(char name[], peak_param *pkPar, double_arr *peakList, int filterInd, error **err)
{
  if (pkPar->outPeakList == 0) return;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Peak list\n");
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiPeakField(file, pkPar);
  outAsciiNoiseInfo(file, pkPar);
  outAsciiFilterInfo(file, pkPar, 3, filterInd); //-- type = 3 (K+N map)
  fprintf(file, "#\n");
  fprintf(file, "# Number of peaks = %d\n", peakList->length); 
  fprintf(file, "#\n");
  fprintf(file, "#      SNR\n");
  fprintf(file, "#      [-]\n");
  
  fprintf(file, "#\n");
  
  int i;
  for (i=0; i<peakList->length; i++) fprintf(file, "  %8.5f\n", peakList->array[i]);
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

#ifdef __CAMELUS_USE_FITS__
void outFitsPeakField(FITS_t *fits, peak_param *pkPar)
{
  addLineSpread(fits);
  addKeyword(fits, TINT, "BUFF",     &pkPar->bufferSize,        "[pix] Border buffer size");
  addKeyword(fits, TINT, "PKRESOLX", &pkPar->peakFieldResol[0], "[pix] Peak field resolution x");
  addKeyword(fits, TINT, "PKRESOLY", &pkPar->peakFieldResol[1], "[pix] Peak field resolution y");
  return;
}
#endif

void outFitsPeakList(char name[], peak_param *pkPar, double_arr *peakList, int filterInd)
{
  if (pkPar->outPeakList == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  addColumn(fits, "SNR", TFLOAT, "-       ");
  
  float fBuff;
  int i;
  
  for (i=0; i<peakList->length; i++) {
    fBuff = (float)peakList->array[i]; writeTableColumn(fits, 0, 1, &fBuff);
    nextRow(fits);
  }
  
  addLineSpread(fits);
  addKeyword(fits, TINT,    "NBPEAKS", &peakList->length, "[-] Number of peaks");
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", &pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsPeakField(fits, pkPar);
  outFitsNoiseInfo(fits, pkPar);
  outFitsFilterInfo(fits, pkPar, 3, filterInd); //-- type = 3 (K+N map)
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

//----------------------------------------------------------------------
//-- Functions related to histogram

void setHist_nu(peak_param *pkPar, hist_t *hist)
{
  int i;
  for (i=0; i<hist->length; i++) {
    hist->x_lower[i] = pkPar->bin_nu[i];
    hist->n[i]       = 0;
  }
  hist->x_max = pkPar->bin_nu[hist->length];
  hist->dx    = -1.0;
  hist->n_tot = 0;
  return;
}

void makeHist(double_arr *peakList, hist_t *hist, int verbose)
{
  reset_hist_t(hist);
  int i;
  for (i=0; i<peakList->length; i++) push_hist_t(hist, peakList->array[i], verbose);
  return;
}

void outAsciiHistInfo(FILE *file, peak_param *pkPar)
{
  fprintf(file, "# doLocalNoise = %d, N_nu = %d\n", pkPar->doLocalNoise, pkPar->N_nu);
  return;
}

void outAsciiHist(char name[], peak_param *pkPar, hist_t *hist, int filterInd, error **err)
{
  if (pkPar->outHist == 0) return;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Peak histogram\n");
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiPeakField(file, pkPar);
  outAsciiNoiseInfo(file, pkPar);
  outAsciiFilterInfo(file, pkPar, 3, filterInd); //-- type = 3 (K+N map)
  outAsciiHistInfo(file, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Total counts = %d\n", hist->n_tot);
  fprintf(file, "#\n");
  fprintf(file, "# x_lower  x_upper       N\n");
  
  int i;
  for (i=0; i<hist->length-1; i++) fprintf(file, "  %7.4f  %7.4f  %6d\n", hist->x_lower[i], hist->x_lower[i+1], hist->n[i]);
  fprintf(file, "  %7.4f  %7.1e  %6d\n", hist->x_lower[hist->length-1], hist->x_max, hist->n[hist->length-1]);
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

#ifdef __CAMELUS_USE_FITS__
void outFitsHistInfo(FITS_t *fits, peak_param *pkPar)
{
  addLineSpread(fits);
  addKeyword(fits, TINT, "LOCNOISE", &pkPar->doLocalNoise, "[-] 0 = uniform global noise level, 1 = local noise");
  addKeyword(fits, TINT, "NNU",      &pkPar->N_nu,         "[-] Number of S/N bins");
  return;
}
#endif

void outFitsHist(char name[], peak_param *pkPar, hist_t *hist, int filterInd)
{
  if (pkPar->outHist == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  addColumn(fits, "XLOWER", TFLOAT, "-       ");
  addColumn(fits, "XUPPER", TFLOAT, "-       ");
  addColumn(fits, "N",      TINT,   "-       ");
  
  float fBuff;
  int i;
  
  for (i=0; i<hist->length-1; i++) {
    fBuff = (float)hist->x_lower[i];              writeTableColumn(fits, 0, 1, &fBuff);
    fBuff = (float)hist->x_lower[i+1];            writeTableColumn(fits, 1, 1, &fBuff);
                                                  writeTableColumn(fits, 2, 1, &hist->n[i]);
    nextRow(fits);
  }
    fBuff = (float)hist->x_lower[hist->length-1]; writeTableColumn(fits, 0, 1, &fBuff);
    fBuff = (float)hist->x_max;                   writeTableColumn(fits, 1, 1, &fBuff);
                                                  writeTableColumn(fits, 2, 1, &hist->n[hist->length-1]);
    nextRow(fits);
  
  addLineSpread(fits);
  addKeyword(fits, TINT,    "NNU",  &pkPar->N_nu,         "[-] Number of S/N bins");
  addKeyword(fits, TDOUBLE, "DNU",     &hist->dx,         "[-] Binwidth");
  addKeyword(fits, TINT,    "NTOT",    &hist->n_tot,      "[-] Total counts");
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsPeakField(fits, pkPar);
  outFitsNoiseInfo(fits, pkPar);
  outFitsFilterInfo(fits, pkPar, 3, filterInd); //-- type = 3 (K+N map)
  outFitsHistInfo(fits, pkPar);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

//----------------------------------------------------------------------

