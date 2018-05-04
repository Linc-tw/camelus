

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

// MKDEBUG: New name in Linc-tw, used to be computeLocalVariance. Used sigma_half_sq instead of 1/totweight.
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

// MKDEBUG: New name in Linc-tw, used to be computeLocalVariance_arr. But modified,
// no longer using peak->sigma_half.!
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

// MKDEBUG New in Linc-tw: FFTScaleInd indicates scale for which global noise is estimated.
// Before only local noise was implemented here. Also new: totweight instead of size.
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

// MKDEBUG See kappaToSNR_TTF for changes wrt Linc-tw
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

// MKDEBUG map_t -> signal_map according to Lin-tw
void computePeaks2(char name[], peak_param *peak, signal_map *kMap, double_arr *peakList, error **err)
{
  //-- kMap should have been turned into S/N map.
  //-- Mask has been taken into account in kappaToSNR.


  int N1    = kMap->N1;
  int N2    = kMap->N2;
  int count = 0;
  int bufferSize = peak->bufferSize;
  double *SNR    = kMap->value1;
  double *fArr   = peakList->array;
  int i, j, jN1;
  double pos[2];

  for (j=bufferSize; j<N2-bufferSize; j++) {
    jN1 = j * N1;
    for (i=bufferSize; i<N1-bufferSize; i++) {
      if (isPeak(SNR, N1, i, j)) {
		//fArr[count] = SNR[i+jN1];
		count++;
      }
    }
  }

   FILE *file = fopen(name, "w");
  fprintf(file, "# Peak list\n");
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  if (peak->FFT_nbFilters) fprintf(file, "# FFT filter = %s, FFT scale = %g [arcmin] = %g [pix]\n", STR_FILTER_T(peak->FFT_filter[0]), peak->FFT_scale[0], peak->FFT_scaleInPix[0]);
  if (peak->DC_nbFilters)  fprintf(file, "# DC filter = %s, DC scale = %g [arcmin] = %g [pix]\n", STR_FILTER_T(peak->DC_filter[0]),  peak->DC_scale[0],  peak->DC_scaleInPix[0]);
  fprintf(file, "# sigma_eps = %g, sigma_pix = %g, sigma_noise = %g\n", peak->sigma_eps, peak->sigma_pix, peak->FFT_sigma_noise[0]);
  fprintf(file, "# Buffer size = %d [pix]\n", peak->bufferSize);
  fprintf(file, "#\n");
  fprintf(file, "# Number of pixels = %d\n", peakList->length); 
  fprintf(file, "#\n");
  fprintf(file, "#      SNR \t theta_x  \t theta_y  \n");
  fprintf(file, "#      [-]  \t [arcmin]  \t [arcmin]  \n");

  pos[0]=1;
  pos[1]=1;
  bufferSize=0;

  	for (j=bufferSize; j<N2-bufferSize; j++) {
    		jN1 = j * N1;
    		for (i=bufferSize; i<N1-bufferSize; i++) {
      		if (isPeak(SNR, N1, i, j)) {
	  		//getPixPos_map_t(kMap, pos, i, j); // MKDEBUG: this call had no effect. Function removed in Linc-tw
			getPixPos(pos, peak->theta_pix, i, j);
  	  		fprintf(file, " %8.5f     %5.3f      %5.3f \n",SNR[i+jN1], pos[0], pos[1]);
      }
    }
  }  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

// MKDEBUG Remove selectPeaks_mrlens, as in Lin-tw

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

// New functions for TablesRondes
     
void doPeakList_withInputs(char fileName[], char fileName2[],char end[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
  char fpeakList[STRING_LENGTH_MAX];
  char fpeakListPos[STRING_LENGTH_MAX];
  char fpeakHist[STRING_LENGTH_MAX];

  sprintf(fpeakList, "peakList_%s",end);
  sprintf(fpeakListPos, "peakListPos_%s",end);
  sprintf(fpeakHist, "peakhist_%s",end);

  // MKDEBUG: New structure in Linc-tw, replaced lines below.
  pipeline_t *pipe = initialize_pipeline_t(cmhm, peak, err); forwardError(*err, __LINE__,);

  /*
  halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp   = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                       forwardError(*err, __LINE__,);
  gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask   = initializeMask(peak, err);                                                 forwardError(*err, __LINE__,);
  FFT_arr *FFTSmoother = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernel(peak, FFTSmoother);
  FFT_arr *DCSmoother  = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  map_t *kMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_arr *variance    = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernelForVariance(peak, variance);
  double_arr *peakList = initialize_double_arr(length);
  hist_t *nuHist       = initialize_hist_t(peak->N_nu);
  setHist_nu(peak, nuHist);
  */

  if ((fileName == NULL)||(fileName2 == NULL)) {
    //-- no input files
    printf("Problem input files missing \n");
  }
  else {
    printf(" Input file for halo : \"%s\" \n", fileName);
    printf(" Input file for galaxies : \"%s\" \n", fileName2);
    read_halo_map(fileName, cmhm, peak, pipe->hMap, err);                                forwardError(*err, __LINE__,);
    read_gal_map2(fileName2, cmhm, peak, pipe->gMap, err);                                forwardError(*err, __LINE__,);
  }
  makeMapAndOutputAll2(fileName, fileName2, cmhm, peak, pipe->gMap, pipe->FFTSmoother, pipe->DCSmoother, pipe->kMap, err); forwardError(*err, __LINE__,);

  // MKDEBUG  used to be computeLocalVariance_arr
  makeLocalVariance(peak, pipe->gMap, pipe->variance);
   
  // MKDEBUG: Replaced with lines from multiscale.c:mapToMultiscale_FFT/DC
  //if (peak->DC_nbFilters) kappaToSNR_DC(peak, pipe->gMap, pipe->DCSmoother->array[0], pipe->kMap);
  //else                         kappaToSNR_FFT(peak, pipe->gMap, pipe->FFTSmoother->array[0], pipe->kMap, pipe->variance->array[0]);
  int j;
  for (j=0; j<peak->FFT_nbFilters; j++) {
    kappaToSNR_FFT(peak, pipe->gMap, pipe->FFTSmoother->array[j], pipe->kMap, pipe->variance->array[j], j);
  }
  for (j=0; j<peak->DC_nbFilters; j++) {
    kappaToSNR_DC(peak, pipe->DCSmoother->array[j], pipe->kMap, j);
  }

  selectPeaks(peak, pipe->kMap, pipe->peakList, err);   forwardError(*err, __LINE__,);
 
  // MKDEBUG replaced. The following will only write peaks from 0th filter!
  //outputPeakList(fpeakList, peak, pipe->peakList);
  outAsciiPeakList(fpeakList, peak, pipe->peakList, 0, err);
  forwardError(*err, __LINE__,);

  computePeaks2(fpeakListPos, peak, pipe->kMap, pipe->peakList, err);
  int silent = 1;
  makeHist(pipe->peakList, pipe->nuHist, silent);

  // MKDEBUG replaced. Output of 0th filter.
  //outputHist(fpeakHist, pipe->nuHist);
  outAsciiHist(fpeakHist, peak, pipe->nuHist, 0, err);
  forwardError(*err, __LINE__,);

  free_pipeline_t(pipe);

  printf("------------------------------------------------------------------------\n");
  return;
}

void doProduce_Catalog(char HaloFileName[],char GalFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);

  // MKDEBUG: New structure in Linc-tw, replaced lines below.
  pipeline_t *pipe = initialize_pipeline_t(cmhm, peak, err); forwardError(*err, __LINE__,);

  /*  
  halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp   = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                       forwardError(*err, __LINE__,);
  gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask   = initializeMask(peak, err);                                                 forwardError(*err, __LINE__,);
  FFT_arr *FFTSmoother = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernel(peak, FFTSmoother);
  FFT_arr *DCSmoother  = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  map_t *kMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_arr *variance    = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernelForVariance(peak, variance);
  double_arr *peakList = initialize_double_arr(length);
  hist_t *nuHist       = initialize_hist_t(peak->N_nu);
  setHist_nu(peak, nuHist);
  */

  //-- Carry out fast simulation
  sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  // MKDEBUG: New in Linc-tw, set doVolume=1
  setMassSamplers(cmhm, peak, sampArr, 1, err);      forwardError(*err, __LINE__,);
  makeFastSimul(cmhm, peak, sampArr, pipe->hMap, err);    forwardError(*err, __LINE__,);

  // MKDEBUG Replaced, new name in Linc-tw
  //outputFastSimul(HaloFileName, cmhm, peak, hMap);
  outAsciiHaloCat(HaloFileName, cmhm, peak, pipe->hMap, err);
  free_sampler_arr(sampArr);

  cleanOrMakeOrResample(cmhm, peak, pipe->gSamp, pipe->gMap, pipe->mask, err);            forwardError(*err, __LINE__,);
  lensingCatalogueAndOutputAll2(GalFileName, cmhm, peak, pipe->hMap, pipe->gMap, pipe->k1Inter, err);

  free_pipeline_t(pipe);

  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------
void doProduce_Catalog_N(int N,char HaloFileName[],char GalFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
    
  char HaloFileName2[STRING_LENGTH_MAX];
  char GalFileName2[STRING_LENGTH_MAX];
  int i;

  for (i=0; i<N; i++) {
    sprintf(HaloFileName2, "%s_%3.3d",HaloFileName, i+1);
    sprintf(GalFileName2, "%s_%3.3d",GalFileName, i+1);
    printf(" >>> Compute halo %s \n", HaloFileName2);
    printf(" >>> Compute gal %s \n", GalFileName2);

    // MKDEBUG: New structure in Linc-tw, replaced lines below.
    pipeline_t *pipe = initialize_pipeline_t(cmhm, peak, err); forwardError(*err, __LINE__,);

    /*
    halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
    sampler_t *galSamp   = initialize_sampler_t(peak->N_z_gal);
    setGalaxySampler(cmhm, peak, galSamp, err);                                                       forwardError(*err, __LINE__,);
    gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
    short_mat *CCDMask   = initializeMask(peak, err);                                                 forwardError(*err, __LINE__,);
    FFT_arr *FFTSmoother = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
    if (peak->FFT_nbFilters) makeKernel(peak, FFTSmoother);
    FFT_arr *DCSmoother  = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
    map_t *kMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
    FFT_arr *variance    = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
    if (peak->FFT_nbFilters) makeKernelForVariance(peak, variance);
    double_arr *peakList = initialize_double_arr(length);
    hist_t *nuHist       = initialize_hist_t(peak->N_nu);
    setHist_nu(peak, nuHist);
    */

    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
    // MKDEBUG: New in Linc-tw, set doVolume=1
    setMassSamplers(cmhm, peak, sampArr, 1, err);      forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, pipe->hMap, err);    forwardError(*err, __LINE__,);

    // MKDEBUG Replaced, new name in Linc-tw
    //outputFastSimul(HaloFileName2, cmhm, peak, hMap);
    outAsciiHaloCat(HaloFileName, cmhm, peak, pipe->hMap, err);
    free_sampler_arr(sampArr);

    cleanOrMakeOrResample(cmhm, peak, pipe->gSamp, pipe->gMap, pipe->mask, err);
    forwardError(*err, __LINE__,);
    lensingCatalogueAndOutputAll2(GalFileName2, cmhm, peak, pipe->hMap, pipe->gMap, pipe->k1Inter, err);

    free_pipeline_t(pipe);

    printf("------------------------------------------------------------------------\n");

  }
  return;
}


void doPeakList_withInputs_N(int N,char fileName[], char fileName2[],char end[],cosmo_hm *cmhm, peak_param *peak, error **err)
{
    
  char HaloFileName2[STRING_LENGTH_MAX];
  char GalFileName2[STRING_LENGTH_MAX];
  char PeakHistfich2[STRING_LENGTH_MAX];
  char PeakListfich2[STRING_LENGTH_MAX];
  int i;

  for (i=0; i<N; i++) {

    sprintf(HaloFileName2, "%s_%3.3d",fileName, 1);
    sprintf(GalFileName2, "%s_%3.3d",fileName2, i+1);
    sprintf(PeakListfich2, "PeakList_%s_%3.3d",end, i+1);
    sprintf(PeakHistfich2, "PeakHist_%s_%3.3d",end, i+1);


    printf(" >>> Read halo %s \n", HaloFileName2);
    printf(" >>> Read gal %s \n", GalFileName2);
    printf(" >>> Write histo peak %s \n", PeakListfich2);
    printf(" >>> Write list peak %s \n",PeakHistfich2 );


    int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);

    // MKDEBUG: New structure in Linc-tw, replaced lines below.
    pipeline_t *pipe = initialize_pipeline_t(cmhm, peak, err); forwardError(*err, __LINE__,);

    /*
       halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
       sampler_t *galSamp   = initialize_sampler_t(peak->N_z_gal);
       setGalaxySampler(cmhm, peak, galSamp, err);                                                       forwardError(*err, __LINE__,);
       gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
       short_mat *CCDMask   = initializeMask(peak, err);                                                 forwardError(*err, __LINE__,);
       FFT_arr *FFTSmoother = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
       if (peak->FFT_nbFilters) makeKernel(peak, FFTSmoother);
       FFT_arr *DCSmoother  = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
       map_t *kMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
       FFT_arr *variance    = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
       if (peak->FFT_nbFilters) makeKernelForVariance(peak, variance);
       double_arr *peakList = initialize_double_arr(length);
       hist_t *nuHist       = initialize_hist_t(peak->N_nu);
       setHist_nu(peak, nuHist);
       */

  // MKDEBUG: Check new Linc-tw funtion readCatOrMakeSimulAndOutput

    if ((fileName == NULL)||(fileName2 == NULL)) {
      //-- no input files
      printf("Problem input files missing \n");
    } else {
      printf(" Input file for halo : \"%s\" \n", HaloFileName2);
      printf(" Input file for galaxies : \"%s\" \n", GalFileName2);
      read_halo_map(HaloFileName2, cmhm, peak, pipe->hMap, err);                                forwardError(*err, __LINE__,);
      read_gal_map(GalFileName2, cmhm, peak, pipe->gMap, err);                                forwardError(*err, __LINE__,);
    }
    makeMapAndOutputAll2(HaloFileName2, GalFileName2, cmhm, peak, pipe->gMap, pipe->FFTSmoother, pipe->DCSmoother, pipe->kMap, err); forwardError(*err, __LINE__,);

    // MKDEBUG used to be computeLocalVariance_arr
    makeLocalVariance(peak, pipe->gMap, pipe->variance);

    // MKDEBUG: Replaced with lines from multiscale.c:mapToMultiscale_FFT/DC
    //if (peak->DC_nbFilters) kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
    //else                         kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);
    int j;
    for (j=0; j<peak->FFT_nbFilters; j++) {
      kappaToSNR_FFT(peak, pipe->gMap, pipe->FFTSmoother->array[j], pipe->kMap, pipe->variance->array[j], j);
    }
    for (j=0; j<peak->DC_nbFilters; j++) {
      kappaToSNR_DC(peak, pipe->DCSmoother->array[j], pipe->kMap, j);
    }

    selectPeaks(peak, pipe->kMap, pipe->peakList, err);   forwardError(*err, __LINE__,);

    // MKDEBUG replaced. The following will only write peaks from 0th filter!
    //outputPeakList(PeakListfich2, peak, peakList);
    outAsciiPeakList(PeakListfich2, peak, pipe->peakList, 0, err);

      // computePeaks2("peakListPos", peak, pipe->kMap, pipe->peakList, err);

    int silent = 1;
    makeHist(pipe->peakList, pipe->nuHist, silent); 

    // MKDEBUG replaced. Output of 0th filter.
    //outputHist(PeakHistfich2, nuHist);
    outAsciiHist(PeakHistfich2, peak, pipe->nuHist, 0, err);
    forwardError(*err, __LINE__,);
    //computePeaks2("TEST_TABLE_PEAK",peak,kMap,peakList,err);

    free_pipeline_t(pipe);

    printf("------------------------------------------------------------------------\n");
  }
  return;
}

void doProduce_Catalog_DM_HOD(int N,char CmhmName[],char HaloFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  char HaloFileName2[STRING_LENGTH_MAX];
  int i;


  printf("-----------------------------  HOD Ngal  -------------------------------\n");
  for (i=0; i<N; i++) {
    halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);
    forwardError(*err,__LINE__,);
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
    setMassSamplers(cmhm, peak, sampArr, 1, err); 
    forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err);
    forwardError(*err, __LINE__,);
    sprintf(HaloFileName2, "%s_%3.3d", HaloFileName, i+1);
    // MKDEBUG Replaced, new name in Linc-tw
    //outputFastSimul_HOD(CmhmName, HaloFileName2, cmhm, peak, hMap);
    outAsciiHaloCat(HaloFileName, cmhm, peak, hMap, err);
    forwardError(*err, __LINE__,);
    free_sampler_arr(sampArr);
    free_halo_map(hMap);
  }
  return;
}


void doPeakList_withInputs_hod(char fileNameHal[], char fileNameGal[], char end[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
  char fpeakList[STRING_LENGTH_MAX];
  char fpeakListPos[STRING_LENGTH_MAX];
  char fpeakHist[STRING_LENGTH_MAX];

  sprintf(fpeakList, "peakList_%s",end);
  sprintf(fpeakListPos, "peakListPos_%s",end);
  sprintf(fpeakHist, "peakhist_%s",end);

  // MKDEBUG: New structure in Linc-tw, replaced lines below.
  pipeline_t *pipe = initialize_pipeline_t(cmhm, peak, err); forwardError(*err, __LINE__,);

  /*
  halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp   = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                       forwardError(*err, __LINE__,);
  gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask   = initializeMask(peak, err);                                                 forwardError(*err, __LINE__,);
  FFT_arr *FFTSmoother = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernel(peak, FFTSmoother);
  FFT_arr *DCSmoother  = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  map_t *kMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);

  FFT_arr *variance    = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernelForVariance(peak, variance);
  double_arr *peakList = initialize_double_arr(length);
  hist_t *nuHist       = initialize_hist_t(peak->N_nu);
  setHist_nu(peak, nuHist);
  */

  if ((fileNameHal == NULL)||(fileNameGal == NULL)) {
    //-- no input files
    printf("Problem input files missing \n");
  }
  else {
    printf(" Input file for halo : \"%s\" \n", fileNameHal);
    printf(" Input file for galaxies : \"%s\" \n", fileNameGal);
    read_halo_map(fileNameHal, cmhm, peak, pipe->hMap, err); 
    forwardError(*err, __LINE__,);
    printf("halo read \n");
    read_gal_map2(fileNameGal, cmhm, peak, pipe->gMap, err);
    forwardError(*err, __LINE__,);
    printf("gal read \n");
  }


  makeMapAndOutputAll2(fileNameHal, fileNameGal, cmhm, peak, pipe->gMap, pipe->FFTSmoother, pipe->DCSmoother, pipe->kMap, err);
  forwardError(*err, __LINE__,);

  // MKDEBUG used to be computeLocalVariance_arr
  makeLocalVariance(peak, pipe->gMap, pipe->variance);

  // MKDEBUG: Replaced with lines from multiscale.c:mapToMultiscale_FFT/DC
  //if (peak->DC_nbFilters) kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
  //else                         kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);
  int j;
  for (j=0; j<peak->FFT_nbFilters; j++) {
    kappaToSNR_FFT(peak, pipe->gMap, pipe->FFTSmoother->array[j], pipe->kMap, pipe->variance->array[j], j);
  }
  for (j=0; j<peak->DC_nbFilters; j++) {
    kappaToSNR_DC(peak, pipe->DCSmoother->array[j], pipe->kMap, j);
  }

  // MKDEBUG replaced after Linc-tw. Needs to be checked.
  //lensingCatalogueAndOutputAll(cmhm, peak, pipe->hMap, pipe->gMap, err);                 forwardError(*err, __LINE__,);
  sprintf(peak->prefix, "galCat");
  lensingCatalogueAndOutput(cmhm, peak, pipe->hMap, pipe->gMap, pipe->k1Inter, err); forwardError(*err, __LINE__,);

  // MKDEBUG replaced to Linc-tw. Question: is noise now added, whereas before it was not?
  // Also: What about gMap2?
  //makeMapAndOutputAll(cmhm, peak, pipe->gMap, pipe->FFTSmoother, pipe->DCSmoother, pipe->kMap, err); forwardError(*err, __LINE__,);
  gal_map *gMap2 = NULL;
  makeMapsAndOutput(cmhm, peak, pipe->gMap, gMap2, pipe->FFTSmoother, pipe->DCSmoother, pipe->kMap, err);   forwardError(*err, __LINE__,);
  // MKDEBUG used to be computeLocalVariance_arr
  makeLocalVariance(peak, pipe->gMap, pipe->variance);

   
  selectPeaks(peak, pipe->kMap, pipe->peakList, err);   forwardError(*err, __LINE__,);


  // MKDEBUG replaced. The following will only write peaks from 0th filter!
  //outputPeakList(fpeakList, peak, peakList);
  outAsciiPeakList(fpeakList, peak, pipe->peakList, 0, err);

  computePeaks2(fpeakListPos, peak, pipe->kMap, pipe->peakList, err);
  int silent = 1;
  makeHist(pipe->peakList, pipe->nuHist, silent);
  //outputHist(fpeakHist, nuHist);
  outAsciiHist(fpeakHist, peak, pipe->nuHist, 0, err);
  forwardError(*err, __LINE__,);


  //computePeaks2("TEST_TABLE_PEAK",peak,kMap,peakList,err);

  free_pipeline_t(pipe);

  printf("------------------------------------------------------------------------\n");
  return;
}

void doProduce_Catalog_DM_galaxies(int N, char CmhmName[], char HaloFileName[], char GalaxyFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  char HaloFileName2[STRING_LENGTH_MAX], GalaxyFileName2[STRING_LENGTH_MAX];
  int i;

  
  printf("-----------------------------  HOD galaxy Ngal  -------------------------------\n");
  for (i=0; i<N; i++) {
    halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);
    forwardError(*err,__LINE__,);
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
    setMassSamplers(cmhm, peak, sampArr, 1, err); 
    forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err);
    forwardError(*err, __LINE__,);
    sprintf(HaloFileName2, "%s_%3.3d",HaloFileName, i+1);
    sprintf(GalaxyFileName2, "%s_%3.3d",GalaxyFileName, i+1);
    //printf("test \n");
    // MKDEBUG Replaced, new name in Linc-tw
    //outputFastSimul_galaxies(CmhmName, HaloFileName2, GalaxyFileName2, cmhm, peak, hMap);
    outAsciiHaloCat(HaloFileName, cmhm, peak, hMap, err);
    forwardError(*err, __LINE__,);
    //printf("test2 \n");
    free_sampler_arr(sampArr);
    free_halo_map(hMap);
  }
  return;
}
