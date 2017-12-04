

  /***************************************
   **  peakSelection.c			**
   **  Chieh-An Lin, François Lanusse	**
   **  Version 2016.03.20		**
   ***************************************/


#include "peakSelection.h"


//----------------------------------------------------------------------
//-- Functions related to local variance

void fillGaussianKernelForVariance(fftw_complex *kernel, int length, int M, double scaleInPix)
{
  //-- If K = sum(w_i kappa_i), then variance = sum(w_i^2 sigma_i^2) / sum(sigma_i) ^2.
  
  //-- Initialization
  reset_fftw_complex(kernel, length);
  
  double scaleInPix_invSq = 1.0 / (SQ(scaleInPix));         //-- [pix^2]
  double relay       = CUTOFF_FACTOR_GAUSSIAN * scaleInPix; //-- CUTOFF_FACTOR_GAUSSIAN set to 2.2 in peakParameters.h, equivalent to a Gaussian trunked at 3-sigma
  int cutInPix       = (int)ceil(relay);
  double cutInPix_sq = SQ(relay);
  double sum         = 0.0;
  double value;
  int i, j, mod_j_M, r_sq;
  
  for (j=-cutInPix; j<=cutInPix; j++) {
    mod_j_M = POS_MOD(M, j) * M;
    for (i=-cutInPix; i<=cutInPix; i++) {
      r_sq = i*i + j*j;
      if (r_sq > cutInPix_sq) continue;
      value = exp(-r_sq * scaleInPix_invSq);
      sum += value;
      kernel[POS_MOD(M, i) + mod_j_M][0] = SQ(value);
    }
  }
  
  //-- Normalization
  rescaleReal_fftw_complex(kernel, length, 1.0 / (sum*sum));
  return;
}

void fillStarletKernelForVariance(fftw_complex *kernel, int length, int M, double scaleInPix)
{
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
      sum  += fabs(value);
      kernel[POS_MOD(M, i) + mod_j_M][0] = SQ(value);
    }
  }
  
  //-- Normalization
  rescaleReal_fftw_complex(kernel, length, 1.0 / (sum*sum));
  return;
}

void makeKernelForVariance(peak_param *peak, FFT_arr *variance)
{
  FFT_t *var;
  int i;
  for (i=0; i<peak->FFT_nbFilters; i++) {
    var = variance->array[i];
    if (peak->FFT_filter[i] == gauss)     fillGaussianKernelForVariance(var->kernel, var->length, var->N, peak->FFT_scaleInPix[i]);
    else if (peak->FFT_filter[i] == star) fillStarletKernelForVariance(var->kernel, var->length, var->N, peak->FFT_scaleInPix[i]);
    fftw_execute(var->kernel_f); //-- To Fourier space
  }
  return;
}

void computeLocalVariance(gal_map *gMap, FFT_t *var, double sigma_half_sq)
{
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int M  = var->N;
  double threshold         = gMap->fillingThreshold;
  fftw_complex *var_before = var->before;
  
  int i, j, jN1, jM, index_FFT, size;
  
  //-- Fill sigma_eps^2 / 2N in var->before
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
	size = gMap->map[i+jN1]->size;
	if (size == 0 || (double)size < threshold) var_before[index_FFT][0] = 0.0;
	else                                       var_before[index_FFT][0] = sigma_half_sq / (double)size;
	var_before[index_FFT][1] = 0.0;
      }
    }
  }
  
  execute_FFT_t(var);
  return;
}

void computeLocalVariance_arr(peak_param *peak, gal_map *gMap, FFT_arr *variance)
{
  if (peak->FFT_nbFilters == 0) return;
  
  double sigma_half_sq = SQ(peak->sigma_half);
  FFT_t *var           = variance->array[0];
  
  computeLocalVariance(gMap, var, sigma_half_sq);
  
  fftw_complex *var_before = var->before;
  int FFTLength = var->length;
  int i;
  
  for (i=1; i<variance->length; i++) {
    var = variance->array[i];
    multiplication_fftw_complex(var_before, var->kernel, var->after, FFTLength); //-- Multiplication
    fftw_execute(var->after_b);                                                  //-- Go to direct space
    rescaleReal_fftw_complex(var->after, FFTLength, peak->FFTNormFactor);        //-- Rescale, only real part is interesting.
  }
  return;
}

void kappaToSNR_FFT(peak_param *peak, gal_map *gMap, FFT_t *FFTSmoo, map_t *kMap, FFT_t *var)
{
  //-- Retrieve kappa from smoo_table, divide it by sigma_noise (becoming SNR),
  //-- and stock SNR in kMap->value1.
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int M  = var->N;
  double threshold = gMap->fillingThreshold;
  double *SNR      = kMap->value1;
  fftw_complex *FFTSmoo_table = FFTSmoo->after;
  fftw_complex *var_after     = var->after; //-- Contain the variance
  
  int i, j, jN1, jM, index_kMap, size;
  
  for (j=0; j<N2; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<N1; i++) {
      index_kMap = i + jN1;
      size = gMap->map[index_kMap]->size;
      SNR[index_kMap] = ((double)size < threshold) ? -DBL_MAX : (FFTSmoo_table[i+jM][0] / sqrt(var_after[i+jM][0]));
    }
  }
  return;
}

void kappaToSNR_DC(peak_param *peak, gal_map *gMap, FFT_t *DCSmoo, map_t *kMap)
{
  //-- Retrieve kappa from smoo_table, divide it by sigma_noise (becoming SNR),
  //-- and stock SNR in kMap->value1.
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  int M  = DCSmoo->N;
  int bufferSize   = peak->bufferSize;
  double threshold = 0.0;
  fftw_complex *DCSmoo_kernel = DCSmoo->kernel; //-- Contain the variance and the filling factor
  
  int i, j, jM;
  
  //-- Compute the filling threshold
  for (j=bufferSize; j<N2-bufferSize; j++) { 
    jM = j * M;
    for (i=bufferSize; i<N1-bufferSize; i++) threshold += DCSmoo_kernel[i+jM][1];
  }
  threshold /= (double)((N1 - 2*bufferSize) * (N2 - 2*bufferSize)); //-- Average filling factor
  threshold *= FILLING_THRESHOLD_RATIO;                             //-- Set threshold to the half of average
  
  fftw_complex *DCSmoo_before = DCSmoo->before;
  double *SNR = kMap->value1;
  double filling;
  int jN1, index_kMap, index_FFT;
  
  for (j=0; j<N2; j++) {
    jN1 = j * N1;
    jM  = j * M;
    for (i=0; i<N1; i++) {
      index_kMap = i + jN1;
      index_FFT  = i + jM;
      filling = DCSmoo_kernel[index_FFT][1];
      SNR[index_kMap] = (filling < threshold) ? -DBL_MAX : (DCSmoo_before[index_FFT][0] / sqrt(DCSmoo_kernel[index_FFT][0]));
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

void selectPeaks(peak_param *peak, map_t *kMap, double_arr *peakList, error **err)
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
  
  for (j=bufferSize; j<N2-bufferSize; j++) {
    jN1 = j * N1;
    for (i=bufferSize; i<N1-bufferSize; i++) {
      if (isPeak(SNR, N1, i, j)) {
	fArr[count] = SNR[i+jN1];
	count++;
      }
    }
  }
  
  peakList->length = count;
  return;
}

void computePeaks2(char name[],peak_param *peak, map_t *kMap, double_arr *peakList, error **err)
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
  double limit[4];

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
  if (peak->doNonlinear)   fprintf(file, "# Filter = %s, number of scales = %d [-], FDR = %g [-]\n", STR_FILTER_T(mrlens), peak->MRLens_nbScales, peak->MRLens_FDR);
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
	  		getPixPos_map_t(kMap, pos, i, j);
			getPixPos(pos, limit,peak->theta_pix,i,j);
  	  		fprintf(file, " %8.5f     %5.3f      %5.3f \n",SNR[i+jN1], pos[0], pos[1]);
      }
    }
  }  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}


void selectPeaks_mrlens(char name[], peak_param *peak, gal_map *gMap, double_arr *peakList)
{
#ifndef __releaseMenu__
  //-- This fuction selects peak from an image stocked in a float*.
  //-- Peaks are kappa peaks not S/N peaks.
  
  //-- Read nonlinear map
  FITS_t *fits = initializeImageReader_FITS_t(name);
  float *kappa = read2DImage(fits);
  free_FITS_t(fits);
  
  double threshold = gMap->fillingThreshold;
  int i, size;
  
  //-- Masking 
  for (i=0; i<gMap->length; i++) {
    size = gMap->map[i]->size;
    if ((double)size < threshold) kappa[i] = -DBL_MAX;
  }
  
  int N1    = peak->resol[0];
  int N2    = peak->resol[1];
  int count = 0;
  int bufferSize = peak->bufferSize;
  double *fArr = peakList->array;
  int j, jN1;
  
  for (j=bufferSize; j<N2-bufferSize; j++) {
    jN1 = j * N1;
    for (i=bufferSize; i<N1-bufferSize; i++) {
      if (isPeak_float(kappa, N1, i, j)) {
	fArr[count] = (double)kappa[i+jN1];
	count++;
      }
    }
  }
  
  peakList->length = count;
  free(kappa);
#endif
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

void outputPeakList(char name[], peak_param *peak, double_arr *peakList)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Peak list\n");
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  if (peak->FFT_nbFilters) fprintf(file, "# FFT filter = %s, FFT scale = %g [arcmin] = %g [pix]\n", STR_FILTER_T(peak->FFT_filter[0]), peak->FFT_scale[0], peak->FFT_scaleInPix[0]);
  if (peak->DC_nbFilters)  fprintf(file, "# DC filter = %s, DC scale = %g [arcmin] = %g [pix]\n", STR_FILTER_T(peak->DC_filter[0]),  peak->DC_scale[0],  peak->DC_scaleInPix[0]);
  if (peak->doNonlinear)   fprintf(file, "# Filter = %s, number of scales = %d [-], FDR = %g [-]\n", STR_FILTER_T(mrlens), peak->MRLens_nbScales, peak->MRLens_FDR);
  fprintf(file, "# sigma_eps = %g, sigma_pix = %g, sigma_noise = %g\n", peak->sigma_eps, peak->sigma_pix, peak->FFT_sigma_noise[0]);
  fprintf(file, "# Buffer size = %d [pix]\n", peak->bufferSize);
  fprintf(file, "#\n");
  fprintf(file, "# Number of pixels = %d\n", peakList->length); 
  fprintf(file, "#\n");
  fprintf(file, "#      SNR\n");
  fprintf(file, "#      [-]\n");
  
  int i;
  for (i=0; i<peakList->length; i++) fprintf(file, "  %8.5f\n", peakList->array[i]);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

void peakListFromMassFct(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask,
			 FFT_arr *FFTSmoother, FFT_arr *DCSmoother, map_t *kMap, FFT_arr *variance, double_arr *peakList, error **err)
{
  //-- Allow only one filter
  //-- Priority: nonlinear > direct convolution > FFT
  
  makeFastSimul(cmhm, peak, sampArr, hMap, err);                  forwardError(*err, __LINE__,);
  cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err); forwardError(*err, __LINE__,);
  lensingCatalogue(cmhm, peak, hMap, gMap, err);                  forwardError(*err, __LINE__,);
  makeMap(peak, gMap, FFTSmoother, DCSmoother, kMap, err);        forwardError(*err, __LINE__,);
  computeLocalVariance_arr(peak, gMap, variance);
  
  if (peak->doNonlinear) selectPeaks_mrlens("kappaMap_mrlens.fits", peak, gMap, peakList);
  else if (peak->DC_nbFilters) {
    kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
    selectPeaks(peak, kMap, peakList, err);                       forwardError(*err, __LINE__,);
  }
  else if (peak->FFT_nbFilters) {
    kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);
    selectPeaks(peak, kMap, peakList, err);                       forwardError(*err, __LINE__,);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to histogram

void setHist_nu(peak_param *peak, hist_t *hist)
{
  int i;
  for (i=0; i<hist->length; i++) {
    hist->x_lower[i] = peak->nu_bin[i];
    hist->n[i]       = 0;
  }
  hist->x_max = peak->nu_bin[hist->length];
  hist->dx    = -1.0;
  hist->n_tot = 0;
  return;
}

void setHist_kappa(peak_param *peak, hist_t *hist)
{
  int i;
  for (i=0; i<hist->length; i++) {
    hist->x_lower[i] = peak->kappa_bin[i];
    hist->n[i]       = 0;
  }
  hist->x_max = peak->kappa_bin[hist->length];
  hist->dx    = -1.0;
  hist->n_tot = 0;
  return;
}

void makeHist(double_arr *peakList, hist_t *hist, int silent)
{
  reset_hist_t(hist);
  int i;
  if (silent) {
    for (i=0; i<peakList->length; i++) silentPush_hist_t(hist, peakList->array[i]);
  }
  else {
    for (i=0; i<peakList->length; i++) push_hist_t(hist, peakList->array[i]);
  }
  return;
}

void outputHist(char name[], hist_t *hist)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Peak histogram\n");
  fprintf(file, "# Number of bins = %d\n", hist->length);
  fprintf(file, "# Bin width = %g\n", hist->dx);
  fprintf(file, "# Total counts = %d\n", hist->n_tot);
  fprintf(file, "#\n");
  fprintf(file, "# x_lower  x_upper       n\n");
  
  int i;
  double upper;
  for (i=0; i<hist->length; i++) {
    if (hist->dx > 0) {
      upper = hist->x_lower[i] + hist->dx;
    } else {
      if (i<hist->length-1) {
        upper = hist->x_lower[i+1];
      } else {
        upper = hist->x_max;
      }
    }
    fprintf(file, "  %7.4g  %7.4g  %6d\n", hist->x_lower[i], upper, hist->n[i]);
  }
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doPeakList(char fileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
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
  
  if (fileName == NULL) {
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
    setMassSamplers(cmhm, peak, sampArr, err);                               forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err);                           forwardError(*err, __LINE__,);
    outputFastSimul("haloCat", cmhm, peak, hMap);
    free_sampler_arr(sampArr);
  }
  else {
  	printf("\"%s\" read inpute \n", fileName);
    read_halo_map(fileName, cmhm, hMap, err);                                forwardError(*err, __LINE__,);
  }
  
  cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err);            forwardError(*err, __LINE__,);
  lensingCatalogueAndOutputAll(cmhm, peak, hMap, gMap, err);                 forwardError(*err, __LINE__,);
  makeMapAndOutputAll(cmhm, peak, gMap, FFTSmoother, DCSmoother, kMap, err); forwardError(*err, __LINE__,);
  computeLocalVariance_arr(peak, gMap, variance);
   
  if (peak->doNonlinear)       selectPeaks_mrlens("kappaMap_mrlens.fits", peak, gMap, peakList);
  else if (peak->DC_nbFilters) kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
  else                         kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);
  
  selectPeaks(peak, kMap, peakList, err);                                    forwardError(*err, __LINE__,);
  outputPeakList("peakList", peak, peakList);
  int silent = 1;
  makeHist(peakList, nuHist, silent);
  outputHist("peakHist", nuHist);
  
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(nuHist);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doPeakList_repeat(cosmo_hm *cmhm, peak_param *peak, int N, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  peak->printMode = 1; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  setMassSamplers(cmhm, peak, sampArr, err);                                                        forwardError(*err, __LINE__,);
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
  
  char name[STRING_LENGTH_MAX];
  int i;
  
  for (i=0; i<N; i++) {
    clock_t start = clock();
    printf("\n-- Making peak lists: %d / %d realizations\n", i+1, N);
    peakListFromMassFct(cmhm, peak, sampArr, hMap, galSamp, gMap, CCDMask, FFTSmoother, DCSmoother, kMap, variance, peakList, err); forwardError(*err, __LINE__,);
    //-- Output
    sprintf(name, "peakList_%d", i+1);
    outputPeakList(name, peak, peakList);
    routineTime(start, clock());
  }
  
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  printf("------------------------------------------------------------------------\n");
  return;
}
     
void doPeakList_withInputs(char fileName[], char fileName2[],char end[],cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
  char fpeakList[STRING_LENGTH_MAX];
  char fpeakListPos[STRING_LENGTH_MAX];
  char fpeakHist[STRING_LENGTH_MAX];

   sprintf(fpeakList, "peakList_%s",end);
   sprintf(fpeakListPos, "peakListPos_%s",end);
   sprintf(fpeakHist, "peakhist_%s",end);

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
 // printf(" read inpute 0 : \"%s\" \n", fileName);
  if ((fileName == NULL)||(fileName2 == NULL)) {
    //-- no input files
    	printf("Problem input files missing \n");
  }
  else {
  	printf(" Input file for halo : \"%s\" \n", fileName);
  	printf(" Input file for galaxies : \"%s\" \n", fileName2);
    read_halo_map(fileName, cmhm, hMap, err);                                forwardError(*err, __LINE__,);
    read_gal_map2(fileName2, cmhm,peak, gMap, err);                                forwardError(*err, __LINE__,);
  }
  makeMapAndOutputAll2(fileName,fileName2, cmhm, peak, gMap, FFTSmoother, DCSmoother, kMap, err); forwardError(*err, __LINE__,);

  computeLocalVariance_arr(peak, gMap, variance);
   
  if (peak->doNonlinear)       selectPeaks_mrlens("kappaMap_mrlens.fits", peak, gMap, peakList);
  else if (peak->DC_nbFilters) kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
  else                         kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);

  selectPeaks(peak, kMap, peakList, err);   forwardError(*err, __LINE__,);

  outputPeakList(fpeakList, peak, peakList);
  computePeaks2(fpeakListPos,peak,kMap,peakList,err);
  int silent = 1;
  makeHist(peakList, nuHist, silent);
  outputHist(fpeakHist, nuHist);


  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(nuHist);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doProduce_Catalog(char HaloFileName[],char GalFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
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


    //-- Carry out fast simulation
   sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
   setMassSamplers(cmhm, peak, sampArr, err);      forwardError(*err, __LINE__,);
   makeFastSimul(cmhm, peak, sampArr, hMap, err);    forwardError(*err, __LINE__,);
   outputFastSimul(HaloFileName, cmhm, peak, hMap);
   free_sampler_arr(sampArr);

   cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err);            forwardError(*err, __LINE__,);
   lensingCatalogueAndOutputAll2(GalFileName,cmhm,peak, hMap,gMap,err);

  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(nuHist);
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


		//-- Carry out fast simulation
	   sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
	   setMassSamplers(cmhm, peak, sampArr, err);      forwardError(*err, __LINE__,);
	   makeFastSimul(cmhm, peak, sampArr, hMap, err);    forwardError(*err, __LINE__,);
	   outputFastSimul(HaloFileName2, cmhm, peak, hMap);
	   free_sampler_arr(sampArr);

	   cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err);
       forwardError(*err, __LINE__,);
	   lensingCatalogueAndOutputAll2(GalFileName2,cmhm,peak, hMap,gMap,err);

	  free_halo_map(hMap);
	  free_sampler_t(galSamp);
	  free_gal_map(gMap);
	  free_short_mat(CCDMask);
	  free_FFT_arr(FFTSmoother);
	  free_FFT_arr(DCSmoother);
	  free_map_t(kMap);
	  free_FFT_arr(variance);
	  free_double_arr(peakList);
	  free_hist_t(nuHist);
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
	 // printf(" read inpute 0 : \"%s\" \n", fileName);
	  if ((fileName == NULL)||(fileName2 == NULL)) {
		//-- no input files
			printf("Problem input files missing \n");
	  }
	  else {
	  	printf(" Input file for halo : \"%s\" \n", HaloFileName2);
	  	printf(" Input file for galaxies : \"%s\" \n", GalFileName2);
		read_halo_map(HaloFileName2, cmhm, hMap, err);                                forwardError(*err, __LINE__,);
		read_gal_map(GalFileName2, cmhm,peak, gMap, err);                                forwardError(*err, __LINE__,);
	  }
	  makeMapAndOutputAll2(HaloFileName2,GalFileName2, cmhm, peak, gMap, FFTSmoother, DCSmoother, kMap, err); forwardError(*err, __LINE__,);

	  computeLocalVariance_arr(peak, gMap, variance);
	   
	  if (peak->doNonlinear)       selectPeaks_mrlens("kappaMap_mrlens.fits", peak, gMap, peakList);
	  else if (peak->DC_nbFilters) kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
	  else                         kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);

	  selectPeaks(peak, kMap, peakList, err);   forwardError(*err, __LINE__,);
	   outputPeakList(PeakListfich2, peak, peakList);

	  // computePeaks2("peakListPos",peak,kMap,peakList,err);

	  int silent = 1;
	  makeHist(peakList, nuHist, silent); 
	  outputHist(PeakHistfich2, nuHist);
	  //computePeaks2("TEST_TABLE_PEAK",peak,kMap,peakList,err);
	  free_halo_map(hMap);
	  free_sampler_t(galSamp);
	  free_gal_map(gMap);
	  free_short_mat(CCDMask);
	  free_FFT_arr(FFTSmoother);
	  free_FFT_arr(DCSmoother);
	  free_map_t(kMap);
	  free_FFT_arr(variance);
	  free_double_arr(peakList);
	  free_hist_t(nuHist);
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
	   setMassSamplers(cmhm, peak, sampArr, err); 
       forwardError(*err, __LINE__,);
	   makeFastSimul(cmhm, peak, sampArr, hMap, err);
       forwardError(*err, __LINE__,);
   	   sprintf(HaloFileName2, "%s_%3.3d",HaloFileName, i+1);
	   outputFastSimul_HOD(CmhmName,HaloFileName2, cmhm, peak, hMap);
	   free_sampler_arr(sampArr);
  	   free_halo_map(hMap);
	}
  return;
}


void doPeakList_withInputs_hod(char fileNameHal[], char fileNameGal[],char end[],cosmo_hm *cmhm, peak_param *peak, error **err)
{
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
  char fpeakList[STRING_LENGTH_MAX];
  char fpeakListPos[STRING_LENGTH_MAX];
  char fpeakHist[STRING_LENGTH_MAX];

  sprintf(fpeakList, "peakList_%s",end);
  sprintf(fpeakListPos, "peakListPos_%s",end);
  sprintf(fpeakHist, "peakhist_%s",end);

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

 // printf(" read inpute 0 : \"%s\" \n", fileName);
  if ((fileNameHal == NULL)||(fileNameGal == NULL)) {
    //-- no input files
    	printf("Problem input files missing \n");
  }
  else {
  	printf(" Input file for halo : \"%s\" \n", fileNameHal);
  	printf(" Input file for galaxies : \"%s\" \n", fileNameGal);
    read_halo_map(fileNameHal, cmhm, hMap, err); 
	forwardError(*err, __LINE__,);
    read_gal_map2(fileNameGal, cmhm,peak, gMap, err); 
    forwardError(*err, __LINE__,);
  }


  makeMapAndOutputAll2(fileNameHal,fileNameGal, cmhm, peak, gMap, FFTSmoother, DCSmoother, kMap, err);
  forwardError(*err, __LINE__,);

  computeLocalVariance_arr(peak, gMap, variance);
   
  if (peak->doNonlinear)       selectPeaks_mrlens("kappaMap_mrlens.fits", peak, gMap, peakList);
  else if (peak->DC_nbFilters) kappaToSNR_DC(peak, gMap, DCSmoother->array[0], kMap);
  else                         kappaToSNR_FFT(peak, gMap, FFTSmoother->array[0], kMap, variance->array[0]);

  lensingCatalogueAndOutputAll(cmhm, peak, hMap, gMap, err);                 forwardError(*err, __LINE__,);
  makeMapAndOutputAll(cmhm, peak, gMap, FFTSmoother, DCSmoother, kMap, err); forwardError(*err, __LINE__,);
  computeLocalVariance_arr(peak, gMap, variance);
   
  selectPeaks(peak, kMap, peakList, err);   forwardError(*err, __LINE__,);


  outputPeakList(fpeakList, peak, peakList);
  computePeaks2(fpeakListPos,peak,kMap,peakList,err);
  int silent = 1;
  makeHist(peakList, nuHist, silent);
  outputHist(fpeakHist, nuHist);

  //computePeaks2("TEST_TABLE_PEAK",peak,kMap,peakList,err);
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(nuHist);
  printf("------------------------------------------------------------------------\n");
  return;
}
