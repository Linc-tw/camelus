

  /***************************************
   **  peakSelection.c			**
   **  Chieh-An Lin, François Lanusse	**
   **  Version 2015.03.25		**
   ***************************************/


#include "peakSelection.h"


//----------------------------------------------------------------------
//-- Functions related to peak selection

int isPeak(double *kappa, int N1, int i, int j)
{
  int m, n;
  double centerValue = kappa[i+j*N1];
  for (m=-1; m<=1; m++) {
    for (n=-1; n<=1; n++) {
      if (m==0 && n==0) continue;
      if (kappa[(i+m)+(j+n)*N1] >= centerValue) return 0;
    }
  }
  return 1;
}

void selectPeaks(peak_param *peak, map_t *kMap, double_arr *peakList, error **err)
{
  int bufferSize = peak->bufferSize;
  testError(bufferSize<=0, peak_badValue, "Buffer size should be at least 1.", *err, __LINE__); forwardError(*err, __LINE__,);
  
  double *kappa = kMap->kappa;
  double *fArr = peakList->array;
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  
  int i, j, count = 0;
  for (j=bufferSize; j<N2-bufferSize; j++) {
    for (i=bufferSize; i<N1-bufferSize; i++) {
      if (isPeak(kappa, N1, i, j)) {
	fArr[count] = kappa[i+j*N1] / peak->sigma_noise;
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

void outputPeakList(char name[], peak_param *peak, double_arr *peakList)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Peak list\n");
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file, "# Filter = %s, theta_G = %g [arcmin], s = %g [pix]\n", STR_FILTER_T(peak->filter), peak->theta_G, peak->s);
  fprintf(file, "# sigma_eps = %g, sigma_pix = %g, sigma_noise = %g\n", peak->sigma_eps, peak->sigma_pix, peak->sigma_noise);
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

void peakListFromMassFct(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, gal_map *gMap, 
			 map_t *kMap, map_t *nMap, FFT_t *transformer, double_arr *peakList, error **err)
{
  makeFastSimul(cmhm, peak, sampArr, hMap, err);     forwardError(*err, __LINE__,);
  lensingForMap(cmhm, peak, hMap, gMap, err);        forwardError(*err, __LINE__,);
  makeMap(peak, gMap, kMap, nMap, transformer, err); forwardError(*err, __LINE__,);
  selectPeaks(peak, kMap, peakList, err);            forwardError(*err, __LINE__,);
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
  for (i=0; i<hist->length; i++) fprintf(file, "  %7.4f  %7.4f  %6d\n", hist->x_lower[i], hist->x_lower[i]+hist->dx, hist->n[i]);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

void inputHist(char name[], hist_t *hist, error **err)
{
  //-- TODO
  return;
}

//----------------------------------------------------------------------
//-- Functions related to chi2_t

chi2_t *initialize_chi2_t(int N, int d, error **err)
{
  //-- N              = number of data
  //-- d              = dimension of observable vector
  //-- *X_model       = observables from model, mean of several realizations for Camelus
  //-- *intermediate  = to stock invCov * (X_model - X_obs)
  //-- *X_obs         = observables from observation 
  //-- *cov           = covariance matrix for X_model, debiased
  //-- *cov2          = just a copy of cov
  //-- *invCov        = inverse of cov, debiased
  //-- *perm          = used for matrix inversion
  chi2_t *chichi       = (chi2_t*)malloc_err(sizeof(chi2_t), err); forwardError(*err, __LINE__,);
  chichi->N            = N;
  chichi->d            = d;
  chichi->X_model      = gsl_vector_alloc(d);
  chichi->intermediate = gsl_vector_alloc(d);
  chichi->X_obs        = gsl_vector_alloc(d);
  chichi->cov          = gsl_matrix_alloc(d, d);
  chichi->cov2         = gsl_matrix_alloc(d, d);
  chichi->invCov       = gsl_matrix_alloc(d, d);
  chichi->perm         = gsl_permutation_alloc(d);
  return chichi;
}

void free_chi2_t(chi2_t *chichi)
{
  gsl_vector_free(chichi->X_model);
  gsl_vector_free(chichi->intermediate);
  gsl_vector_free(chichi->X_obs);
  gsl_matrix_free(chichi->cov);
  gsl_matrix_free(chichi->cov2);
  gsl_matrix_free(chichi->invCov);
  gsl_permutation_free(chichi->perm);
  free(chichi);
  return;
}

void update_chi2_t(chi2_t *chichi, hist_t *obsHist, double *dataMat)
{
  //-- data should be n*d matrix
  int N = chichi->N;
  int d = chichi->d;
  gsl_vector *X_model = chichi->X_model;
  gsl_matrix *cov     = chichi->cov;
  gsl_matrix *cov2    = chichi->cov2;
  double value;
  
  //-- Update mean
  int i;
  for (i=0; i<d; i++) {
    value = gsl_stats_mean(dataMat+i*N, 1, N);
    gsl_vector_set(X_model, i, value);
    gsl_vector_set(chichi->X_obs, i, (double)(obsHist->n[i]));
  }
  
  //-- Update covariance
  int j;
  for (j=0; j<d; j++) {
    for (i=0; i<d; i++) {
      if (i < j) {
	value = gsl_matrix_get(cov, j, i);
	gsl_matrix_set(cov, i, j, value);
	continue;
      }
      value = gsl_stats_covariance_m(dataMat+i*N, 1, dataMat+j*N, 1, N, gsl_vector_get(X_model, i), gsl_vector_get(X_model, j));
      gsl_matrix_set(cov, i, j, value);
    } 
  }
  
  //-- Make a copy, because the invertion will destroy cov
  gsl_matrix_memcpy(cov2, cov); //-- Copy cov to cov2
  
  //-- Make LU decomposition and invert
  int s;
  gsl_linalg_LU_decomp(cov2, chichi->perm, &s);
  gsl_linalg_LU_invert(cov2, chichi->perm, chichi->invCov);
  
  //-- Debias the C_{ij}^{-1}
  //-- C^-1_nonbias = (N - d - 2) / (N - 1) * C^-1_bias
  double factor = (N - d - 2) / (double)(N - 1);
  gsl_matrix_scale(chichi->invCov, factor); //-- invCov *= factor
  return;
}

double execute_chi2_t(chi2_t *chichi)
{
  //-- Let Delta X = X_model - X_obs,
  //-- L = 1 / sqrt[(2 pi)^d * det(Cov)] * exp[-0.5 *(Delta X)^T * Cov^-1 * (Delta X)]
  //-- -2 ln L = -2 * [ -0.5 * ln (2 pi)^d - 0.5 * ln det(Cov) - 0.5 * (Delta X)^T * Cov^-1 * (Delta X) ]
  //--         = cst + ln det(Cov) + (Delta X)^T * Cov^-1 * (Delta X)
  //-- We set chi2 = ln det(Cov) + (Delta X)^T * Cov^-1 * (Delta X)
  
  //-- data should be N*d matrix
  int N = chichi->N;
  int d = chichi->d;
  gsl_vector *X_model = chichi->X_model;
  double value;
  
  gsl_vector_sub(X_model, chichi->X_obs); //-- X_model -= X_obs
  gsl_blas_dsymv(CblasUpper, 1.0, chichi->invCov, X_model, 0.0, chichi->intermediate); //-- intermediate = invCov * (X_model - X_obs)
  gsl_blas_ddot(X_model, chichi->intermediate, &value);
  value += gsl_linalg_LU_lndet(chichi->cov);
  return value;
}

//----------------------------------------------------------------------
//-- Main functions

void doPeakList(char KNMap[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  map_t *kMap = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  
  if (KNMap == NULL) {
    //-- Carry out KMap from fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
    setMassSamplers(cmhm, peak, sampArr, err);                                                        forwardError(*err, __LINE__,);
    halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
    gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
    makeGalaxies(cmhm, peak, gMap, err);                                                              forwardError(*err, __LINE__,);
    map_t *nMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
    FFT_t *transformer   = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
    fillGaussianKernel(transformer, peak->s);
    
    makeFastSimul(cmhm, peak, sampArr, hMap, err);     forwardError(*err, __LINE__,);
    lensingForMap(cmhm, peak, hMap, gMap, err);        forwardError(*err, __LINE__,);
    makeMap(peak, gMap, kMap, nMap, transformer, err); forwardError(*err, __LINE__,);
    
    outputFastSimul("haloList", cmhm, peak, hMap);
    outputGalaxies("galList", cmhm, peak, gMap);
    outputMap("KNMap", cmhm, peak, kMap);
    outputMap("NMap", cmhm, peak, nMap);
    
    free_sampler_arr(sampArr);
    free_halo_map(hMap);
    free_gal_map(gMap);
    free_map_t(nMap);
    free_FFT_t(transformer);
  }
  else {
    read_map_t(KNMap, peak, kMap, KN_map, err); forwardError(*err, __LINE__,);
  }
  
  int length   = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  int nBins    = 35;
  double lower = -5.00;
  double upper = 12.50;
  double_arr *peakList = initialize_double_arr(length);
  hist_t *hist         = initialize_hist_t(nBins);
  set_hist_t(hist, lower, upper);
  
  selectPeaks(peak, kMap, peakList, err); forwardError(*err, __LINE__,);
  outputPeakList("peakList", peak, peakList);
  int silent = 0;
  makeHist(peakList, hist, silent);
  outputHist("peakHist", hist);
  
  free_map_t(kMap);
  free_double_arr(peakList);
  free_hist_t(hist);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doPeakList_repeat(cosmo_hm *cmhm, peak_param *peak, int N, error **err)
{
  char name[STRING_LENGTH_MAX];
  int length = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  peak->printMode = 0; //-- Use 2 on planer
  
  sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
  setMassSamplers(cmhm, peak, sampArr, err);                                                        forwardError(*err, __LINE__,);
  halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  makeGalaxies(cmhm, peak, gMap, err);                                                              forwardError(*err, __LINE__,);
  map_t *kMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  map_t *nMap          = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_t *transformer   = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
  fillGaussianKernel(transformer, peak->s);
  double_arr *peakList = initialize_double_arr(length);
  
  int i;
  for (i=0; i<N; i++) {
    clock_t start = clock();
    printf("\n-- Making peak lists: %d / %d realizations\n", i+1, N);
    peakListFromMassFct(cmhm, peak, sampArr, hMap, gMap, kMap, nMap, transformer, peakList, err);   forwardError(*err, __LINE__,);
    //-- Output
    sprintf(name, "peakList_fast%d_gauss1.0", i+1);
    outputPeakList(name, peak, peakList);
    routineTime(start, clock());
  }
  
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  free_gal_map(gMap);
  free_map_t(kMap);
  free_map_t(nMap);
  free_FFT_t(transformer);
  free_double_arr(peakList);
  printf("------------------------------------------------------------------------\n");
  return;
}
/*
void doPeakHistForPMC(cosmo_hm *cmhm, peak_param *peak, hist_t *hist, double *dataMat, int N, error **err)
{
  //-- This function is called in chi2_wrapper for CosmoPMC.
  //-- Binning is given by hist. Results for N realizations are stocked in dataMat.
  
  char name[STRING_LENGTH_MAX];
  //int nbGal = (int)round(peak->n_gal * peak->area);
  int N1    = peak->resol[0] - 2 * peak->bufferSize;
  int N2    = peak->resol[1] - 2 * peak->bufferSize;
  
  sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
  setMassSamplers(cmhm, peak, sampArr, err);                                                       forwardError(*err, __LINE__,);
  halo_arr *halo       = initialize_halo_arr(NUMBER_HALOS_MAX, err);                               forwardError(*err, __LINE__,); //-- NUMBER_HALOS_MAX defined in peakParameters.h
  gal_arr *gal         = initialize_gal_arr(nbGal, err);                                           forwardError(*err, __LINE__,);
  makeGalaxies(cmhm, peak, gal, err);                                                              forwardError(*err, __LINE__,);
  gal_ptr_arr *ptr     = initialize_gal_ptr_arr(gal, err);                                         forwardError(*err, __LINE__,);
  gal_tree *tree       = initialize_gal_tree(ptr, 0, ptr->length, err);                            forwardError(*err, __LINE__,);
  pix_arr *map         = initialize_pix_arr(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  pix_arr *noise       = initialize_pix_arr(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  FFT_t *transformer   = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
  fillGaussianKernel(transformer, peak->s);
  double_arr *peakList = initialize_double_arr(N1 * N2);
  
  int i, j;
  for (i=0; i<N; i++) {
    clock_t start = clock();
    printf("\n-- Making peak histograms: %d / %d realizations\n", i+1, N);
    makePeakListFromBegin(cmhm, peak, sampArr, halo, gal, tree, map, noise, transformer, peakList, err); forwardError(*err, __LINE__,);
    int silent = 1;
    makeHist(peakList, hist, silent);
    for (j=0; j<hist->length; j++) dataMat[i+j*N] = (double)(hist->n[j]);
    routineTime(start, clock());
  }
  
  free_sampler_arr(sampArr);
  free_halo_arr(halo);
  free_gal_arr(gal);
  free_gal_ptr_arr(ptr);
  free_gal_tree(tree);
  free_pix_arr(map);
  free_pix_arr(noise);
  free_FFT_t(transformer);
  free_double_arr(peakList);
  printf("------------------------------------------------------------------------\n");
  return;
}
*/
//----------------------------------------------------------------------

