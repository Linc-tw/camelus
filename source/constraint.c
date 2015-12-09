

  /*******************************
   **  constraint.c		**
   **  Chieh-An Lin		**
   **  Version 2015.12.09	**
   *******************************/


#include "constraint.h"


//----------------------------------------------------------------------
//-- Functions related to data matrix

void fillMultiscale(peak_param *peak, hist_t *hist, double_mat *multiscale)
{
  int jNbin = peak->scaleInd * multiscale->N1;
  int i;
  for (i=0; i<hist->length; i++) multiscale->matrix[i+jNbin] = (double)hist->n[i];
  return;
}

void multiscaleFromMassFct(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask, 
			   FFT_arr *smoother, map_t *kMap, FFT_arr *variance, double_arr *peakList, hist_t *hist, hist_t *hist2, double_mat *multiscale, error **err)
{
  makeFastSimul(cmhm, peak, sampArr, hMap, err);                  forwardError(*err, __LINE__,);
  cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err); forwardError(*err, __LINE__,);
  lensingCatalogue(cmhm, peak, hMap, gMap, err);                  forwardError(*err, __LINE__,);
  makeMap(peak, gMap, smoother, kMap, err);                       forwardError(*err, __LINE__,);
  computeLocalVariance_arr(peak, gMap, variance);
  
  char name[STRING_LENGTH_MAX];
  int j;
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) {
    for (j=0; j<smoother->length; j++) {
      peak->scaleInd = j;
      kappaToSNR(peak, gMap, smoother->array[j], kMap, variance->array[j]);
      selectPeaks(peak, kMap, peakList, err);                     forwardError(*err, __LINE__,);
      makeHist(peakList, hist, 1);
      fillMultiscale(peak, hist, multiscale);
    }
  }
  if (peak->doSmoothing == 3 || peak->doSmoothing == 4) {
    peak->scaleInd = peak->nbLinFilters;
    sprintf(name, "kappaMap_mrlens_MPI%d.fits", peak->MPIInd);
    selectPeaks_mrlens(name, peak, gMap, peakList);
    makeHist(peakList, hist2, 1);
    fillMultiscale(peak, hist2, multiscale);
  }
  return;
}

void fillRealization(peak_param *peak, double_mat *multiscale, double *matrix)
{
  int j     = peak->realizationInd;
  int d_tot = multiscale->length;
  int i;
  for (i=0; i<d_tot; i++) matrix[i+j*d_tot] = multiscale->matrix[i];
  printf("data forwarded\n");
  return;
}

void outputMultiscale(char name[], double_mat *multiscale)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Multiscale data matrix\n");
  fprintf(file, "# N_nu      = %d\n", multiscale->N1);
  fprintf(file, "# nbFilters = %d\n", multiscale->N2);
  fprintf(file, "# d_tot     = %d\n", multiscale->length);
  fprintf(file, "#\n");
  fprintf(file, "# F0X0  F0X1  F0X2 ...\n");
  fprintf(file, "# F1X0  F1X1  F1X2 ...\n");
  fprintf(file, "# F2X0  F2X1  F2X2 ...\n");
  fprintf(file, "# ...   ...   ...  ...\n");
  
  int i, j;
  for (j=0; j<multiscale->N2; j++) {
    for (i=0; i<multiscale->N1; i++) {
      fprintf(file, " %5.1f ", multiscale->matrix[i+j*multiscale->N1]);
    }
    fprintf(file, "\n");
  }
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

void outputDataMat(char name[], double_mat *dataMat)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Multiscale data matrix\n");
  fprintf(file, "# d_tot = %d\n", dataMat->N1);
  fprintf(file, "# N     = %d\n", dataMat->N2);
  fprintf(file, "#\n");
  fprintf(file, "# R1F0X0  R1F0X1  ...  R1F1X0  R1F1X1  ...\n");
  fprintf(file, "# R2F1X0  R2F1X1  ...  R2F1X0  R2F1X1  ...\n");
  fprintf(file, "#  ...     ...    ...   ...     ...    ...\n");
  
  int i, j;
  for (j=0; j<dataMat->N2; j++) {
    for (i=0; i<dataMat->N1; i++) {
      fprintf(file, " %5.1f ", dataMat->matrix[i+j*dataMat->N1]);
    }
    fprintf(file, "\n");
  }
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}
/*
void outfitsCosmoParam(FITS_t *fits, cosmo_hm *cmhm, peak_param *peak)
{
  //-- Same as _paperI__outfitsCosmoParam
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "SIMUL",   peak->simulName,             "");
  addKeyword(fits, TDOUBLE, "OMEGAM",  &cmhm->cosmo->Omega_m,       "");
  addKeyword(fits, TDOUBLE, "OMEGADE", &cmhm->cosmo->Omega_de,      "");
  addKeyword(fits, TDOUBLE, "OMEGAB",  &cmhm->cosmo->Omega_b,       "");
  addKeyword(fits, TDOUBLE, "NS",      &cmhm->cosmo->n_spec,        "");
  addKeyword(fits, TDOUBLE, "H100",    &cmhm->cosmo->h_100,         "");
  addKeyword(fits, TDOUBLE, "SIGMA8",  &cmhm->cosmo->normalization, "");
  addKeyword(fits, TDOUBLE, "W0DE",    &cmhm->cosmo->w0_de,         "");
  addKeyword(fits, TDOUBLE, "W1DE",    &cmhm->cosmo->w1_de,         "");
  return;
}

void outfits_double_mat(FITS_t *fits, double_mat *dataMat, int N_bin, int nbFilters)
{
  char name[64];
  char name2[64];
  int d_tot = dataMat->N1; //-- d_tot = N_bin * nbFilters, N_bin = MAX(N_nu, N_kappa)
  int N     = dataMat->N2;
  int i, j, jNbin;
  
  for (j=0; j<nbFilters; j++) {
    jNbin = j * N_bin;
    for (i=0; i<N_bin; i++) {
      sprintf(name, "F%dX%d", j, i);
      addColumn(fits, name, TFLOAT, "-       ");
      sprintf(name, "TTYPE%d", i+1+jNbin);
      sprintf(name2, "Filter %d, bin %d", j, i);
      updateComment(fits, name, name2);
    }
  }
  
  addLineSpread(fits);
  addKeyword(fits, TINT, "DIM", &d_tot, "Dimension of data vector");
  addKeyword(fits, TINT, "N",   &N,     "Number of realizations");
  
  int jdtot;
  float relayE;
  
  for (j=0; j<N; j++) {
    jdtot = j * d_tot;
    for (i=0; i<d_tot; i++) {
      relayE = (float)(dataMat->matrix[i+jdtot]);
      writeTableColumn(fits, i, 1, &relayE);
    }
    nextRow(fits);
  }
  return;
}

void outfitsDataMat(char name[], cosmo_hm *cmhm, peak_param *peak, double_mat *dataMat)
{
  int N_bin     = peak->doSmoothing < 3 ? peak->N_nu : (peak->doSmoothing == 3 ? peak->N_kappa : MAX(peak->N_nu, peak->N_kappa));
  FITS_t *fits  = initializeTableWriter_FITS_t(name);
  outfits_double_mat(fits, dataMat, N_bin, peak->nbFilters);
  
  int resol2 = peak->resol[0] - 2 * peak->bufferSize;
  char name2[STRING_LENGTH_MAX];
  int i;
  
  addLineSpread(fits);
  addKeyword(fits, TDOUBLE,   "ZS",       &peak->z_s,                    "If ZS > 0, all source galaxies fixed at ZS");
  addKeyword(fits, TINT,      "RANDPOS",  &peak->doRandGalPos,           "0 = regular, 1 = random");
  addKeyword(fits, TDOUBLE,   "NGAL",     &peak->n_gal,                  "[arcmin^-2] Galaxy number density");
  addKeyword(fits, TDOUBLE,   "SIGEPS",   &peak->sigma_eps,              "[-] Ellipticity dispersion");
  addKeyword(fits, TINT,      "DOKAPPA",  &peak->doKappa,                "0 = gamma, 1 = kappa, 2 = g, 3 = g with linear KS");
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING,   "FIELD",    STR_FIELD_T(peak->field),      "");
  addKeyword(fits, TDOUBLE,   "OMEGA0",   &peak->Omega[0],               "[arcmin] Field width");
  addKeyword(fits, TDOUBLE,   "OMEGA1",   &peak->Omega[1],               "[arcmin] Field height");
  addKeyword(fits, TDOUBLE,   "THPIX",    &peak->theta_pix,              "[arcmin] Pixel size");
  
  addKeyword(fits, TINT,      "RESOL",    &peak->resol[0],               "[pix] Map resolution");
  addKeyword(fits, TINT,      "BUFF",     &peak->bufferSize,             "[pix] Border buffer size");
  addKeyword(fits, TINT,      "RESOL2",   &resol2,                       "[pix] Peak field resolution");
  
  addLineSpread(fits);
  addKeyword(fits, TINT,      "NBFILTER", &peak->nbFilters,              "Number of filters");
  for (i=0; i<peak->nbFilters; i++) {
    sprintf(name2, "FILTER%d", i);
    addKeyword(fits, TSTRING, name2,      STR_FILTER_T(peak->filter[i]), "Filter");
    sprintf(name2, "SCALE%d", i);
    addKeyword(fits, TDOUBLE, name2,      &peak->scale[i],               "[arcmin|-] Filter size in arcmin for linear, maximum scale number for nonlinear");
  }
  addKeyword(fits, TINT,      "NNU",      &dataMat->N1,                  "Number of bins");
  
  outfitsCosmoParam(fits, cmhm, peak);
    
  free_FITS_t(fits);
  printf("\"%s\" made\n", name);
  return;
}
*/
//----------------------------------------------------------------------
//-- Main function

void doMultiscale(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  char name[STRING_LENGTH_MAX];
  int N1_mask = (int)(peak->Omega[0] * peak->theta_CCD_inv);
  int N2_mask = (int)(peak->Omega[1] * peak->theta_CCD_inv);
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  int N_bin   = peak->doSmoothing < 3 ? peak->N_nu : (peak->doSmoothing == 3 ? peak->N_kappa : MAX(peak->N_nu, peak->N_kappa));
  peak->printMode = 2; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  sampler_arr *sampArr   = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  setMassSamplers(cmhm, peak, sampArr, err);                                                          forwardError(*err, __LINE__,);
  halo_map *hMap         = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp     = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                         forwardError(*err, __LINE__,);
  gal_map *gMap          = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask     = initialize_short_mat(N1_mask, N2_mask);
  if (peak->doMask == 1) fillMask_CFHTLenS_W1(peak, CCDMask);
  FFT_arr *smoother      = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernel(peak, smoother);
  map_t *kMap            = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_arr *variance      = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernelForVariance(peak, variance);
  double_arr *peakList   = initialize_double_arr(length);
  hist_t *hist           = initialize_hist_t(peak->N_nu);
  setHist(peak, hist);
  hist_t *hist2          = initialize_hist_t(peak->N_kappa);
  setHist2(peak, hist2);
  double_mat *multiscale = initialize_double_mat(N_bin, peak->nbFilters);
  
  //-- Peak histogram
  printf("\n-- Making a realization\n");
  multiscaleFromMassFct(cmhm, peak, sampArr, hMap, galSamp, gMap, CCDMask,
			smoother, kMap, variance, peakList, hist, hist2, multiscale, err); forwardError(*err, __LINE__,);
  
  //-- Output
  sprintf(name, "multiscale_nbFilters%d", peak->nbFilters);
  outputMultiscale(name, multiscale);
  
  //-- Free
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(smoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(hist);
  free_hist_t(hist2);
  free_double_mat(multiscale);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doDataMatrix(cosmo_hm *cmhm, peak_param *peak, int N, error **err)
{
  char name[STRING_LENGTH_MAX];
  int N1_mask = (int)(peak->Omega[0] * peak->theta_CCD_inv);
  int N2_mask = (int)(peak->Omega[1] * peak->theta_CCD_inv);
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  int N_bin   = peak->doSmoothing < 3 ? peak->N_nu : (peak->doSmoothing == 3 ? peak->N_kappa : MAX(peak->N_nu, peak->N_kappa));
  peak->printMode = 2; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  sampler_arr *sampArr  = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  setMassSamplers(cmhm, peak, sampArr, err);                                                          forwardError(*err, __LINE__,);
  halo_map *hMap         = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp     = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                         forwardError(*err, __LINE__,);
  gal_map *gMap          = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask     = initialize_short_mat(N1_mask, N2_mask);
  if (peak->doMask == 1) fillMask_CFHTLenS_W1(peak, CCDMask);
  FFT_arr *smoother      = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernel(peak, smoother);
  map_t *kMap            = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_arr *variance      = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernelForVariance(peak, variance);
  double_arr *peakList   = initialize_double_arr(length);
  hist_t *hist           = initialize_hist_t(peak->N_nu);
  setHist(peak, hist);
  hist_t *hist2          = initialize_hist_t(peak->N_kappa);
  setHist2(peak, hist2);
  double_mat *multiscale = initialize_double_mat(N_bin, peak->nbFilters);
  double_mat *dataMat    = initialize_double_mat(N_bin * peak->nbFilters, N);
  
  clock_t start = clock();
  int i;
  
  //-- Peak histogram loop
  printf("\n-- Making N = %d realizations\n", N);
  for (i=0; i<N; i++) {
    printf("Realization %4d: ", i+1);
    peak->realizationInd = i;
    multiscaleFromMassFct(cmhm, peak, sampArr, hMap, galSamp, gMap, CCDMask,
			  smoother, kMap, variance, peakList, hist, hist2, multiscale, err); forwardError(*err, __LINE__,);
    fillRealization(peak, multiscale, dataMat->matrix);
  }
  
  double duration = (double)(clock() - start) / CLOCKS_PER_SEC;
  printf("%d secs for %d realizations, %.2f secs per run\n", (int)duration, N, duration/N);
  printf("\n");
  
  //-- Output
  sprintf(name, "dataMat_nbFilters%d_N%d", peak->nbFilters, N);
  outputDataMat(name, dataMat);
  
  //-- Free
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(smoother);
  free_map_t(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(hist);
  free_hist_t(hist2);
  free_double_mat(multiscale);
  free_double_mat(dataMat);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------


/*
//----------------------------------------------------------------------
//-- Functions related to likelihood_t

likelihood_t *initialize_likelihood_t(int d, int N, error **err)
{
  //-- d             = dimension of data vector
  //-- N             = number of data
  //-- matrix        = d * N data matrix
  //-- *x_mod        = Model prediction, mean of several realizations of Camelus
  //-- *x_obs        = Data vector from observation
  //-- *intermediate = to stock invCov * (x^mod - x^obs)
  //-- *cov          = covariance matrix for x^mod, debiased
  //-- *cov2         = just a copy of cov
  //-- *invCov       = inverse of cov, debiased
  //-- *perm         = used for matrix inversion
  //-- **estArr;     = Array of kernel density estimators
  likelihood_t *LLH = (likelihood_t*)malloc_err(sizeof(likelihood_t), err); forwardError(*err, __LINE__,);
  LLH->d            = d;
  LLH->N            = N;
  LLH->matrix       = (double*)malloc_err(d * N * sizeof(double), err);     forwardError(*err, __LINE__,);
  LLH->x_mod        = gsl_vector_alloc(d);
  LLH->x_obs        = gsl_vector_alloc(d);
  LLH->cov          = gsl_matrix_alloc(d, d);
  LLH->cov2         = gsl_matrix_alloc(d, d);
  LLH->invCov       = gsl_matrix_alloc(d, d);
  LLH->perm         = gsl_permutation_alloc(d);
  LLH->estArr       = initialize_KDE_arr(d, N);
  
  LLH->intermediate = gsl_vector_alloc(d);
  LLH->B_arr        = gsl_vector_alloc(d);
  return LLH;
}

void free_likelihood_t(likelihood_t *LLH)
{
  if (LLH->matrix) {free(LLH->matrix); LLH->matrix = NULL;}
  gsl_vector_free(LLH->x_mod);
  gsl_vector_free(LLH->x_obs);
  gsl_matrix_free(LLH->cov);
  gsl_matrix_free(LLH->cov2);
  gsl_matrix_free(LLH->invCov);
  gsl_permutation_free(LLH->perm);
  if (LLH->estArr) {free_KDE_arr(LLH->estArr); LLH->estArr = NULL;}
  
  gsl_vector_free(LLH->intermediate);
  gsl_vector_free(LLH->B_arr);
  
  free(LLH); LLH = NULL;
  return;
}

void set_likelihood_t(likelihood_t *LLH)
{
  //-- Set x_mod, cov, invCov
  
  int d = LLH->d;
  int N = LLH->N;
  double *matrix    = LLH->matrix;
  gsl_vector *x_mod = LLH->x_mod;
  gsl_matrix *cov   = LLH->cov;
  gsl_matrix *cov2  = LLH->cov2;
  double value;
  
  //-- Set x^mod
  int i;
  for (i=0; i<d; i++) {
    value = gsl_stats_mean(matrix+i, d, N);
    gsl_vector_set(x_mod, i, value);
  }
  
  //-- Set cov
  int j;
  for (j=0; j<d; j++) {
    for (i=0; i<d; i++) {
      if (i < j) {
	value = gsl_matrix_get(cov, j, i);
	gsl_matrix_set(cov, i, j, value);
	continue;
      }
      value = gsl_stats_covariance_m(matrix+i, d, matrix+j, d, N, gsl_vector_get(x_mod, i), gsl_vector_get(x_mod, j));
      gsl_matrix_set(cov, i, j, value);
    } 
  }
  
  //-- Make a copy, because the invertion will destroy cov
  gsl_matrix_memcpy(cov2, cov); //-- Copy cov to cov2
  
  //-- Make LU decomposition and invert
  int s;
  gsl_linalg_LU_decomp(cov2, LLH->perm, &s);
  gsl_linalg_LU_invert(cov2, LLH->perm, LLH->invCov);
  
  //-- Debias the Cov_{ij}^{-1}
  //-- Cov^-1_nonbias = (N - d - 2) / (N - 1) * Cov^-1_bias
  double factor = (N - d - 2) / (double)(N - 1);
  gsl_matrix_scale(LLH->invCov, factor); //-- invCov *= factor
  return;
}

void setKDE_likelihood_t(likelihood_t *LLH)
{
  int d = LLH->d;
  int N = LLH->N;
  int i, j;
  for (i=0; i<d; i++) {
    for (j=0; j<N; j++) LLH->estArr->array[i]->sample[j] = LLH->matrix[i+j*d];
    set_KDE_t(LLH->estArr->array[i], -1); //-- Silverman's rule
  }
  return;
}

void fillObservation_likelihood_t(char name[], likelihood_t *LLH, error **err)
{
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2;
  
  //-- Read comments
  int c = fgetc(file);
  while (c == (int)'#') {
    buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    c = fgetc(file);
  }
  ungetc(c, file);
  
  //-- Read the data vector dimension
  int nbBins;
  buffer2 = fscanf(file, "%d\n", &nbBins);
  testErrorRet(nbBins!=LLH->d, peak_match, "Data size error", *err, __LINE__,);
  
  //-- Read observation data
  double value;
  int i;
  for (i=0; i<nbBins; i++) {
    buffer2 = fscanf(file, "%lf\n", &value);
    gsl_vector_set(LLH->x_obs, i, value);
  }
  
  fclose(file);
  printf("Observation data read\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to chi^2 computation

double computeGaussianChi2(likelihood_t *LLH)
{
  //-- Let DeltaX = x^mod - x^obs,
  //-- L = 1 / sqrt[(2 pi)^d * det(Cov)] * exp[-0.5 * DeltaX^T * Cov^-1 * DeltaX]
  //-- -2 ln L = -2 * [ -0.5 * ln (2 pi)^d - 0.5 * ln det(Cov) - 0.5 * DeltaX^T * Cov^-1 * DeltaX ]
  //--         = cst + ln det(Cov) + DeltaX^T * Cov^-1 * DeltaX
  //-- We set chi2 = ln det(Cov) + DeltaX^T * Cov^-1 * DeltaX
  
  int d = LLH->d;
  int N = LLH->N;
  gsl_vector *x_mod = LLH->x_mod;
  double value;
  
  gsl_vector_sub(x_mod, LLH->x_obs); //-- x^mod -= x^obs
  gsl_blas_dsymv(CblasUpper, 1.0, LLH->invCov, x_mod, 0.0, LLH->intermediate); //-- intermediate = invCov * (x^mod - x^obs)
  gsl_blas_ddot(x_mod, LLH->intermediate, &value);
  value += gsl_linalg_LU_lndet(LLH->cov);
  return value;
}

double computeCopulaChi2(likelihood_t *LLH)
{
  //-- Let 
  //--    mu_i = mean_{i=1..d} x^model_i
  //--     q_i = gauss_cdf^-1( cdf_i(x^obs_i) )
  //--   pdf_i = marginalized pdf for x^obs_i
  //--   cdf_i = marginalized cdf for x^obs_i
  //--  
  //-- -2 ln L = cst + ln det(cov) + [sum_ij (q_i - mu_i) invCov_ij (q_j - mu_j)] - 2 [sum_i ln(sigma_i)] - [sum_i (q_i - mu_i)^2 / sigma_i^2] - 2 [sum_i ln( pdf_i(x_obs_i) )]
  //--         = cst + A + [sum_ij B_i invCov_ij B_j] - 2 [sum_i C_i] - [sum_i D_i^2] - 2 [sum_i E_i]
  //-- where
  //--   A   = ln det(cov)
  //--   B_i = q_i - mu_i
  //--   C_i = ln(sigma_i)
  //--   D_i = (q_i - mu_i) / sigma_i
  //--   E_i = ln( pdf_i(x^obs_i) )
  
  int d = LLH->d;
  int N = LLH->N;
  gsl_vector *x_mod = LLH->x_mod;
  gsl_vector *x_obs = LLH->x_obs;
  
  double value, sigma, x;
  double A = gsl_linalg_LU_lndet(LLH->cov);
  double C = 0.0;
  double D = 0.0;
  double E = 0.0;
  
  int i;
  for (i=0; i<d; i++) {
    //-- Compute C
    sigma = gsl_matrix_get(LLH->cov, i, i); //-- sigma_i^2
    sigma = sqrt(sigma);                    //-- sigma_i
    C    += log(sigma);                     //-- ln(sigma_i) => sum_i C_i
    
    //-- Compute D
    x     = gsl_vector_get(x_obs, i);
    value = integrate_KDE_t(LLH->estArr->array[i], x); //-- cdf(x^obs_i)
    value = gsl_cdf_ugaussian_Pinv(value);             //-- (q_i - mu_i) / sigma_i
    D    += pow(value, 2);                             //-- (q_i - mu_i)^2 / sigma_i^2 => sum_i D_i^2
    
    //-- Compute B
    gsl_vector_set(LLH->B_arr, i, value * sigma); //-- q_i - mu_i
    
    //-- Compute E
    value = execute_KDE_t(LLH->estArr->array[i], x); //-- pdf_i(x^obs_i)
    E    += log(value);                              //-- ln( pdf_i(x^obs_i) ) => sum_i E_i
  }
  
  gsl_blas_dsymv(CblasUpper, 1.0, LLH->invCov, LLH->B_arr, 0.0, LLH->intermediate); //-- intermediate = invCov * (x_mod - x_obs)
  gsl_blas_ddot(LLH->B_arr, LLH->intermediate, &value);                             //-- value = [sum_ij B_i invCov_ij B_j]
  
  value += A - 2*C - D - 2*E;
  return value;
}

void doChi2(cosmo_hm *cmhm, peak_param *peak, int N, error **err)
{
  char name[STRING_LENGTH_MAX];
  int N1_mask = (int)(peak->Omega[0] * peak->theta_CCD_inv);
  int N2_mask = (int)(peak->Omega[1] * peak->theta_CCD_inv);
  int length  = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  peak->printMode = 2; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  //-- Initialization
  sampler_arr *sampArr   = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  setMassSamplers(cmhm, peak, sampArr, err);                                                          forwardError(*err, __LINE__,);
  halo_map *hMap         = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp     = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                         forwardError(*err, __LINE__,);
  gal_map *gMap          = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask     = initialize_short_mat(N1_mask, N2_mask);
  if (peak->doMask == 1) fillMask_CFHTLenS_W1(peak, CCDMask);
  FFT_arr *smoother      = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernel(peak, smoother);
  map_t *kMap            = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_arr *varArr        = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernelForVariance(peak, varArr);
  double_arr *peakList   = initialize_double_arr(length);
  hist_t *hist           = initialize_hist_t(peak->N_nu);
  setHist(peak, hist);
  hist_t *hist2          = initialize_hist_t(peak->N_nu);
  setHist(peak, hist2); //TODO
  double_mat *multiscale = initialize_double_mat(peak->N_nu, peak->nbFilters);
  likelihood_t *LLH      = initialize_likelihood_t(peak->N_nu * peak->nbFilters, N, err);             forwardError(*err, __LINE__,);
  //sprintf(name, ); //TODO
  fillObservation_likelihood_t("../demo/peakHist", LLH, err);                                         forwardError(*err, __LINE__,);
  
  //-- Peak histogram loop
  int silent = 1;
  printf("\n-- Making N = %d realizations\n", N);
  int i;
  for (i=0; i<N; i++) {
    printf("Realization %4d: ", i+1);
    peak->realizationInd = i;
    multiscalePeakHistFromMassFct(cmhm, peak, sampArr, hMap, galSamp, gMap, CCDMask,
				  smoother, kMap, varArr, peakList, hist, hist2, multiscale, err); forwardError(*err, __LINE__,);
    fillRealization(peak, multiscale, LLH->matrix);
    printf("\n");
  }
  printf("\n");
  
  //-- Computation
  set_likelihood_t(LLH);
  setKDE_likelihood_t(LLH);
  double L_vg = computeGaussianChi2(LLH);
  double L_vc = computeCopulaChi2(LLH);
  printf("-2 ln L_vg = %e\n", L_vg);
  printf("-2 ln L_vc = %e\n", L_vc);
  
  //-- Free
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  free_FFT_arr(smoother);
  free_map_t(kMap);
  free_FFT_arr(varArr);
  free_double_arr(peakList);
  free_hist_t(hist);
  free_double_mat(multiscale);
  free_likelihood_t(LLH);
  printf("------------------------------------------------------------------------\n");
  return;
}
*/
