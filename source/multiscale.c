

  /*******************************************************
   **  multiscale.c					**
   **  Version 2018.03.14				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "multiscale.h"


//----------------------------------------------------------------------
//-- Functions related to data matrix

void addToMultiscale(hist_t *hist, double_mat *multiscale, int scaleInd)
{
  //-- Need to initialize multiscale
  int jNbin = scaleInd * multiscale->N1;
  int i;
  for (i=0; i<hist->length; i++) multiscale->matrix[i+jNbin] += (double)hist->n[i];
  return;
}

void mapToMultiscale_FFT(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err)
{
  int j;
  for (j=0; j<pkPar->FFT_nbFilters; j++) {
    kappaToSNR_FFT(pkPar, gMap, FFTSmoother->array[j], kMap, variance->array[j], j);
    selectPeaks(pkPar, kMap, peakList, err); forwardError(*err, __LINE__,);
    makeHist(peakList, nuHist, 0); //-- verbose = 0
    addToMultiscale(nuHist, multiscale, j);
  }
  return;
}

void mapToMultiscale_DC(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err)
{
  int W = pkPar->FFT_nbFilters;
  int j;
  for (j=0; j<pkPar->DC_nbFilters; j++) {
    //-- DC
    kappaToSNR_DC(pkPar, DCSmoother->array[j], kMap, j);
    selectPeaks(pkPar, kMap, peakList, err); forwardError(*err, __LINE__,);
    makeHist(peakList, nuHist, 0); //-- verbose = 0
    addToMultiscale(nuHist, multiscale, j+W);
  }
  return;
}

void mapToMultiscaleAndOutput_FFT(peak_param *pkPar, gal_map *gMap, FFT_arr *FFTSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err)
{
  char name[STRING_LENGTH_MAX];
  int j;
  
  for (j=0; j<pkPar->FFT_nbFilters; j++) {
    //-- FFT
    kappaToSNR_FFT(pkPar, gMap, FFTSmoother->array[j], kMap, variance->array[j], j);
    selectPeaks(pkPar, kMap, peakList, err); forwardError(*err, __LINE__,);
    
    if (pkPar->doFITS == 0) {
      sprintf(name, "%speakList_FFT%d", pkPar->prefix, j);
      outAsciiPeakList(name, pkPar, peakList, j, err); forwardError(*err, __LINE__,);
    }
    else {
      sprintf(name, "%speakList_FFT%d.fits", pkPar->prefix, j);
      outFitsPeakList(name, pkPar, peakList, j);
    }
    
    makeHist(peakList, nuHist, 0); //-- verbose = 0
    
    if (pkPar->doFITS == 0) {
      sprintf(name, "%speakHist_FFT%d", pkPar->prefix, j);
      outAsciiHist(name, pkPar, nuHist, j, err); forwardError(*err, __LINE__,);
    }
    else {
      sprintf(name, "%speakHist_FFT%d.fits", pkPar->prefix, j);
      outFitsHist(name, pkPar, nuHist, j);
    }
    
    addToMultiscale(nuHist, multiscale, j);
  }
  return;
}

void mapToMultiscaleAndOutput_DC(peak_param *pkPar, gal_map *gMap, FFT_arr *DCSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, double_mat *multiscale, error **err)
{
  int W = pkPar->FFT_nbFilters;
  char name[STRING_LENGTH_MAX];
  int j;
  
  for (j=0; j<pkPar->DC_nbFilters; j++) {
    //-- DC
    kappaToSNR_DC(pkPar, DCSmoother->array[j], kMap, j);
    selectPeaks(pkPar, kMap, peakList, err); forwardError(*err, __LINE__,);
    
    if (pkPar->doFITS == 0) {
      sprintf(name, "%speakList_DC%d", pkPar->prefix, j);
      outAsciiPeakList(name, pkPar, peakList, j+W, err); forwardError(*err, __LINE__,);
    }
    else {
      sprintf(name, "%speakList_DC%d.fits", pkPar->prefix, j);
      outFitsPeakList(name, pkPar, peakList, j+W);
    }
    
    makeHist(peakList, nuHist, 0); //-- verbose = 0
    
    if (pkPar->doFITS == 0) {
      sprintf(name, "%speakHist_DC%d", pkPar->prefix, j);
      outAsciiHist(name, pkPar, nuHist, j+W, err); forwardError(*err, __LINE__,);
    }
    else {
      sprintf(name, "%speakHist_DC%d.fits", pkPar->prefix, j);
      outFitsHist(name, pkPar, nuHist, j+W);
    }
  
    addToMultiscale(nuHist, multiscale, j+W);
  }
  return;
}

void massFctToMultiscale(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, interpolator_t *k1Inter, halo_map *hMap, sampler_t *gSamp, gal_map *gMap, mask_map *mask, 
			 FFT_arr *FFTSmoother, FFT_arr *DCSmoother, signal_map *kMap, FFT_arr *variance, double_arr *peakList, hist_t *nuHist, 
			 double_mat *multiscale, int nbPatches, error **err)
{
  reset_double(multiscale->matrix, multiscale->length); //-- Reset
  int i;
  
  for (i=0; i<nbPatches; i++) {
    makeFastSimul(chPar, pkPar, hSampArr, hMap, err);                                                 forwardError(*err, __LINE__,);
    cleanOrMakeOrResample(chPar, pkPar, gSamp, gMap, mask, err);                                      forwardError(*err, __LINE__,);
    lensingCatalogue(chPar, pkPar, hMap, gMap, k1Inter, err);                                         forwardError(*err, __LINE__,);
    makeMaps(pkPar, gMap, FFTSmoother, DCSmoother, kMap, err);                                        forwardError(*err, __LINE__,);
    makeLocalVariance(pkPar, gMap, variance);
    mapToMultiscale_FFT(pkPar, gMap, FFTSmoother, kMap, variance, peakList, nuHist, multiscale, err); forwardError(*err, __LINE__,);
    mapToMultiscale_DC(pkPar, gMap, DCSmoother, kMap, variance, peakList, nuHist, multiscale, err);   forwardError(*err, __LINE__,);
  }
  return;
}

void forwardResult(peak_param *pkPar, double_mat *multiscale, double *matrix)
{
  int j     = pkPar->realizationInd;
  int d_tot = multiscale->length;
  int i;
  for (i=0; i<d_tot; i++) matrix[i+j*d_tot] = multiscale->matrix[i];
  if (pkPar->verbose < 3) printf("Forwarded result\n");
  else if (pkPar->verbose <= 5) printf("forwarded result\n");
  return;
}

void inverseCovarianceMatrix(double_mat *dataMat, gsl_matrix *invCov)
{
  int d_tot = dataMat->N1;
  int N     = dataMat->N2;
  
  gsl_vector *mean      = gsl_vector_alloc(d_tot);
  gsl_matrix *cov       = gsl_matrix_alloc(d_tot, d_tot);
  gsl_matrix *cholesky  = gsl_matrix_alloc(d_tot, d_tot);
  gsl_permutation *perm = gsl_permutation_alloc(d_tot);
  
  int i, j, jdTot;
  
  for (i=0; i<d_tot; i++) mean->data[i] = gsl_stats_mean(dataMat->matrix, d_tot, N);
  for (j=0; j<d_tot; j++) {
    jdTot = j * d_tot;
    for (i=0; i<d_tot; i++) {
      if (i < j) cov->data[i+jdTot] = cov->data[j+i*d_tot];
      cov->data[i+jdTot] = gsl_stats_covariance_m(dataMat->matrix+i, d_tot, dataMat->matrix+j, d_tot, N, mean->data[i], mean->data[j]);
    }
  }
  
  //-- Make LU decomposition and invert
  int s;
  gsl_linalg_LU_decomp(cov, perm, &s);
  gsl_linalg_LU_invert(cov, perm, invCov);
  
  //-- Debias the C_{ij}^{-1}
  //-- C^-1_unbiased = (N - d_tot - 2) / (N - 1) * C^-1_biased
  double factor = (double)(N - d_tot - 2) / (double)(N - 1);
  gsl_matrix_scale(invCov, factor);
  
  gsl_vector_free(mean);
  gsl_matrix_free(cov);
  gsl_matrix_free(cholesky);
  gsl_permutation_free(perm);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to output

void outAsciiMultiscale(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *multiscale, error **err)
{
  if (pkPar->outMultiscale == 0) return;
  
  int N_nu = pkPar->N_nu;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Multiscale counts\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiPeakField(file, pkPar);
  outAsciiAllFilterInfo(file, pkPar);
  outAsciiHistInfo(file, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# d_tot = %d\n", pkPar->N_nu * pkPar->nbFilters);
  fprintf(file, "#\n");
  fprintf(file, "# F0X0  F0X1  F0X2 ...\n");
  fprintf(file, "# F1X0  F1X1  F1X2 ...\n");
  fprintf(file, "# F2X0  F2X1  F2X2 ...\n");
  fprintf(file, "# ...   ...   ...  ...\n");
  
  int i, j, jN_nu;
  
  for (j=0; j<pkPar->nbFilters; j++) {
    jN_nu = j * N_nu;
    for (i=0; i<N_nu; i++) {
      fprintf(file, " %6.1f ", multiscale->matrix[i+jN_nu]);
    }
    fprintf(file, "\n");
  }
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

void outFitsMultiscale(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *multiscale)
{
  if (pkPar->outMultiscale == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  
  int N_nu = pkPar->N_nu;
  char name2[STRING_LENGTH_MAX];
  int i;
  
  for (i=0; i<N_nu; i++) {
    sprintf(name2, "X%d", i);
    addColumn(fits, name2, TFLOAT, "-       ");
  }
  
  int j, jN_nu;
  float fBuff;
  
  for (j=0; j<pkPar->nbFilters; j++) {
    jN_nu = j * N_nu;
    for (i=0; i<N_nu; i++) {
      fBuff = (float)(multiscale->matrix[i+jN_nu]); writeTableColumn(fits, i, 1, &fBuff);
    }
    nextRow(fits);
  }
  
  int iBuff;
  
  addLineSpread(fits);
  iBuff = N_nu * pkPar->nbFilters;
  addKeyword(fits, TINT, "DTOT", &iBuff, "Dimension of data vector");
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsPeakField(fits, pkPar);
  outFitsAllFilterInfo(fits, pkPar);
  outFitsHistInfo(fits, pkPar);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

void outputMultiscale(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *multiscale, error **err)
{
  char name2[STRING_LENGTH_MAX];
  if (pkPar->doFITS == 0) {
    outAsciiMultiscale(name, chPar, pkPar, multiscale, err);
    forwardError(*err, __LINE__,);
  }
  else {
    sprintf(name2, "%s.fits", name);
    outFitsMultiscale(name2, chPar, pkPar, multiscale);
  }
  return;
}

void outAsciiDataMat(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *dataMat, error **err)
{
  int d_tot = dataMat->N1; //-- d_tot = N_nu * nbFilters
  int N     = dataMat->N2;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Data matrix\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiPeakField(file, pkPar);
  outAsciiAllFilterInfo(file, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# d_tot = %d\n", d_tot);
  fprintf(file, "# N     = %d\n", N);
  fprintf(file, "#\n");
  fprintf(file, "# R1F0X0  R1F0X1  ...  R1F1X0  R1F1X1  ...\n");
  fprintf(file, "# R2F1X0  R2F1X1  ...  R2F1X0  R2F1X1  ...\n");
  fprintf(file, "#  ...     ...    ...   ...     ...    ...\n");
  
  int i, j, jd_tot;
  
  for (j=0; j<N; j++) {
    jd_tot = j * d_tot;
    for (i=0; i<d_tot; i++) {
      fprintf(file, " %5.1f ", dataMat->matrix[i+jd_tot]);
    }
    fprintf(file, "\n");
  }
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

void outFitsDataMat(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *dataMat)
{
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits  = initializeTableWriter_FITS_t(name);
  
  char name2[STRING_LENGTH_MAX];
  int d_tot = dataMat->N1; //-- d_tot = N_nu * nbFilters
  int i;
  
  for (i=0; i<d_tot; i++) {
    sprintf(name2, "X%d", i);
    addColumn(fits, name2, TFLOAT, "-       ");
  }
  
  int N = dataMat->N2;
  int j, jdtot;
  float fBuff;
  
  for (j=0; j<N; j++) {
    jdtot = j * d_tot;
    for (i=0; i<d_tot; i++) {
      fBuff = (float)(dataMat->matrix[i+jdtot]);
      writeTableColumn(fits, i, 1, &fBuff);
    }
    nextRow(fits);
  }
  
  addLineSpread(fits);
  addKeyword(fits, TINT, "DTOT", &d_tot, "Dimension of data vector");
  addKeyword(fits, TINT, "N",    &N,     "Number of realizations");
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsPeakField(fits, pkPar);
  outFitsAllFilterInfo(fits, pkPar);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

void outputDataMat(char name[], cosmo_hm *chPar, peak_param *pkPar, double_mat *dataMat, error **err)
{
  char name2[STRING_LENGTH_MAX];
  if (pkPar->doFITS == 0) {
    outAsciiDataMat(name, chPar, pkPar, dataMat, err);
    forwardError(*err, __LINE__,);
  }
  else {
    sprintf(name2, "%s.fits", name);
    outFitsDataMat(name2, chPar, pkPar, dataMat);
  }
  return;
}

void outAsciiInvCov(char name[], cosmo_hm *chPar, peak_param *pkPar, gsl_matrix *invCov, error **err)
{
  int d_tot = invCov->size1;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Inverse covariance matrix\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiPeakField(file, pkPar);
  outAsciiAllFilterInfo(file, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# d_tot = %d\n", d_tot);
  fprintf(file, "#\n");
  
  int i, j, jd_tot;
  
  for (j=0; j<d_tot; j++) {
    jd_tot = j * d_tot;
    for (i=0; i<d_tot; i++) {
      fprintf(file, " %12.5e", invCov->data[i+jd_tot]);
    }
    fprintf(file, "\n");
  }
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

void outFitsInvCov(char name[], cosmo_hm *chPar, peak_param *pkPar, gsl_matrix *invCov)
{
#ifdef __CAMELUS_USE_FITS__
  int d_tot = invCov->size1;
  FITS_t *fits = initializeImageWriter_FITS_t(name);
  write2DImage_double(fits, TDOUBLE, DOUBLE_IMG, d_tot, d_tot, (void*)invCov->data);
  
  addLineSpread(fits);
  addKeyword(fits, TINT, "DTOT", &d_tot, "Dimension of data vector");
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath, "[-] Path of the peak parameter file");
  outFitsFieldInfo(fits, pkPar);
  outFitsPeakField(fits, pkPar);
  outFitsAllFilterInfo(fits, pkPar);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

void outputInvCov(char name[], cosmo_hm *chPar, peak_param *pkPar, gsl_matrix *invCov, error **err)
{
  char name2[STRING_LENGTH_MAX];
  if (pkPar->doFITS == 0) {
    outAsciiInvCov(name, chPar, pkPar, invCov, err);
    forwardError(*err, __LINE__,);
  }
  else {
    sprintf(name2, "%s.fits", name);
    outFitsInvCov(name2, chPar, pkPar, invCov);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to pipeline_t

pipeline_t *initialize_pipeline_t(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  //-- Camelus pipeline
  
  pipeline_t *pipe  = (pipeline_t*)malloc_err(sizeof(pipeline_t), err);                               forwardError(*err, __LINE__, NULL);
  pipe->hSampArr    = initialize_sampler_arr(pkPar->N_z_halo, pkPar->N_M+1);
  pipe->k1Inter     = initialize_interpolator_t(pkPar->N_z_gal+1);
  setMassSampAndK1Inter(chPar, pkPar, pipe->hSampArr, NULL, pipe->k1Inter, err);                      forwardError(*err, __LINE__, NULL); //-- lambda = NULL
  pipe->hMap        = initialize_halo_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);   forwardError(*err, __LINE__, NULL);
  pipe->gSamp       = initialize_sampler_t(pkPar->N_z_gal);
  setGalaxySampler(chPar, pkPar, pipe->gSamp, err);                                                   forwardError(*err, __LINE__, NULL);
  pipe->gMap        = initialize_gal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);    forwardError(*err, __LINE__, NULL);
  pipe->mask        = initializeMask(pkPar, err);                                                     forwardError(*err, __LINE__, NULL);
  pipe->FFTSmoother = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  makeKernel(pkPar, pipe->FFTSmoother);
  pipe->DCSmoother  = initialize_FFT_arr(MAX(pkPar->DC_nbFilters, 1), pkPar->FFTSize);
  pipe->kMap        = initialize_signal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err); forwardError(*err, __LINE__, NULL);
  pipe->variance    = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  makeKernelForVariance(pkPar, pipe->variance);
  pipe->peakList    = initialize_double_arr(pkPar->peakListLength);
  pipe->nuHist      = initialize_hist_t(pkPar->N_nu);
  setHist_nu(pkPar, pipe->nuHist);
  pipe->multiscale  = initialize_double_mat(pkPar->N_nu, pkPar->nbFilters);
  
  if (pkPar->MPIInd == 0) printf("Initialized the pipeline\n");
  return pipe;
}

void free_pipeline_t(pipeline_t *pipe)
{
  if (pipe) {
    if (pipe->hSampArr)    {free_sampler_arr(pipe->hSampArr);   pipe->hSampArr    = NULL;}
    if (pipe->k1Inter)     {free_interpolator_t(pipe->k1Inter); pipe->k1Inter     = NULL;}
    if (pipe->hMap)        {free_halo_map(pipe->hMap);          pipe->hMap        = NULL;}
    if (pipe->gSamp)       {free_sampler_t(pipe->gSamp);        pipe->gSamp       = NULL;}
    if (pipe->gMap)        {free_gal_map(pipe->gMap);           pipe->gMap        = NULL;}
    if (pipe->mask)        {free_mask_map(pipe->mask);          pipe->mask        = NULL;}
    if (pipe->FFTSmoother) {free_FFT_arr(pipe->FFTSmoother);    pipe->FFTSmoother = NULL;}
    if (pipe->DCSmoother)  {free_FFT_arr(pipe->DCSmoother);     pipe->DCSmoother  = NULL;}
    if (pipe->kMap)        {free_signal_map(pipe->kMap);        pipe->kMap        = NULL;}
    if (pipe->variance)    {free_FFT_arr(pipe->variance);       pipe->variance    = NULL;}
    if (pipe->peakList)    {free_double_arr(pipe->peakList);    pipe->peakList    = NULL;}
    if (pipe->nuHist)      {free_hist_t(pipe->nuHist);          pipe->nuHist      = NULL;}
    if (pipe->multiscale)  {free_double_mat(pipe->multiscale);  pipe->multiscale  = NULL;}
    
    free(pipe); pipe = NULL;
  }
  return;
}

//----------------------------------------------------------------------
//-- Main function

void doMultiscale(cosmo_hm *chPar, peak_param *pkPar, error **err)
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
  FFT_arr *variance       = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  makeKernelForVariance(pkPar, variance);
  double_arr *peakList    = initialize_double_arr(pkPar->peakListLength);
  hist_t *nuHist          = initialize_hist_t(pkPar->N_nu);
  setHist_nu(pkPar, nuHist);
  double_mat *multiscale  = initialize_double_mat(pkPar->N_nu, pkPar->nbFilters);
  
  char name[STRING_LENGTH_MAX];
  int j;
  
  readCatOrMakeSimulAndOutput(chPar, pkPar, hSampArr, hMap, err);                                            forwardError(*err, __LINE__,);
  cleanOrMakeOrResample(chPar, pkPar, gSamp, gMap, mask, err);                                               forwardError(*err, __LINE__,);
  lensingCatalogueAndOutput(chPar, pkPar, hMap, gMap, k1Inter, err);                                         forwardError(*err, __LINE__,);
  makeMapsAndOutput(chPar, pkPar, gMap, gMap2, FFTSmoother, DCSmoother, kMap, err);                          forwardError(*err, __LINE__,);
  makeLocalVariance(pkPar, gMap, variance);
  mapToMultiscaleAndOutput_FFT(pkPar, gMap, FFTSmoother, kMap, variance, peakList, nuHist, multiscale, err); forwardError(*err, __LINE__,);
  mapToMultiscaleAndOutput_DC(pkPar, gMap, DCSmoother, kMap, variance, peakList, nuHist, multiscale, err);   forwardError(*err, __LINE__,);
  
  //-- Output
  sprintf(name, "%smultiscale_nbFilters%d", pkPar->prefix, pkPar->nbFilters);
  outputMultiscale(name, chPar, pkPar, multiscale, err); forwardError(*err, __LINE__,);
  
  //-- Free
  free_sampler_arr(hSampArr);
  free_interpolator_t(k1Inter);
  free_halo_map(hMap);
  free_sampler_t(gSamp);
  free_gal_map(gMap);
  free_mask_map(mask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_signal_map(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(nuHist);
  free_double_mat(multiscale);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doDataMatrix(cosmo_hm *chPar, peak_param *pkPar, int N, error **err)
{
  if (pkPar->doInHaloCat > 0)   printf("Found doInHaloCat = 1, forced to 0\n");
  pkPar->doInHaloCat = 0;
  if (pkPar->doInGalCat > 0)    printf("Found doInGalCat = 1, forced to 0\n");
  pkPar->doInGalCat = 0;
  if (pkPar->outHaloCat > 0)    printf("Found outHaloCat = 1, forced to 0\n");
  pkPar->outHaloCat = 0;
  if (pkPar->outGalCat > 0)     printf("Found outGalCat = 1, forced to 0\n");
  pkPar->outGalCat = 0;
  if (pkPar->outMaps > 0)       printf("Found outMaps = 1, forced to 0\n");
  pkPar->outMaps = 0;
  if (pkPar->outTruth > 0)      printf("Found outTruth = 1, forced to 0\n");
  pkPar->outTruth = 0;
  if (pkPar->outMask > 0)       printf("Found outMask = 1, forced to 0\n");
  pkPar->outMask = 0;
  if (pkPar->outPeakList > 0)   printf("Found outPeakList = 1, forced to 0\n");
  pkPar->outPeakList = 0;
  if (pkPar->outHist > 0)       printf("Found outHist = 1, forced to 0\n");
  pkPar->outHist = 0;
  if (pkPar->outMultiscale > 0) printf("Found outMultiscale = 1, forced to 0\n");
  pkPar->outMultiscale = 0;
  
  sampler_arr *hSampArr   = initialize_sampler_arr(pkPar->N_z_halo, pkPar->N_M+1);
  interpolator_t *k1Inter = initialize_interpolator_t(pkPar->N_z_gal+1);
  setMassSampAndK1Inter(chPar, pkPar, hSampArr, NULL, k1Inter, err);                                        forwardError(*err, __LINE__,); //-- lambda = NULL
  halo_map *hMap          = initialize_halo_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);   forwardError(*err, __LINE__,);
  sampler_t *gSamp        = initialize_sampler_t(pkPar->N_z_gal+1);
  setGalaxySampler(chPar, pkPar, gSamp, err);                                                               forwardError(*err, __LINE__,);
  gal_map *gMap           = initialize_gal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);    forwardError(*err, __LINE__,);
  mask_map *mask          = initializeMask(pkPar, err);                                                     forwardError(*err, __LINE__,);
  FFT_arr *FFTSmoother    = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  makeKernel(pkPar, FFTSmoother);
  FFT_arr *DCSmoother     = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  signal_map *kMap        = initialize_signal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err); forwardError(*err, __LINE__,);
  FFT_arr *variance       = initialize_FFT_arr(pkPar->smootherSize, pkPar->FFTSize);
  makeKernelForVariance(pkPar, variance);
  double_arr *peakList    = initialize_double_arr(pkPar->peakListLength);
  hist_t *nuHist          = initialize_hist_t(pkPar->N_nu);
  setHist_nu(pkPar, nuHist);
  double_mat *multiscale  = initialize_double_mat(pkPar->N_nu, pkPar->nbFilters);
  double_mat *dataMat     = initialize_double_mat(pkPar->N_nu * pkPar->nbFilters, N);
  gsl_matrix *invCov      = gsl_matrix_alloc(pkPar->N_nu * pkPar->nbFilters, pkPar->N_nu * pkPar->nbFilters);
  
  clock_t start = clock();
  int i;
  
  //-- Peak histogram loop
  printf("\n-- Making N = %d realizations\n", N);
  for (i=0; i<N; i++) {
    printf("Realization %4d: ", i+1); fflush(stdout);
    pkPar->realizationInd = i;
    massFctToMultiscale(chPar, pkPar, hSampArr, k1Inter, hMap, gSamp, gMap, mask, FFTSmoother, DCSmoother, kMap, 
			variance, peakList, nuHist, multiscale, 1, err); forwardError(*err, __LINE__,); //-- nbPatches = 1
    forwardResult(pkPar, multiscale, dataMat->matrix);
  }
  
  double duration = (double)(clock() - start) / CLOCKS_PER_SEC;
  printf("%d secs for %d realizations, %.2f secs per run\n", (int)duration, N, duration/N);
  printf("\n");
  
  //-- Output
  char name[STRING_LENGTH_MAX];
  sprintf(name, "%sdataMat_nbFilters%d_N%d", pkPar->prefix, pkPar->nbFilters, N);
  outputDataMat(name, chPar, pkPar, dataMat, err); forwardError(*err, __LINE__,);
  
  inverseCovarianceMatrix(dataMat, invCov);
  sprintf(name, "%sinvCov_nbFilters%d_N%d", pkPar->prefix, pkPar->nbFilters, N);
  outputInvCov(name, chPar, pkPar, invCov, err); forwardError(*err, __LINE__,);
  
  //-- Free
  free_sampler_arr(hSampArr);
  free_interpolator_t(k1Inter);
  free_halo_map(hMap);
  free_sampler_t(gSamp);
  free_gal_map(gMap);
  free_mask_map(mask);
  free_FFT_arr(FFTSmoother);
  free_FFT_arr(DCSmoother);
  free_signal_map(kMap);
  free_FFT_arr(variance);
  free_double_arr(peakList);
  free_hist_t(nuHist);
  free_double_mat(multiscale);
  free_double_mat(dataMat);
  gsl_matrix_free(invCov);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

