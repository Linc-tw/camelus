

  /***********************************************
   **  ABC.c					**
   **  Chieh-An Lin				**
   **  Version 2016.03.26			**
   **						**
   **  References:				**
   **  - Marin et al. (2011)			**
   **  - Weyant et al. (2013) - ApJ, 764, 116	**
   ***********************************************/


#include "ABC.h"


//----------------------------------------------------------------------
//-- Functions related to print

void valueStringForPrint(char buffer2[], param_t p, double value)
{
  //-- To change if the parameter space is enlarged
  
  switch (p) {
    case param_Omega_b:
      sprintf(buffer2, "%6.4f", value);  break;
    case param_Omega_m: case param_Omega_de: case param_n_s: case param_h_100: case param_sigma_8: 
      sprintf(buffer2, "%5.3f", value);  break;
    case param_w0_de: case param_w1_de: case param_beta_NFW:
      sprintf(buffer2, "% 5.3f", value); break;
    case param_c_0:
      sprintf(buffer2, "%5.2f", value);  break;
  }
  return;
}

void paramStringForOutput(char buffer2[], param_t p)
{
  //-- To change if the parameter space is enlarged
  
  switch (p) {
    case param_Omega_m:
      sprintf(buffer2, "%s", "   Omega_m"); break;
    case param_Omega_de: 
      sprintf(buffer2, "%s", "  Omega_de"); break;
    case param_Omega_b:
      sprintf(buffer2, "%s", "   Omega_b"); break;
    case param_n_s:
      sprintf(buffer2, "%s", "       n_s"); break;
    case param_h_100:
      sprintf(buffer2, "%s", "     h_100"); break;
    case param_sigma_8:
      sprintf(buffer2, "%s", "   sigma_8"); break;
    case param_w0_de:
      sprintf(buffer2, "%s", "     w0_de"); break;
    case param_w1_de:
      sprintf(buffer2, "%s", "     w1_de"); break;
    case param_c_0:
      sprintf(buffer2, "%s", "       c_0"); break;
    case param_beta_NFW:
      sprintf(buffer2, "%s", "  beta_NFW"); break;
  }
  return;
}

void valueStringForOutput(char buffer2[], param_t p, double value)
{
  //-- To change if the parameter space is enlarged
  
  switch (p) {
    case param_Omega_b:
      sprintf(buffer2, "%8.6f", value);
      break;
    case param_Omega_m: case param_Omega_de: case param_n_s: case param_h_100: case param_sigma_8:
      sprintf(buffer2, "%8.5f", value);
      break;
    case param_w0_de: case param_w1_de: case param_beta_NFW:
      sprintf(buffer2, "% 7.5f", value);
      break;
    case param_c_0:
      sprintf(buffer2, "%8.4f", value);
      break;
  }
  return;
}

void particleString(char buffer1[], char buffer2[], int f, param_t *doParam, double *part)
{
  int count = 0;
  int i;
  
  sprintf(buffer1, "(%s", STR_PARAM_T(doParam[0])); 
  for (i=1; i<f; i++) sprintf(buffer1, "%s, %s", buffer1, STR_PARAM_T(doParam[i]));
  
  valueStringForPrint(buffer2, doParam[0], part[0]);
  sprintf(buffer1, "%s) = (%s", buffer1, buffer2);
  for (i=1; i<f; i++) {
    valueStringForPrint(buffer2, doParam[i], part[i]);
    sprintf(buffer1, "%s, %s", buffer1, buffer2);
  }
  sprintf(buffer1, "%s)", buffer1);
  return;
}

void completeParticleString(char buffer1[], char buffer2[], double *value)
{
  int count = 0;
  int i;
  
  sprintf(buffer1, "(%s", STR_PARAM_T(0)); 
  for (i=1; i<NB_PARAM_T; i++) sprintf(buffer1, "%s, %s", buffer1, STR_PARAM_T(i));
  
  valueStringForPrint(buffer2, (param_t)0, value[0]);
  sprintf(buffer1, "%s) = (%s", buffer1, buffer2);
  for (i=1; i<NB_PARAM_T; i++) {
    valueStringForPrint(buffer2, (param_t)i, value[i]);
    sprintf(buffer1, "%s, %s", buffer1, buffer2);
  }
  sprintf(buffer1, "%s)", buffer1);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to iteration_t

iteration_t *initialize_iteration_t(int f, int *doParam, int Q, int MPISize, error **err)
{
  iteration_t *iter = (iteration_t*)malloc_err(sizeof(iteration_t), err); forwardError(*err, __LINE__, NULL);
  iter->f           = f;
  iter->doParam     = (param_t*)malloc_err(f * sizeof(param_t), err);     forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<f; i++) iter->doParam[i] = (param_t)doParam[i];
  iter->Q           = Q;
  iter->Q_MPI       = (int)ceil((double)Q / (double)MPISize);
  iter->nbAttempts  = 0;
  iter->epsilon     = ABC_TOLERANCE_MAX;
  iter->mean        = gsl_vector_alloc(f);
  iter->buffer      = gsl_vector_alloc(f);
  iter->cov         = gsl_matrix_alloc(f, f);
  iter->cov2        = gsl_matrix_alloc(f, f);
  iter->invCov      = gsl_matrix_alloc(f, f);
  iter->cholesky    = gsl_matrix_alloc(f, f);
  iter->perm        = gsl_permutation_alloc(f);
  iter->matrix      = initialize_double_mat(f + 2, iter->Q_MPI * MPISize);
  return iter;
}

void free_iteration_t(iteration_t *iter)
{
  if (iter->doParam) {free(iter->doParam); iter->doParam = NULL;}
  gsl_vector_free(iter->mean);
  gsl_vector_free(iter->buffer);
  gsl_matrix_free(iter->cov);
  gsl_matrix_free(iter->cov2);
  gsl_matrix_free(iter->invCov);
  gsl_matrix_free(iter->cholesky);
  gsl_permutation_free(iter->perm);
  free_double_mat(iter->matrix);
  free(iter);
  return;
}

void print_gsl_matrix(gsl_matrix *mat)
{
  int f = mat->size1;
  int i, j;
  
  printf("  [[% 9.5f", mat->data[0+0*f]);
  for (i=1; i<f; i++) printf(", % 9.5f", mat->data[i+0*f]);
  printf("]");
  
  for (j=1; j<f; j++) {
    printf(",\n   [% 9.5f", mat->data[0+j*f]);
    for (i=1; i<f; i++) printf(", % 9.5f", mat->data[i+j*f]);
    printf("]");
  }
  
  printf("]\n");
  return;
}

void print_iteration_t(iteration_t *iter)
{
  printf("# iteration_t (first 20 elements)\n");
  printf("nbAttempts = %d\n", iter->nbAttempts);
  printf("epsilon    = %f\n", iter->epsilon);
  printf("mean       = ");  printDoubleArray(iter->mean->data, iter->mean->size, 3); printf("\n");
  printf("cov        =\n"); print_gsl_matrix(iter->cov);
  printf("invCov     =\n"); print_gsl_matrix(iter->invCov);
  
  int f2 = iter->f + 2;
  int L  = MIN(iter->Q, 20);
  
  double *begin = iter->matrix->matrix;
  double *end   = begin + f2 * L;
  char buffer1[STRING_LENGTH_MAX], buffer2[16];
  double *part;
  for (part=begin; part<end; part+=f2) {
    particleString(buffer1, buffer2, iter->f, iter->doParam, part);
    printf("%s, delta = %9.5f, weight = %.3f\n", buffer1, part[f2-2], part[f2-1]);
  }
  return;
}

void output_gsl_matrix(FILE *file, gsl_matrix *mat)
{
  int f = mat->size1;
  int i, j;
  
  fprintf(file, "#   [[% 9.5f", mat->data[0+0*f]);
  for (i=1; i<f; i++) fprintf(file, ", % 9.5f", mat->data[i+0*f]);
  fprintf(file, "]");
  
  for (j=1; j<f; j++) {
    fprintf(file, ",\n#    [% 9.5f", mat->data[0+j*f]);
    for (i=1; i<f; i++) fprintf(file, ", % 9.5f", mat->data[i+j*f]);
    fprintf(file, "]");
  }
  
  fprintf(file, "]\n");
  return;
}

void output_iteration_t(FILE *file, iteration_t *iter)
{
  int Q = iter->Q;
  fprintf(file, "# iteration_t\n");
  fprintf(file, "# f                  = %d\n", iter->f);
  fprintf(file, "# Q                  = %d\n", Q);
  fprintf(file, "# Number of attempts = %d\n", iter->nbAttempts);
  fprintf(file, "# Success rate       = %.5f\n", Q / (double)iter->nbAttempts);
  fprintf(file, "# epsilon            = %.5f\n", iter->epsilon);
  fprintf(file, "# cov                =\n");
  output_gsl_matrix(file, iter->cov);
  fprintf(file, "# invCov             =\n");
  output_gsl_matrix(file, iter->invCov);
  fprintf(file, "#\n");
  
  int f            = iter->f;
  param_t *doParam = iter->doParam;
  char buffer2[16];
  int i;
  
  fprintf(file, "#");
  for (i=0; i<f; i++) {
    paramStringForOutput(buffer2, doParam[i]);
    fprintf(file, "%s", buffer2);
  }
  fprintf(file, "      delta        weight\n");
  
  int f2        = f + 2;
  double *begin = iter->matrix->matrix;
  double *end   = begin + f2 * Q;
  double *part;
  
  for (part=begin; part<end; part+=f2) {
    fprintf(file, " ");
    for (i=0; i<f; i++) {
      valueStringForOutput(buffer2, doParam[i], part[i]);
      fprintf(file, "  %s", buffer2);
    }
    fprintf(file, "  %10.5e",   part[f]);
    fprintf(file, "  %10.5e\n", part[f+1]);
  }
  return;
}

void updateMean_iteration_t(iteration_t *iter)
{
  int f  = iter->f;
  int f2 = iter->f + 2;
  int Q  = iter->Q;
  double *weight = iter->matrix->matrix + f2 - 1;
  double *param;
  int i;
  
  for (i=0; i<f; i++) {
    param = iter->matrix->matrix + i;
    iter->mean->data[i] = gsl_stats_wmean(weight, f2, param, f2, Q);
  }
  return;
}

void updateCovariance_iteration_t(iteration_t *iter)
{
  int f2 = iter->f + 2;
  int Q  = iter->Q;
  double w_tot    = 0.0;
  double w_sq_tot = 0.0;
  double *begin   = iter->matrix->matrix;
  double *end     = begin + f2 * Q;
  double *part;
  
  for (part=begin; part<end; part+=f2) {
    w_tot    += part[f2-1];
    w_sq_tot += pow(part[f2-1], 2);
  }
  
  int f = iter->f;
  double *meanArr = iter->mean->data;
  double *covMat  = iter->cov->data;
  
  double mean_i, mean_j, sum;
  int i, j;
  
  for (j=0; j<f; j++) {
    for (i=0; i<f; i++) {
      sum = 0.0;
      if (i < j) {
	covMat[i+j*f] = covMat[j+i*f];
	continue;
      }
      mean_i = meanArr[i];
      mean_j = meanArr[j];
      for (part=begin; part<end; part+=f2) sum += part[f2-1] * (part[i] - mean_i) * (part[j] - mean_j);
      sum *= w_tot / (pow(w_tot, 2) - w_sq_tot);
      covMat[i+j*f] = sum;
    }
  }
  
  //-- Make a copy, because the invertion will destroy cov
  gsl_matrix_memcpy(iter->cov2, iter->cov);
  
  //-- Make LU decomposition and invert
  int s;
  gsl_linalg_LU_decomp(iter->cov2, iter->perm, &s);
  gsl_linalg_LU_invert(iter->cov2, iter->perm, iter->invCov);
  
  //-- Debias the C_{ij}^{-1}
  //-- C^-1_unbiased = (Q - f - 2) / (Q - 1) * C^-1_biased
  double factor = (Q - f - 2) / (double)(Q - 1);
  gsl_matrix_scale(iter->invCov, factor);
  return;
}

void updateCholesky_iteration_t(iteration_t *iter)
{
  gsl_matrix_memcpy(iter->cholesky, iter->cov);
  gsl_linalg_cholesky_decomp(iter->cholesky);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to PMC_ABC_t

PMC_ABC_t *initialize_PMC_ABC_t(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  PMC_ABC_t *ABC = (PMC_ABC_t*)malloc_err(sizeof(PMC_ABC_t), err);                                forwardError(*err, __LINE__, NULL);
  ABC->f         = peak->ABC_f;
  ABC->Q         = peak->ABC_Q;
  ABC->r_stop    = peak->ABC_r_stop;
  int i;
  ABC->summ      = -1;
  for (i=0; i<NB_SUMM_T; i++) if (strcmp(peak->ABC_summ, STR_SUMM_T((summ_t)i))==0) ABC->summ = (summ_t)i;
  testErrorVA((int)ABC->summ==-1, conf_undef, "Parametrization '%s' not found in type %s, cannot assign value of type %s", \
		*err, __LINE__, peak->ABC_summ, "STR_SUMM_T", "summ_t");                          forwardError(*err, __LINE__, NULL);
  
  ABC->Q_MPI     = (int)ceil((double)ABC->Q / (double)peak->MPISize);
  
  ABC->t         = -1; //-- To be increment later
  ABC->oldIter   = initialize_iteration_t(ABC->f, peak->ABC_doParam, ABC->Q, peak->MPISize, err); forwardError(*err, __LINE__, NULL);
  ABC->newIter   = initialize_iteration_t(ABC->f, peak->ABC_doParam, ABC->Q, peak->MPISize, err); forwardError(*err, __LINE__, NULL);
  ABC->deltaList = initialize_double_arr(ABC->Q);
  
  ABC->priorFct  = prior_rectangle;
  ABC->limit     = initialize_double_mat(2, ABC->f);                                              forwardError(*err, __LINE__, NULL);
  ABC->flatness  = garanteeFlatness(ABC->f, ABC->newIter->doParam);
  fillLimitAndValue(cmhm, ABC);
  
  int N_bin      = (peak->doNonlinear == 0) ? peak->N_nu : (peak->doSmoothing == 4) ? peak->N_kappa : MAX(peak->N_nu, peak->N_kappa);
  ABC->x_obs     = initialize_double_arr(N_bin * peak->nbFilters);
  ABC->x_mod     = initialize_double_arr(N_bin * peak->nbFilters);
  ABC->summFct   = summary_multiscale;
  
  //-- For Paper III, this will be overwritten by another function
  if (ABC->summ == summ_gauss)       ABC->distFct = dist_gauss;
  else if (ABC->summ == summ_star)   ABC->distFct = dist_star;
  else if (ABC->summ == summ_tanh) ;
  else if (ABC->summ == summ_mrlens) ;
  else if (ABC->summ == summ_FDR02) ;
  else if (ABC->summ == summ_FDR05) ;
  else {*err = addError(peak_unknown, "Unknown summary type", *err, __LINE__);                  forwardError(*err, __LINE__, NULL);}
  
  int length     = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  
  //-- Camelus pipeline
  ABC->sampArr     = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  ABC->hMap        = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__, NULL);
  ABC->galSamp     = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, ABC->galSamp, err);                                              forwardError(*err, __LINE__, NULL);
  ABC->gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__, NULL);
  ABC->CCDMask     = NULL; //-- Initialize seperately //-- initializeMask(peak, err);                forwardError(*err, __LINE__, NULL); 
  ABC->FFTSmoother = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernel(peak, ABC->FFTSmoother);
  ABC->DCSmoother  = initialize_FFT_arr(MAX(peak->DC_nbFilters, 1), peak->FFTSize);
  ABC->kMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__, NULL);
  ABC->variance    = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->FFT_nbFilters) makeKernelForVariance(peak, ABC->variance);
  ABC->peakList    = initialize_double_arr(length);
  ABC->nuHist      = initialize_hist_t(peak->N_nu);
  setHist_nu(peak, ABC->nuHist);
  ABC->kappaHist   = initialize_hist_t(peak->N_kappa);
  setHist_kappa(peak, ABC->kappaHist);
  ABC->multiscale  = initialize_double_mat(N_bin, peak->nbFilters);
  
  if (peak->MPIInd == 0) printf("ABC initialization done\n");
  return ABC;
}

void free_PMC_ABC_t(PMC_ABC_t *ABC)
{
  if (ABC->oldIter)     {free_iteration_t(ABC->oldIter);   ABC->oldIter     = NULL;}
  if (ABC->newIter)     {free_iteration_t(ABC->newIter);   ABC->newIter     = NULL;}
  if (ABC->deltaList)   {free_double_arr(ABC->deltaList);  ABC->deltaList   = NULL;}
  
  if (ABC->limit)       {free_double_mat(ABC->limit);      ABC->limit       = NULL;}
  
  if (ABC->x_obs)       {free_double_arr(ABC->x_obs);      ABC->x_obs       = NULL;}
  if (ABC->x_mod)       {free_double_arr(ABC->x_mod);      ABC->x_mod       = NULL;}
  
  if (ABC->sampArr)     {free_sampler_arr(ABC->sampArr);   ABC->sampArr     = NULL;}
  if (ABC->hMap)        {free_halo_map(ABC->hMap);         ABC->hMap        = NULL;}
  if (ABC->galSamp)     {free_sampler_t(ABC->galSamp);     ABC->galSamp     = NULL;}
  if (ABC->gMap)        {free_gal_map(ABC->gMap);          ABC->gMap        = NULL;}
  if (ABC->CCDMask)     {free_short_mat(ABC->CCDMask);     ABC->CCDMask     = NULL;}
  if (ABC->FFTSmoother) {free_FFT_arr(ABC->FFTSmoother);   ABC->FFTSmoother = NULL;}
  if (ABC->DCSmoother)  {free_FFT_arr(ABC->DCSmoother);    ABC->DCSmoother  = NULL;}
  if (ABC->kMap)        {free_map_t(ABC->kMap);            ABC->kMap        = NULL;}
  if (ABC->variance)    {free_FFT_arr(ABC->variance);      ABC->variance    = NULL;}
  if (ABC->peakList)    {free_double_arr(ABC->peakList);   ABC->peakList    = NULL;}
  if (ABC->nuHist)      {free_hist_t(ABC->nuHist);         ABC->nuHist      = NULL;}
  if (ABC->kappaHist)   {free_hist_t(ABC->kappaHist);      ABC->kappaHist   = NULL;}
  if (ABC->multiscale)  {free_double_mat(ABC->multiscale); ABC->multiscale  = NULL;}
  
  free(ABC); ABC = NULL;
  return;
}

void print_PMC_ABC_t(peak_param *peak, PMC_ABC_t *ABC)
{
  int i;
  printf("\nABC parameters\n");
  printf("Dimension of parameter   = %d\n", ABC->f);
  printf("Parameter space          = %s", STR_PARAM_T(ABC->newIter->doParam[0])); for (i=1; i<ABC->f; i++) printf(", %s", STR_PARAM_T(ABC->newIter->doParam[i])); printf("\n");
  printf("Number of particles      = %d = %d procs * %d - %d\n", ABC->Q, peak->MPISize, ABC->Q_MPI, peak->MPISize * ABC->Q_MPI - ABC->Q);
  printf("Shutoff success rate     = %.3f\n", ABC->r_stop);
  printf("Summary type             = %s\n", peak->ABC_summ);
  printf("Dimension of data vector = %d\n", ABC->x_obs->length);
  printf("------------------------------------------------------------------------\n");
  return;
}

int garanteeFlatness(int f, param_t *doParam)
{
  //-- If Omega_m is free parameter, but not Omega_de, the Universe is flat and return +1.
  //-- If Omega_de is free parameter, but not Omega_m, the Universe is flat and return -1.
  //-- If both are free parameters, the Universe is not flat and return 0.
  //-- If none of both is free parameter, return 0.
  
  int flatness = 0;
  int i;
  for (i=0; i<f; i++) {
    if (doParam[i] == param_Omega_m) {
      flatness += 1;
      break;
    }
  }
  for (i=0; i<f; i++) {
    if (doParam[i] == param_Omega_de) {
      flatness -= 1;
      break;
    }
  }
  return flatness;
}

void fillLimitAndValue(cosmo_hm *cmhm, PMC_ABC_t *ABC)
{
  //-- To change if the parameter space is enlarged
  
  double buffer[2*NB_PARAM_T];
  double *value  = ABC->value;
  int i;
  
  for (i=0; i<NB_PARAM_T; i++) {
    switch ((param_t)i) {
      case param_Omega_m:
	buffer[0+2*i] = OMEGA_M_MIN;
	buffer[1+2*i] = OMEGA_M_MAX;
	value[i]      = cmhm->cosmo->Omega_m;
	break;
      case param_Omega_de:
	buffer[0+2*i] = OMEGA_DE_MIN;
	buffer[1+2*i] = OMEGA_DE_MAX;
	value[i]      = cmhm->cosmo->Omega_de;
	break;
      case param_Omega_b:
	buffer[0+2*i] = OMEGA_B_MIN;
	buffer[1+2*i] = OMEGA_B_MAX;
	value[i]      = cmhm->cosmo->Omega_b;
	break;
      case param_n_s:
	buffer[0+2*i] = N_S_MIN;
	buffer[1+2*i] = N_S_MAX;
	value[i]      = cmhm->cosmo->n_spec;
	break;
      case param_h_100:
	buffer[0+2*i] = H_100_MIN;
	buffer[1+2*i] = H_100_MAX;
	value[i]      = cmhm->cosmo->h_100;
	break;
      case param_sigma_8:
	buffer[0+2*i] = SIGMA_8_MIN;
	buffer[1+2*i] = SIGMA_8_MAX;
	value[i]      = cmhm->cosmo->normalization;
	break;
      case param_w0_de:
	buffer[0+2*i] = W0_DE_MIN;
	buffer[1+2*i] = W0_DE_MAX;
	value[i]      = cmhm->cosmo->w0_de;
	break;
      case param_w1_de:
	buffer[0+2*i] = W1_DE_MIN;
	buffer[1+2*i] = W1_DE_MAX;
	value[i]      = cmhm->cosmo->w1_de;
	break;
      case param_c_0:
	buffer[0+2*i] = C_0_MIN;
	buffer[1+2*i] = C_0_MAX;
	value[i]      = cmhm->c0;
	break;
      case param_beta_NFW:
	buffer[0+2*i] = BETA_NFW_MIN;
	buffer[1+2*i] = BETA_NFW_MAX;
	value[i]      = cmhm->beta_NFW;
	break;
    }
  }
  
  double *matrix = ABC->limit->matrix;
  param_t *doParam = ABC->newIter->doParam;
  int j;
  
  for (i=0; i<ABC->f; i++) {
    j = (int)doParam[i];
    matrix[0+2*i] = buffer[0+2*j];
    matrix[1+2*i] = buffer[1+2*j];
  }
  return;
}

void fillObservation(char name[], peak_param *peak, PMC_ABC_t *ABC, error **err)
{
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  int c = fgetc(file);
  int count = 0;
  
  char buffer[STRING_LENGTH_MAX];
  char *buffer1;
  int buffer2;
  double buffer3[45];
  
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf", &buffer3[count]);
      count++;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  int i;
  
  if (ABC->summ == summ_gauss) {
    for (i=0; i<18; i++) ABC->x_obs->array[i] = buffer3[i];
  }
  else if (ABC->summ == summ_star) {
    for (i=0; i<27; i++) ABC->x_obs->array[i] = buffer3[i+18];
  }
  
  if (peak->MPIInd == 0) printf("\"%s\" read\n", name);
  return;
}

void fillIteration(char name[], peak_param *peak, PMC_ABC_t *ABC, int t, error **err)
{
  ABC->t = t;
  if (t < 0) return;
  
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  int f2 = ABC->f + 2;
  double *part = ABC->newIter->matrix->matrix;
  int i;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
    else {
      testErrorRet(count>=ABC->Q, peak_overflow, "Too many particles", *err, __LINE__,);
      ungetc(c, file);
      for (i=0; i<f2; i++) buffer2 = fscanf(file, "%lf", &part[i]);
      buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
      count++;
      part += f2;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  testErrorRet(count<ABC->Q, peak_match, "Number of particles not matched", *err, __LINE__,);
  
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
  
  if (peak->MPIInd == 0) {
    printf("\"%s\" read\n", name);
    printf("%d particles generated\n", ABC->Q);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to accept-reject algorithm

void generateParam(peak_param *peak, PMC_ABC_t *ABC, double *newPart, prior_fct *prior)
{
  int f  = ABC->f;
  int f2 = f + 2;
  double *oldBegin     = ABC->oldIter->matrix->matrix;
  gsl_rng *generator   = peak->generator;
  gsl_vector *buffer   = ABC->oldIter->buffer;
  gsl_matrix *cholesky = ABC->oldIter->cholesky;
  
  double *target;
  double q;
  int i;
  
  do {
    target = oldBegin;
    q = gsl_ran_flat(generator, 0.0, 1.0);
    while (q > target[f2-1]) {
      q      -= target[f2-1];
      target += f2;
    }
    
    //-- Compute pdf from multivariate Gaussian pdf
    for (i=0; i<f; i++) buffer->data[i] = gsl_ran_ugaussian(generator);
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, cholesky, buffer);
    for (i=0; i<f; i++) newPart[i] = target[i] + buffer->data[i];
  } while (!prior(ABC->limit, newPart));
  
  return;
}

cosmo_hm *initialize_cosmo_hm_ABC(cosmo_hm *oldCmhm, PMC_ABC_t *ABC, double *newPart, error **err)
{
  //-- If the number of parameters changes this need to be changed everywhere.
  
  //-- OMEGAM, OMEGADE, W0_DE, W1_DE, *W_POLY_DE, N_POLY_DE,
  //-- H100, OMEGAB, OMEGANUMASS, NEFFNUMASS, NORM, NSPEC,
  //-- NZBIN, *Nnz, *nofz, *par_nz, zmin, zmax,
  //-- NONLINEAR, TRANSFER, GROWTH, DEPARAM, normmode,
  //-- C0, ALPHANFW, BETANFW, MASSFCT, HALO_BIAS,
  //-- M_min, M1, M0, sigma_log_M, alpha,
  //-- Mstar0, beta, delta, gamma, B_cut, B_sat, 
  //-- beta_cut, beta_sat, Mstellar_min, Mstellar_max, eta,
  //-- fcen1, fcen2,
  //-- HOD, pi_max, err
  
  param_t *doParam = ABC->newIter->doParam;
  double *value    = ABC->value;
  int i;
  
  for (i=0; i<ABC->f; i++) value[doParam[i]] = newPart[i];
  if      (ABC->flatness == 1)  value[param_Omega_de] = 1.0 - value[param_Omega_m];
  else if (ABC->flatness == -1) value[param_Omega_m]  = 1.0 - value[param_Omega_de];
  
  cosmo *oldCm       = oldCmhm->cosmo;
  redshift_t * oldRs = oldCmhm->redshift;
  cosmo_hm *newCmhm  = init_parameters_hm(  value[param_Omega_m], value[param_Omega_de], value[param_w0_de],
    value[param_w1_de],    oldCm->w_poly_de,  oldCm->N_poly_de, value[param_h_100],   value[param_Omega_b],
	  oldCm->Omega_nu_mass,  oldCm->Neff_nu_mass,   value[param_sigma_8],  value[param_n_s],
      oldRs->Nzbin, oldRs->Nnz,   oldRs->nofz, oldRs->photz, oldRs->par_nz,   oldCmhm->zmin,  oldCmhm->zmax,
      oldCm->nonlinear,     oldCm->transfer,     oldCm->growth,  oldCm->de_param, oldCmhm->cosmo->normmode,  
      value[param_c_0],     oldCmhm->alpha_NFW,  value[param_beta_NFW], oldCmhm->massfct, oldCmhm->halo_bias,   
      oldCmhm->log10M_min,  oldCmhm->log10M1,    oldCmhm->log10M0,   oldCmhm->sigma_log_M,  oldCmhm->alpha, 
      oldCmhm->log10Mstar0, oldCmhm->beta,  		 oldCmhm->delta,  oldCmhm->gamma,  oldCmhm->B_cut, oldCmhm->B_sat, 
      oldCmhm->beta_cut,    oldCmhm->beta_sat,   oldCmhm->log10Mstar_min, oldCmhm->log10Mstar_min, oldCmhm->eta,
      oldCmhm->fcen1,       oldCmhm->fcen2, oldCmhm->hod,         oldCmhm->pi_max,    err);
  forwardError(*err, __LINE__, NULL);
  return newCmhm;
}

void generateModel(cosmo_hm *cmhm, peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err)
{
  cosmo_hm *cmhmABC = initialize_cosmo_hm_ABC(cmhm, ABC, newPart, err);            forwardError(*err, __LINE__,);
  setMassSamplers(cmhmABC, peak, ABC->sampArr, err);                               forwardError(*err, __LINE__,);
  if (peak->doRandGalPos == 0) updateCosmo_gal_map(cmhmABC, peak, ABC->gMap, err); forwardError(*err, __LINE__,);
  
  multiscaleFromMassFct(cmhmABC, peak, ABC->sampArr, ABC->hMap, ABC->galSamp, ABC->gMap, ABC->CCDMask, 
			ABC->FFTSmoother, ABC->DCSmoother, ABC->kMap, ABC->variance, ABC->peakList, ABC->nuHist, ABC->kappaHist, 
			ABC->multiscale, err); forwardError(*err, __LINE__,);
  
  ABC->summFct(ABC->peakList, NULL, ABC->multiscale, ABC->x_mod->array);
  newPart[ABC->f] = ABC->distFct(ABC->x_obs->array, ABC->x_mod->array); 
  free_parameters_hm(&cmhmABC);
  return;
}

void acceptParticleFromPrior(cosmo_hm *cmhm, peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err)
{
  int f = ABC->f;
  int f2 = f + 2;
  int reject = 1;
  param_t *doParam = ABC->newIter->doParam;
  char buffer1[STRING_LENGTH_MAX], buffer2[16];
  
  do {
    do priorGenerator(peak->generator, ABC->limit, newPart);
    while (!ABC->priorFct(ABC->limit, newPart));
    
    particleString(buffer1, buffer2, f, doParam, newPart);
    if (peak->printMode == 2) printf("%s, ", buffer1);
    
    generateModel(cmhm, peak, ABC, newPart, err);
    
    if (isError(*err)) {
      if      (peak->printMode == 2) printf("get error and resample\n");
      else if (peak->printMode == 3) printf("Proc %2d: %s, get error and resample\n", peak->MPIInd, buffer1);
      purgeError(err);
    }
    
    else {
      ABC->newIter->nbAttempts += 1;
      reject = 0;
      if      (peak->printMode == 2) printf("delta = %7.3f, accepted\n", newPart[f]);
      else if (peak->printMode == 3) printf("Proc %2d: %s, %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, accepted\n", 
					    peak->MPIInd, buffer1, ABC->hMap->total, ABC->gMap->total, ABC->gMap->kappa_mean, ABC->gMap->fillingThreshold, newPart[f]);
    }
  }
  while (reject);
  return;
}

void acceptParticle(cosmo_hm *cmhm, peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err)
{
  int f = ABC->f;
  int reject = 1;
  param_t *doParam = ABC->newIter->doParam;
  char buffer1[STRING_LENGTH_MAX], buffer2[16];
  
  do {
    generateParam(peak, ABC, newPart, ABC->priorFct);
    particleString(buffer1, buffer2, f, doParam, newPart);
    if (peak->printMode == 2) printf("%s, ", buffer1);
    
    generateModel(cmhm, peak, ABC, newPart, err);
    
    if (isError(*err)) {
      if      (peak->printMode == 2) printf("get error and resample\n");
      else if (peak->printMode == 3) printf("Proc %2d: %s, get error and resample\n", peak->MPIInd, buffer1);
      purgeError(err);
    }
    
    else {
      ABC->newIter->nbAttempts += 1;
      
      if (newPart[f] <= ABC->newIter->epsilon) {
	reject = 0;
	if      (peak->printMode == 2) printf("delta = %6.3f, accepted\n", newPart[f]);
	else if (peak->printMode == 3) printf("Proc %2d: %s, %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, accepted\n", 
					      peak->MPIInd, buffer1, ABC->hMap->total, ABC->gMap->total, ABC->gMap->kappa_mean, ABC->gMap->fillingThreshold, newPart[f]);
      }
      else {
	if      (peak->printMode == 2) printf("delta = %6.3f, rejected\n", newPart[f]);
	else if (peak->printMode == 3) printf("Proc %2d: %s, %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, rejected\n", 
					      peak->MPIInd, buffer1, ABC->hMap->total, ABC->gMap->total, ABC->gMap->kappa_mean, ABC->gMap->fillingThreshold, newPart[f]);
      }
    }
  }
  while (reject);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to iteration

void swapOldAndNew(PMC_ABC_t *ABC)
{
  iteration_t *buffer = ABC->oldIter;
  ABC->oldIter = ABC->newIter;
  ABC->newIter = buffer;
  return;
}

double argOfExp(double *oldPart, double *newPart, gsl_matrix *invCov, gsl_vector *Delta_param, gsl_vector *intermediate)
{
  double value;
  int i;
  for (i=0; i<Delta_param->size; i++) Delta_param->data[i] = newPart[i] - oldPart[i];
  gsl_blas_dsymv(CblasUpper, 1.0, invCov, Delta_param, 0.0, intermediate); //-- intermediate = invCov * Delta_param
  gsl_blas_ddot(Delta_param, intermediate, &value);
  value *= -0.5;
  return value;
}

void setWeights(PMC_ABC_t *ABC)
{
  int f2 = ABC->f + 2;
  int Q  = ABC->Q;
  double weight    = 1.0 / (double)Q;
  double *oldBegin = ABC->oldIter->matrix->matrix;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *oldEnd   = oldBegin + f2 * Q;
  double *newEnd   = newBegin + f2 * Q;
  
  double *newPart;
  
  if (ABC->t == 0) {
    for (newPart=newBegin; newPart<newEnd; newPart+=f2) newPart[f2-1] = weight;
    return;
  }
  
  gsl_matrix *invCov       = ABC->oldIter->invCov;
  gsl_vector *Delta_param  = ABC->oldIter->buffer;
  gsl_vector *intermediate = ABC->newIter->buffer;
  double sum = 0.0;
  
  double *oldPart;
  double arg;
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    weight = 0.0;
    for (oldPart=oldBegin; oldPart<oldEnd; oldPart+=f2) {
      arg = argOfExp(oldPart, newPart, invCov, Delta_param, intermediate);
      weight += oldPart[f2-1] * exp(arg);
    }
    weight = 1.0 / weight;
    newPart[f2-1] = weight;
    sum += weight;
  }
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) newPart[f2-1] /= sum;
  return;
}

void setEpsilon(PMC_ABC_t *ABC)
{
  //-- After swapping
  if (ABC->t == 0) {
    ABC->newIter->epsilon = ABC_TOLERANCE_MAX;
    return;
  }
  
  int f2 = ABC->f + 2;
  int Q  = ABC->Q;
  double *oldBegin = ABC->oldIter->matrix->matrix + f2 - 2;
  double *deltaArr = ABC->deltaList->array;
  
  double *oldPart;
  int i;
  
  for (i=0, oldPart=oldBegin; i<Q; i++, oldPart+=f2) deltaArr[i] = oldPart[0];
  gsl_sort(deltaArr, 1, Q);
  ABC->newIter->epsilon = gsl_stats_median_from_sorted_data(deltaArr, 1, Q);
  return;
}

void runIteration(cosmo_hm *cmhm, peak_param *peak, PMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  swapOldAndNew(ABC);
  ABC->t++;
  ABC->newIter->nbAttempts = 0;
  setEpsilon(ABC);
  
  if (peak->MPIInd == 0) {
    printf("-- Begin iteration %d\n\n", ABC->t);
    if (ABC->t == 0) printf("Tolerance = DBL_MAX\n\n");
    else             printf("Tolerance = %.5f\n\n", ABC->newIter->epsilon);
  }
#ifndef __releaseMenu__
  MPI_Barrier(MPI_COMM_WORLD); //-- Wait for print
#endif
  
  int f2    = ABC->f + 2;
  int Q     = ABC->Q;
  int Q_MPI = ABC->Q_MPI;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *newEnd   = newBegin + f2 * Q_MPI;
  double *newPart;
  
  //-- Size ajustment for MPI
  if (peak->MPIInd == peak->MPISize - 1) newEnd -= f2 * (peak->MPISize * Q_MPI - Q);
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    if (ABC->t == 0) acceptParticleFromPrior(cmhm, peak, ABC, newPart, err);
    else             acceptParticle(cmhm, peak, ABC, newPart, err);
    forwardError(*err, __LINE__,);
  }
  printf("Proc %2d finished\n", peak->MPIInd);
  
#ifndef __releaseMenu__
  MPI_Allgather(newBegin, f2 * Q_MPI, MPI_DOUBLE, newBegin, f2 * Q_MPI, MPI_DOUBLE, MPI_COMM_WORLD);        //-- Communicate new particles
  MPI_Allreduce(&ABC->newIter->nbAttempts, &ABC->newIter->nbAttempts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //-- Communicate the number of attempts
#endif
  
  setWeights(ABC);
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
  
  if (peak->MPIInd == 0) {
    printf("\nSuccess rate = %.5f\n", Q / (double)ABC->newIter->nbAttempts);
    printf("\n-- End iteration %d, ", ABC->t); routineTime(start, clock());
    printf("------------------------------------------------------------------------\n");
  }
  return;
}

void outputIteration(char name[], PMC_ABC_t *ABC)
{
  FILE *file = fopen(name, "w");
  output_iteration_t(file, ABC->newIter);
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to subset

void setRealQ(PMC_ABC_t *ABC, int Q_real)
{
  //-- ABC->Q used in setEpsilon
  ABC->Q_MPI = Q_real; //-- Used in print_PMC_ABC_t_subset
  ABC->oldIter->Q     = Q_real;
  ABC->oldIter->Q_MPI = Q_real;
  ABC->newIter->Q     = Q_real; //-- Need to be peak->ABC_Q in readIteration, and Q_real in output_iteration_t
  ABC->newIter->Q_MPI = Q_real;
  return;
}

void runIterationWithoutUpdate(cosmo_hm *cmhm, peak_param *peak, PMC_ABC_t *ABC, error **err)
{
  //-- Assume that MPISize = 1
  
  clock_t start = clock();
  swapOldAndNew(ABC);
  ABC->t++;
  ABC->newIter->nbAttempts = 0;
  setEpsilon(ABC);
  
  printf("-- Begin iteration %d\n\n", ABC->t);
  if (ABC->t == 0) printf("Tolerance = DBL_MAX\n\n");
  else             printf("Tolerance = %.5f\n\n", ABC->newIter->epsilon);
  
  int f2    = ABC->f + 2;
  int Q_MPI = ABC->Q_MPI;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *newEnd   = newBegin + f2 * Q_MPI;
  double *newPart;
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    if (ABC->t == 0) acceptParticleFromPrior(cmhm, peak, ABC, newPart, err);
    else             acceptParticle(cmhm, peak, ABC, newPart, err);
    forwardError(*err, __LINE__,);
  }
  
  printf("\nSuccess rate = %.5f\n", Q_MPI / (double)ABC->newIter->nbAttempts);
  printf("\n-- End iteration %d, ", ABC->t); routineTime(start, clock());
  printf("------------------------------------------------------------------------\n");
  return;
}

void readSubset(char name[], iteration_t *iter, error **err)
{
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count;
  
  buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
  buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
  buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
  buffer2 = fscanf(file, "# Number of attempts = %d\n", &count); //-- Read the 4th line which contains nbAttempts
  iter->nbAttempts += count;
  count = 0;
  
  int f2       = iter->f + 2;
  double *part = iter->matrix->matrix + f2 * iter->Q_MPI;
  int i;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
    else {
      testErrorRet(count>=iter->Q, peak_overflow, "Too many particles", *err, __LINE__,);
      ungetc(c, file);
      for (i=0; i<f2; i++) buffer2 = fscanf(file, "%lf", &part[i]);
      buffer1 = fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
      count++;
      part += f2;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  iter->Q_MPI += count;
  printf("\"%s\" read\n", name);
  return;
}

void print_PMC_ABC_t_subset(peak_param *peak, PMC_ABC_t *ABC)
{
  int i;
  printf("\nABC parameters\n");
  printf("Dimension of parameter   = %d\n", ABC->f);
  printf("Parameter space          = %s", STR_PARAM_T(ABC->newIter->doParam[0])); for (i=1; i<ABC->f; i++) printf(", %s", STR_PARAM_T(ABC->newIter->doParam[i])); printf("\n");
  printf("Number of particles      = %d = %d subsets * %d\n", peak->ABC_Q, peak->ABC_Q/ABC->Q_MPI, ABC->Q_MPI);
  printf("Shutoff success rate     = %.3f\n", ABC->r_stop);
  printf("Summary type             = %s\n", peak->ABC_summ);
  printf("Dimension of data vector = %d\n", ABC->x_obs->length);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to prior

void priorGenerator(gsl_rng *generator, double_mat *limit, double *part)
{
  double *matrix = limit->matrix;
  int i;
  for (i=0; i<limit->N2; i++) part[i] = gsl_ran_flat(generator, matrix[0+2*i], matrix[1+2*i]);
  return;
}

int prior_rectangle(double_mat *limit, double *part)
{
  double *matrix = limit->matrix;
  int i;
  for (i=0; i<limit->N2; i++) {
    if (part[i] < matrix[0+2*i] || part[i] > matrix[1+2*i]) return 0;
  }
  return 1;
}

int prior_pentagon(double_mat *limit, double *part)
{
  double *matrix = limit->matrix;
  int i;
  for (i=0; i<limit->N2; i++) {
    if (part[i] < matrix[0+2*i] || part[i] > matrix[1+2*i]) return 0;
  }
  if (part[0] + part[1] > 2.0) return 0;
  return 1;
}

void stressTest()
{
  return;
}

//----------------------------------------------------------------------
//-- Functions related to summary statistics

void summary_multiscale(double_arr *peakList, hist_t *hist, double_mat *multiscale, double *summ)
{
  int i;
  for (i=0; i<multiscale->length; i++) summ[i] = multiscale->matrix[i];
  return;
}

//----------------------------------------------------------------------
//-- Functions related to distance

double dist_gauss(double *a, double *b)
{
  double invCovDiag[18] = {
    0.000695, 0.000879, 0.001538, 0.002896, 0.006451,
    0.018447, 0.047001, 0.119864, 0.097484, 0.001681,
    0.002385, 0.004400, 0.007965, 0.017598, 0.031376,
    0.077548, 0.160422, 0.099050
  };
  double sum = 0.0;
  int i;
  for (i=0; i<18; i++) sum += pow(a[i]-b[i], 2) * invCovDiag[i];
  return sqrt(sum);
}

double dist_star(double *a, double *b)
{
  double invCovDiag[27] = {
    0.000310, 0.000386, 0.000665, 0.001532, 0.004057,
    0.016789, 0.060073, 0.219620, 0.539488, 0.000918,
    0.001014, 0.001815, 0.003225, 0.007953, 0.019257,
    0.052212, 0.144693, 0.206096, 0.002380, 0.003144,
    0.005421, 0.009167, 0.019631, 0.037264, 0.098792,
    0.184801, 0.150118
  };
  double sum = 0.0;
  int i;
  for (i=0; i<27; i++) sum += pow(a[i]-b[i], 2) * invCovDiag[i];
  return sqrt(sum);
}

double dist_6D(double *a, double *b)
{
  return sqrt(pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2) + pow(a[2]-b[2], 2) + pow(a[3]-b[3], 2) + pow(a[4]-b[4], 2) + pow(a[5]-b[5], 2));
}

double dist_5D(double *a, double *b)
{
  return sqrt(pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2) + pow(a[2]-b[2], 2) + pow(a[3]-b[3], 2) + pow(a[4]-b[4], 2));
}

double dist_4D(double *a, double *b)
{
  return sqrt(pow(a[0]-b[0], 2) + pow(a[1]-b[1], 2) + pow(a[2]-b[2], 2) + pow(a[3]-b[3], 2));
}

double dist_1D(double *a, double *b)
{
  return fabs(a[0] - b[0]);
}

//----------------------------------------------------------------------
//-- Main functions

void doABC(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  peak->printMode = 3; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  //-- Initialization
  PMC_ABC_t *ABC = initialize_PMC_ABC_t(cmhm, peak, err); forwardError(*err, __LINE__,);
  fillObservation(peak->ABC_obsPath, peak, ABC, err); forwardError(*err, __LINE__,);
  if (peak->MPIInd == 0) print_PMC_ABC_t(peak, ABC);
  
  char name[STRING_LENGTH_MAX];
  
  //-- Iteration
  while (ABC->newIter->nbAttempts * ABC->r_stop <= ABC->Q) {
    runIteration(cmhm, peak, ABC, err); forwardError(*err, __LINE__,); //-- Switch old and new inside
    if (peak->MPIInd == 0) {
      sprintf(name, "iteration_%s_Q%d_t%d", STR_SUMM_T(ABC->summ), ABC->Q, ABC->t);
      outputIteration(name, ABC);
    }
#ifndef __releaseMenu__
    MPI_Barrier(MPI_COMM_WORLD); //-- Wait for output
#endif
  }
  
  if (peak->MPIInd == 0) printf("%d iterations done\n", ABC->t+1);
  free_PMC_ABC_t(ABC);
  return;
}

void doABC_subset(cosmo_hm *cmhm, peak_param *peak, int Q_real, int t, int subsetInd, error **err)
{
  testErrorRet(peak->MPISize>1, peak_badValue, "MPISize > 1", *err, __LINE__,);
  testErrorRet(peak->ABC_Q%Q_real!=0, peak_badValue, "Q % Q_real != 0", *err, __LINE__,);
  testErrorRet(peak->ABC_Q<=Q_real*subsetInd, peak_badValue, "subsetInd too large", *err, __LINE__,);
  peak->printMode = 3; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  char name[STRING_LENGTH_MAX];
  
  //-- Initialization
  PMC_ABC_t *ABC = initialize_PMC_ABC_t(cmhm, peak, err); forwardError(*err, __LINE__,);
  fillObservation(peak->ABC_obsPath, peak, ABC, err);     forwardError(*err, __LINE__,);
  sprintf(name, "iteration_%s_Q%d_t%d", STR_SUMM_T(ABC->summ), ABC->Q, t-1);
  fillIteration(name, peak, ABC, t-1, err);               forwardError(*err, __LINE__,);
  
  setRealQ(ABC, Q_real); //-- After fillIteration
  print_PMC_ABC_t_subset(peak, ABC);
  
  //-- Iteration
  peak->MPIInd = subsetInd; //-- To have different names for mrlens maps
  runIterationWithoutUpdate(cmhm, peak, ABC, err); forwardError(*err, __LINE__,);
  peak->MPIInd = 0;
  
  //-- Output
  sprintf(name, "iteration_%s_Q%d_t%d_subset%d", STR_SUMM_T(ABC->summ), Q_real, t, subsetInd);
  outputIteration(name, ABC);
  printf("------------------------------------------------------------------------\n");
  free_PMC_ABC_t(ABC);
  return;
}

void doABC_gather(cosmo_hm *cmhm, peak_param *peak, int Q_real, int t, error **err)
{
  testErrorRet(peak->MPISize>1, peak_badValue, "MPISize > 1", *err, __LINE__,);
  testErrorRet(peak->ABC_Q%Q_real!=0, peak_badValue, "Q % Q_real != 0", *err, __LINE__,);
  peak->printMode = 3; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  
  char name[STRING_LENGTH_MAX];
  
  //-- Initialization
  PMC_ABC_t *ABC = initialize_PMC_ABC_t(cmhm, peak, err); forwardError(*err, __LINE__,);
  sprintf(name, "iteration_%s_Q%d_t%d", STR_SUMM_T(ABC->summ), ABC->Q, t-1);
  fillIteration(name, peak, ABC, t-1, err); forwardError(*err, __LINE__,);
  print_PMC_ABC_t_subset(peak, ABC);
  
  swapOldAndNew(ABC);
  ABC->t++;
  setEpsilon(ABC);
  
  int nbSubsets = ABC->Q / Q_real;
  ABC->newIter->Q_MPI      = 0; //-- Considered as a counter
  ABC->newIter->nbAttempts = 0;
  int i;
  
  //-- Read
  for (i=0; i<nbSubsets; i++) {
    sprintf(name, "iteration_%s_Q%d_t%d_subset%d", STR_SUMM_T(ABC->summ), Q_real, t, i);
    readSubset(name, ABC->newIter, err); forwardError(*err, __LINE__,);
  }
  
  //-- Update
  setWeights(ABC);
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  
  //-- Output the gathered iteration
  sprintf(name, "iteration_%s_Q%d_t%d", STR_SUMM_T(ABC->summ), ABC->Q, t);
  outputIteration(name, ABC);
  printf("\n");
  printf("nbAttempts   = %d\n", ABC->newIter->nbAttempts);
  printf("success rate = %f\n", (double)ABC->Q / (double)ABC->newIter->nbAttempts);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

