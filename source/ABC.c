

  /*******************************************************
   **  ABC.c						**
   **  Version 2018.03.11				**
   **							**
   **  References:					**
   **  - Marin et al. (2011)				**
   **  - Weyant et al. (2013) - ApJ, 764, 116		**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


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

iteration_t *initialize_iteration_t(int f, int *doParam, int Q, error **err)
{
  iteration_t *iter = (iteration_t*)malloc_err(sizeof(iteration_t), err); forwardError(*err, __LINE__, NULL);
  iter->f           = f;
  iter->doParam     = (param_t*)malloc_err(f * sizeof(param_t), err);     forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<f; i++) iter->doParam[i] = (param_t)doParam[i];
  iter->nbParticles = 0;
  iter->nbAttempts  = 0;
  iter->epsilon     = ABC_TOLERANCE_MAX;
  iter->mean        = gsl_vector_alloc(f);
  iter->buffer      = gsl_vector_alloc(f);
  iter->cov         = gsl_matrix_alloc(f, f);
  iter->cov2        = gsl_matrix_alloc(f, f);
  iter->invCov      = gsl_matrix_alloc(f, f);
  iter->cholesky    = gsl_matrix_alloc(f, f);
  iter->perm        = gsl_permutation_alloc(f);
  iter->matrix      = initialize_double_mat(f+2, Q);
  return iter;
}

void free_iteration_t(iteration_t *iter)
{
  if (iter) {
    if (iter->doParam) {free(iter->doParam); iter->doParam = NULL;}
    gsl_vector_free(iter->mean);
    gsl_vector_free(iter->buffer);
    gsl_matrix_free(iter->cov);
    gsl_matrix_free(iter->cov2);
    gsl_matrix_free(iter->invCov);
    gsl_matrix_free(iter->cholesky);
    gsl_permutation_free(iter->perm);
    free_double_mat(iter->matrix);
    free(iter); iter = NULL;
  }
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
  printf("mean       = ");  printDoubleArray(iter->mean->data, iter->mean->size, 1.0, 3); printf("\n");
  printf("cov        =\n"); print_gsl_matrix(iter->cov);
  printf("invCov     =\n"); print_gsl_matrix(iter->invCov);
  
  int f2 = iter->f + 2;
  int L  = MIN(iter->nbParticles, 20);
  
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
  int Q = iter->nbParticles;
  fprintf(file, "# iteration_t\n");
  fprintf(file, "# f                   = %d\n", iter->f);
  fprintf(file, "# Number of particles = %d\n", Q);
  fprintf(file, "# Number of attempts  = %d\n", iter->nbAttempts);
  fprintf(file, "# Success rate        = %.5f\n", Q / (double)iter->nbAttempts);
  fprintf(file, "# epsilon             = %.5f\n", iter->epsilon);
  fprintf(file, "# cov                 =\n");
  output_gsl_matrix(file, iter->cov);
  fprintf(file, "# invCov              =\n");
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
  int Q  = iter->nbParticles;
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
  int Q  = iter->nbParticles;
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
  //-- TODO what is the unbiased estimate for weighted sum?
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

PMC_ABC_t *initialize_PMC_ABC_t(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  PMC_ABC_t *ABC      = (PMC_ABC_t*)malloc_err(sizeof(PMC_ABC_t), err);                                  forwardError(*err, __LINE__, NULL);
  ABC->f              = pkPar->ABC_f;
  ABC->Q              = pkPar->ABC_Q;
  ABC->r_stop         = pkPar->ABC_r_stop;
  
  ABC->nbAttempts_max = INT_MAX;
  ABC->Q_MPI          = (int)round((double)ABC->Q / (double)pkPar->MPISize);
  
  ABC->t              = -1; //-- To be increment later
  ABC->oldIter        = initialize_iteration_t(ABC->f, pkPar->ABC_doParam, ABC->Q_MPI*pkPar->MPISize, err); forwardError(*err, __LINE__, NULL);
  ABC->newIter        = initialize_iteration_t(ABC->f, pkPar->ABC_doParam, ABC->Q_MPI*pkPar->MPISize, err); forwardError(*err, __LINE__, NULL);
  ABC->deltaList      = initialize_double_arr(ABC->Q);
  
  ABC->priorFct       = prior_rectangle;
  ABC->limit          = initialize_double_mat(2, ABC->f);                                                forwardError(*err, __LINE__, NULL);
  ABC->flatness       = garanteeFlatness(ABC->f, ABC->newIter->doParam);
  fillLimitAndValue(chPar, ABC);
  
  int d_tot  = pkPar->N_nu * pkPar->nbFilters;
  ABC->x_obs          = initialize_double_arr(d_tot);
  ABC->x_mod          = initialize_double_arr(d_tot);
  ABC->summFct        = summary_double;
  
  ABC->distFct        = (pkPar->ABC_doCorr == 1) ? dist_chi : dist_uncorrChi;
  ABC->invCov         = gsl_matrix_alloc(d_tot, d_tot);
  ABC->Delta_x        = gsl_vector_alloc(d_tot);
  ABC->intermediate   = gsl_vector_alloc(d_tot);
  
  if (pkPar->MPIInd == 0) printf("ABC initialization done\n");
  return ABC;
}

void free_PMC_ABC_t(PMC_ABC_t *ABC)
{
  if (ABC) {
    if (ABC->oldIter)   {free_iteration_t(ABC->oldIter);  ABC->oldIter   = NULL;}
    if (ABC->newIter)   {free_iteration_t(ABC->newIter);  ABC->newIter   = NULL;}
    if (ABC->deltaList) {free_double_arr(ABC->deltaList); ABC->deltaList = NULL;}
    
    if (ABC->limit)     {free_double_mat(ABC->limit);     ABC->limit     = NULL;}
    
    if (ABC->x_obs)     {free_double_arr(ABC->x_obs);     ABC->x_obs     = NULL;}
    if (ABC->x_mod)     {free_double_arr(ABC->x_mod);     ABC->x_mod     = NULL;}
    gsl_matrix_free(ABC->invCov);
    gsl_vector_free(ABC->Delta_x);
    gsl_vector_free(ABC->intermediate);
    
    free(ABC); ABC = NULL;
  }
  return;
}

void print_PMC_ABC_t(peak_param *pkPar, PMC_ABC_t *ABC)
{
  int i;
  printf("\nABC parameters\n");
  printf("Dimension of parameter   = %d\n", ABC->f);
  printf("Parameter space          = %s", STR_PARAM_T(ABC->newIter->doParam[0])); for (i=1; i<ABC->f; i++) printf(", %s", STR_PARAM_T(ABC->newIter->doParam[i])); printf("\n");
  printf("Number of particles      = %d = %d procs * %d - %d\n", ABC->Q, pkPar->MPISize, ABC->Q_MPI, pkPar->MPISize * ABC->Q_MPI - ABC->Q);
  printf("Shutoff success rate     = %.3f\n", ABC->r_stop);
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

void fillLimitAndValue(cosmo_hm *chPar, PMC_ABC_t *ABC)
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
	value[i]      = chPar->cosmo->Omega_m;
	break;
      case param_Omega_de:
	buffer[0+2*i] = OMEGA_DE_MIN;
	buffer[1+2*i] = OMEGA_DE_MAX;
	value[i]      = chPar->cosmo->Omega_de;
	break;
      case param_Omega_b:
	buffer[0+2*i] = OMEGA_B_MIN;
	buffer[1+2*i] = OMEGA_B_MAX;
	value[i]      = chPar->cosmo->Omega_b;
	break;
      case param_n_s:
	buffer[0+2*i] = N_S_MIN;
	buffer[1+2*i] = N_S_MAX;
	value[i]      = chPar->cosmo->n_spec;
	break;
      case param_h_100:
	buffer[0+2*i] = H_100_MIN;
	buffer[1+2*i] = H_100_MAX;
	value[i]      = chPar->cosmo->h_100;
	break;
      case param_sigma_8:
	buffer[0+2*i] = SIGMA_8_MIN;
	buffer[1+2*i] = SIGMA_8_MAX;
	value[i]      = chPar->cosmo->normalization;
	break;
      case param_w0_de:
	buffer[0+2*i] = W0_DE_MIN;
	buffer[1+2*i] = W0_DE_MAX;
	value[i]      = chPar->cosmo->w0_de;
	break;
      case param_w1_de:
	buffer[0+2*i] = W1_DE_MIN;
	buffer[1+2*i] = W1_DE_MAX;
	value[i]      = chPar->cosmo->w1_de;
	break;
      case param_c_0:
	buffer[0+2*i] = C_0_MIN;
	buffer[1+2*i] = C_0_MAX;
	value[i]      = chPar->c0;
	break;
      case param_beta_NFW:
	buffer[0+2*i] = BETA_NFW_MIN;
	buffer[1+2*i] = BETA_NFW_MAX;
	value[i]      = chPar->beta_NFW;
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

void readObservation(char name[], peak_param *pkPar, PMC_ABC_t *ABC, error **err)
{
  //-- Open
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  int d_tot  = ABC->x_obs->length;
  int count = 0;
  
  char buffer[STRING_LENGTH_MAX];
  double value;
  int buffer2;
  
  //-- Read
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf", &value);
      if (buffer2 < 1) break;
      testErrorRetVA(count>=d_tot, peak_overflow, "Too many values in \"%s\"", *err, __LINE__, , name);
      
      ABC->x_obs->array[count] = value;
      count++;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  //-- Check
  testErrorRetVA(count<d_tot, peak_match, "Too few values in \"%s\"", *err, __LINE__, , name);
  if (pkPar->MPIInd == 0) printf("\"%s\" read\n", name);
  return;
}

void readIteration(char name[], peak_param *pkPar, PMC_ABC_t *ABC, int t, error **err)
{
  ABC->t = t;
  if (t < 0) return;
  
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX];
  int f2 = ABC->f + 2;
  double *part = ABC->newIter->matrix->matrix;
  int i;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
    else {
      testErrorRet(ABC->newIter->nbParticles>=ABC->Q, peak_overflow, "Too many particles", *err, __LINE__,);
      ungetc(c, file);
      for (i=0; i<f2; i++) fscanf(file, "%lf", &part[i]);
      fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
      ABC->newIter->nbParticles++;
      part += f2;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  testErrorRet(ABC->newIter->nbParticles<ABC->Q, peak_match, "Too few particles", *err, __LINE__,);
  
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
  
  if (pkPar->MPIInd == 0) {
    printf("\"%s\" read\n", name);
    printf("%d particles generated\n", ABC->Q);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to accept-reject algorithm

void generateParam(PMC_ABC_t *ABC, gsl_rng *generator, double *newPart, prior_fct *prior)
{
  int f  = ABC->f;
  int f2 = f + 2;
  double *oldBegin     = ABC->oldIter->matrix->matrix;
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

cosmo_hm *initialize_cosmo_hm_ABC(cosmo_hm *oldChPar, PMC_ABC_t *ABC, double *newPart, error **err)
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
  
  cosmo *oldCm       = oldChPar->cosmo;
  redshift_t *oldRs  = oldChPar->redshift;
  cosmo_hm *newChPar = init_parameters_hm(
      value[param_Omega_m],  value[param_Omega_de], value[param_w0_de],       value[param_w1_de],       oldCm->w_poly_de,          oldCm->N_poly_de,
      value[param_h_100],    value[param_Omega_b],  oldCm->Omega_nu_mass,     oldCm->Neff_nu_mass,      value[param_sigma_8],      value[param_n_s],
      oldRs->Nzbin,          oldRs->Nnz,            oldRs->nofz,              oldRs->photz,             oldRs->par_nz,             oldChPar->zmin,     oldChPar->zmax,
      oldCm->nonlinear,      oldCm->transfer,       oldCm->growth,            oldCm->de_param,          oldChPar->cosmo->normmode,  
      value[param_c_0],      oldChPar->alpha_NFW,   value[param_beta_NFW],    oldChPar->massfct,        oldChPar->halo_bias,   
      oldChPar->log10M_min,  oldChPar->log10M1,     oldChPar->log10M0,        oldChPar->sigma_log_M,    oldChPar->alpha, 
      oldChPar->log10Mstar0, oldChPar->beta,        oldChPar->delta,          oldChPar->gamma,          oldChPar->B_cut,           oldChPar->B_sat, 
      oldChPar->beta_cut,    oldChPar->beta_sat,    oldChPar->log10Mstar_min, oldChPar->log10Mstar_max, oldChPar->eta,
      oldChPar->fcen1,       oldChPar->fcen2,
      oldChPar->hod,         oldChPar->pi_max,      err);
  forwardError(*err, __LINE__, NULL);
  return newChPar;
}

void generateModel(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, double *newPart, error **err)
{
  cosmo_hm *chParABC = initialize_cosmo_hm_ABC(chPar, ABC, newPart, err);                       forwardError(*err, __LINE__,);
  setMassSamplers(chParABC, pkPar, pipe->hSampArr, 1, err);                                     forwardError(*err, __LINE__,); //-- doVolume = 1
  if (pkPar->doRandGalPos == 0) updateCosmo_gal_map(chParABC, pkPar, pipe->gMap, err);          forwardError(*err, __LINE__,);
  
  massFctToMultiscale(chParABC, pkPar, pipe->hSampArr, pipe->k1Inter, pipe->hMap, pipe->gSamp, pipe->gMap, pipe->mask, pipe->FFTSmoother, pipe->DCSmoother, pipe->kMap, 
		      pipe->variance, pipe->peakList, pipe->nuHist, pipe->multiscale, 1, err); forwardError(*err, __LINE__,); //-- nbPatches = 1
  
  ABC->summFct(NULL, (void*)pipe->multiscale->matrix, pipe->multiscale->length, ABC->x_mod->array);
  newPart[ABC->f] = ABC->distFct(ABC->x_obs->array, ABC->x_mod->array, ABC->invCov, ABC->Delta_x, ABC->intermediate); 
  free_parameters_hm(&chParABC);
  return;
}

int oneSampleTest(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, double *newPart, error **err)
{
  int f = ABC->f;
  param_t *doParam = ABC->newIter->doParam;
  char buffer1[STRING_LENGTH_MAX], buffer2[16];
  
  //-- Generate particle
  if (ABC->t == 0) {
    do priorGenerator(pkPar->generator, ABC->limit, newPart);
    while (!ABC->priorFct(ABC->limit, newPart));
  }
  else generateParam(ABC, pkPar->generator, newPart, ABC->priorFct);
  
  //-- Make string
  particleString(buffer1, buffer2, f, doParam, newPart);
  if (pkPar->verbose < 5) {printf("%s, ", buffer1); fflush(stdout);}
  
  //-- Generate model
  generateModel(chPar, pkPar, pipe, ABC, newPart, err);
  
  //-- Throw error
  if (isError(*err)) {
    if      (pkPar->verbose < 5)  printf("get error and resample\n");
    else if (pkPar->verbose == 5) printf("Proc %2d: %s, get error and resample\n", pkPar->MPIInd, buffer1);
    purgeError(err);
    return -1; //-- To resample
  }
  
  //-- Accept-reject
  ABC->newIter->nbAttempts += 1;
  
  if (newPart[f] <= ABC->newIter->epsilon) {
    if      (pkPar->verbose < 5)  printf("delta = %7.3f, accepted\n", newPart[f]);
    else if (pkPar->verbose == 5) printf("Proc %2d: %s, %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, accepted\n", 
      pkPar->MPIInd, buffer1, pipe->hMap->total, pipe->gMap->total, pipe->gMap->kappa_mean, pipe->gMap->fillingThreshold, newPart[f]);
    return 1; //-- Accepted
  }
  
  if      (pkPar->verbose < 5)  printf("delta = %7.3f, rejected\n", newPart[f]);
  else if (pkPar->verbose == 5) printf("Proc %2d: %s, %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, rejected\n", 
    pkPar->MPIInd, buffer1, pipe->hMap->total, pipe->gMap->total, pipe->gMap->kappa_mean, pipe->gMap->fillingThreshold, newPart[f]);
  return 0; //-- Rejected
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

double chiSquared(double *x1, double *x2, gsl_matrix *invCov, gsl_vector *Delta_x, gsl_vector *intermediate)
{
  double value;
  int i;
  for (i=0; i<Delta_x->size; i++) Delta_x->data[i] = x2[i] - x1[i];
  gsl_blas_dsymv(CblasUpper, 1.0, invCov, Delta_x, 0.0, intermediate); //-- intermediate = invCov * Delta_x
  gsl_blas_ddot(Delta_x, intermediate, &value);
  return value;
}

void setWeights(PMC_ABC_t *ABC)
{
  int f2 = ABC->f + 2;
  int Q  = ABC->Q;
  double *oldBegin = ABC->oldIter->matrix->matrix;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *oldEnd   = oldBegin + f2 * Q;
  double *newEnd   = newBegin + f2 * Q;
  
  double *newPart;
  
  if (ABC->t == 0) {
    for (newPart=newBegin; newPart<newEnd; newPart+=f2) newPart[f2-1] = 1.0 / (double)Q;
    return;
  }
  
  gsl_matrix *invCov       = ABC->oldIter->invCov;
  gsl_vector *Delta_param  = ABC->oldIter->buffer;
  gsl_vector *intermediate = ABC->newIter->buffer;
  double sum = 0.0;
  
  double *oldPart;
  double arg, weight;
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    weight = 0.0;
    for (oldPart=oldBegin; oldPart<oldEnd; oldPart+=f2) {
      arg = -0.5 * chiSquared(oldPart, newPart, invCov, Delta_param, intermediate);
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

void runIteration(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  swapOldAndNew(ABC);
  ABC->t++;
  ABC->newIter->nbParticles = 0;
  ABC->newIter->nbAttempts  = 0;
  setEpsilon(ABC);
  
  if (pkPar->MPIInd == 0) {
    printf("-- Begin iteration %d\n\n", ABC->t);
    if (ABC->t == 0) printf("Tolerance = DBL_MAX\n\n");
    else             printf("Tolerance = %.5f\n\n", ABC->newIter->epsilon);
  }
#ifdef __CAMELUS_USE_MPI__
  MPI_Barrier(MPI_COMM_WORLD); //-- Wait for print
#endif
  
  int f2    = ABC->f + 2;
  int Q     = ABC->Q;
  int Q_MPI = ABC->Q_MPI;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *newPart  = newBegin;
  int accept;
  
  //-- Size ajustment for MPI
  if (pkPar->MPIInd == pkPar->MPISize - 1) Q_MPI -= pkPar->MPISize * Q_MPI - Q;
  
  do {
    accept = oneSampleTest(chPar, pkPar, pipe, ABC, newPart, err); forwardError(*err, __LINE__,);
    if (accept == 1) {
      ABC->newIter->nbParticles += 1;
      newPart += f2;
    }
  }
  while (ABC->newIter->nbParticles < Q_MPI);
  printf("Proc %2d finished\n", pkPar->MPIInd);
  
#ifdef __CAMELUS_USE_MPI__
  MPI_Allgather(newBegin, f2 * Q_MPI, MPI_DOUBLE, newBegin, f2 * Q_MPI, MPI_DOUBLE, MPI_COMM_WORLD);        //-- Communicate new particles
  MPI_Allreduce(&ABC->newIter->nbAttempts, &ABC->newIter->nbAttempts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //-- Communicate the number of attempts
#endif
  
  setWeights(ABC);
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
  
  if (pkPar->MPIInd == 0) {
    printf("\nSuccess rate = %.5f\n", (double)ABC->newIter->nbParticles / (double)ABC->newIter->nbAttempts);
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

void runIterationWithoutUpdate(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, error **err)
{
  //-- Assume that MPISize = 1
  
  clock_t start = clock();
  swapOldAndNew(ABC);
  ABC->t++;
  ABC->newIter->nbParticles = 0;
  ABC->newIter->nbAttempts  = 0;
  setEpsilon(ABC);
  
  printf("-- Begin iteration %d\n\n", ABC->t);
  if (ABC->t == 0) printf("Tolerance = DBL_MAX\n\n");
  else             printf("Tolerance = %.5f\n\n", ABC->newIter->epsilon);
  
  int f2 = ABC->f + 2;
  double *newPart = ABC->newIter->matrix->matrix;
  int accept;
  
  do {
    accept = oneSampleTest(chPar, pkPar, pipe, ABC, newPart, err); forwardError(*err, __LINE__,);
    if (accept == 1) {
      ABC->newIter->nbParticles += 1;
      newPart += f2;
    }
  }
  while (ABC->newIter->nbAttempts < ABC->nbAttempts_max);
  
  printf("\nSuccess rate = %.5f\n", (double)ABC->newIter->nbParticles / (double)ABC->newIter->nbAttempts);
  printf("\n-- End iteration %d, ", ABC->t); routineTime(start, clock());
  printf("------------------------------------------------------------------------\n");
  return;
}

void readSubset(char name[], iteration_t *iter, error **err)
{
  int Q = iter->matrix->N2;
  if (iter->nbParticles == Q) {
    printf("\"%s\" skipped (already got %d particles)\n", name, Q);
    return;
  }
  
  FILE *file = fopen_err(name, "r", err); 
  if (isError(*err)) {
    printf("\"%s\" skipped (file not existed)\n", name);
    purgeError(err);
    return;
  }
  
  char buffer[STRING_LENGTH_MAX];
  int count;
  
  fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
  fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
  fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
  fscanf(file, "# Number of attempts  = %d\n", &count); //-- Read the 4th line which contains nbAttempts
  iter->nbAttempts += count;
  count = 0;
  
  int f2 = iter->f + 2;
  double *part = iter->matrix->matrix + f2 * iter->nbParticles;
  int i;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#')               fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
    else if (iter->nbParticles >= Q) fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
    else {
      ungetc(c, file);
      for (i=0; i<f2; i++) fscanf(file, "%lf", &part[i]);
      fgets(buffer, STRING_LENGTH_MAX, file); //-- Clean the line
      iter->nbParticles++;
      count++;
      part += f2;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  printf("\"%s\" read (%d particles)\n", name, count);
  return;
}

void print_PMC_ABC_t_subset(peak_param *pkPar, PMC_ABC_t *ABC)
{
  int i;
  printf("\nABC parameters\n");
  printf("Dimension of parameter        = %d\n", ABC->f);
  printf("Parameter space               = %s", STR_PARAM_T(ABC->newIter->doParam[0])); for (i=1; i<ABC->f; i++) printf(", %s", STR_PARAM_T(ABC->newIter->doParam[i])); printf("\n");
  printf("Number of particles           = %d\n", ABC->Q);
  printf("Number of attempts per subset = ");
  if (ABC->nbAttempts_max == INT_MAX) printf("+infty\n");
  else printf("%d\n", ABC->nbAttempts_max);
  printf("Dimension of data vector      = %d\n", ABC->x_obs->length);
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

//----------------------------------------------------------------------
//-- Functions related to summary statistics

void summary_int(void *arg0, void *arg1, int length, double *summ)
{
  int *array = (int*)arg1;
  int i;
  for (i=0; i<length; i++) summ[i] = (double)array[i];
  return;
}

void summary_double(void *arg0, void *arg1, int length, double *summ)
{
  double *array = (double*)arg1;
  int i;
  for (i=0; i<length; i++) summ[i] = array[i];
  return;
}

//----------------------------------------------------------------------
//-- Functions related to distance

double dist_chi(double *x1, double *x2, gsl_matrix *invCov, gsl_vector *Delta_x, gsl_vector *intermediate)
{
  return sqrt(chiSquared(x1, x2, invCov, Delta_x, intermediate));
}

double dist_uncorrChi(double *x1, double *x2, gsl_matrix *invCov, gsl_vector *Delta_x, gsl_vector *intermediate)
{
  int d = Delta_x->size;
  double sum = 0.0;
  int i;
  for (i=0; i<d; i++) sum += SQ(x2[i] - x1[i]) * invCov->data[i+i*d];
  return sqrt(sum);
}

void readInvCov(char name[], peak_param *pkPar, gsl_matrix *invCov, error **err)
{
  //-- Open
  FILE *file   = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  int d_tot    = invCov->size1;
  int nbValues = (pkPar->ABC_doCorr == 1) ? SQ(d_tot) : d_tot;
  int count    = 0;
  
  char buffer[STRING_LENGTH_MAX];
  double value;
  int ind, buffer2;
  
  //-- Read
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf", &value);
      if (buffer2 < 1) break;
      testErrorRetVA(count>=nbValues, peak_overflow, "Too many values in \"%s\"", *err, __LINE__, , name);
      
      ind = (pkPar->ABC_doCorr == 1) ? count : (count + count * d_tot);
      invCov->data[ind] = value;
      count++;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  //-- Check
  testErrorRetVA(count<nbValues, peak_match, "Too few values in \"%s\"", *err, __LINE__, , name);
  if (pkPar->MPIInd == 0) printf("\"%s\" read\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doABC(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  pkPar->verbose = 5;
  
  //-- Initialization
  pipeline_t *pipe = initialize_pipeline_t(chPar, pkPar, err);      forwardError(*err, __LINE__,);
  PMC_ABC_t *ABC   = initialize_PMC_ABC_t(chPar, pkPar, err);       forwardError(*err, __LINE__,);
  readObservation(pkPar->ABC_obsPath, pkPar, ABC, err);             forwardError(*err, __LINE__,);
  readInvCov(pkPar->ABC_invCovPath, pkPar, ABC->invCov, err);       forwardError(*err, __LINE__,);
  if (pkPar->MPIInd == 0) print_PMC_ABC_t(pkPar, ABC);
  
  char name[STRING_LENGTH_MAX];
  
  //-- Iteration
  while (ABC->newIter->nbAttempts * ABC->r_stop <= ABC->Q) {
    runIteration(chPar, pkPar, pipe, ABC, err);                     forwardError(*err, __LINE__,); //-- Switch old and new inside
    if (pkPar->MPIInd == 0) {
      sprintf(name, "iteration_Q%d_t%d", ABC->Q, ABC->t);
      outputIteration(name, ABC);
    }
#ifdef __CAMELUS_USE_MPI__
    MPI_Barrier(MPI_COMM_WORLD); //-- Wait for output
#endif
  }
  
  if (pkPar->MPIInd == 0) printf("%d iterations done\n", ABC->t+1);
  free_pipeline_t(pipe);
  free_PMC_ABC_t(ABC);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doABC_subset(cosmo_hm *chPar, peak_param *pkPar, int t, int nbAttempts, int subsetInd, error **err)
{
  testErrorRet(pkPar->MPISize>1, peak_badValue, "MPISize > 1", *err, __LINE__,);
  pkPar->verbose = 3;
  
  char name[STRING_LENGTH_MAX];
  
  //-- Initialization
  pipeline_t *pipe = initialize_pipeline_t(chPar, pkPar, err);      forwardError(*err, __LINE__,);
  PMC_ABC_t *ABC   = initialize_PMC_ABC_t(chPar, pkPar, err);       forwardError(*err, __LINE__,);
  ABC->nbAttempts_max = nbAttempts;
  
  readObservation(pkPar->ABC_obsPath, pkPar, ABC, err);             forwardError(*err, __LINE__,);
  readInvCov(pkPar->ABC_invCovPath, pkPar, ABC->invCov, err);       forwardError(*err, __LINE__,);
  sprintf(name, "iteration_Q%d_t%d", ABC->Q, t-1);
  readIteration(name, pkPar, ABC, t-1, err);                        forwardError(*err, __LINE__,);
  print_PMC_ABC_t_subset(pkPar, ABC);
  
  //-- Iteration
  pkPar->MPIInd = subsetInd; //-- To have different names for mrlens maps
  runIterationWithoutUpdate(chPar, pkPar, pipe, ABC, err);          forwardError(*err, __LINE__,);
  pkPar->MPIInd = 0;
  
  //-- Output
  sprintf(name, "iteration_t%d_subset%d", t, subsetInd);
  outputIteration(name, ABC);
  
  free_pipeline_t(pipe);
  free_PMC_ABC_t(ABC);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doABC_gather(cosmo_hm *chPar, peak_param *pkPar, int t, int nbSubsets, error **err)
{
  testErrorRet(pkPar->MPISize>1, peak_badValue, "MPISize > 1", *err, __LINE__,);
  pkPar->verbose = 3;
  
  char name[STRING_LENGTH_MAX];
  
  //-- Initialization
  PMC_ABC_t *ABC = initialize_PMC_ABC_t(chPar, pkPar, err); forwardError(*err, __LINE__,);
  sprintf(name, "iteration_Q%d_t%d", ABC->Q, t-1);
  readIteration(name, pkPar, ABC, t-1, err);                forwardError(*err, __LINE__,);
  print_PMC_ABC_t_subset(pkPar, ABC);
  
  swapOldAndNew(ABC);
  ABC->t++;
  setEpsilon(ABC);
  
  ABC->newIter->nbParticles = 0;
  ABC->newIter->nbAttempts  = 0;
  int i;
  
  //-- Read
  for (i=0; i<nbSubsets; i++) {
    sprintf(name, "iteration_t%d_subset%d", t, i);
    readSubset(name, ABC->newIter, err);                    forwardError(*err, __LINE__,);
  }
  if (ABC->newIter->nbParticles < ABC->Q) {
    printf("Still miss %d particles, nothing is done!\n", ABC->Q-ABC->newIter->nbParticles);
    printf("\n");
    printf("Number of attempts = %d\n", ABC->newIter->nbAttempts);
    printf("Success rate = %.5f\n", (double)ABC->newIter->nbParticles / (double)ABC->newIter->nbAttempts);
    free_PMC_ABC_t(ABC);
    printf("------------------------------------------------------------------------\n");
    return;
  }
  
  //-- Update
  setWeights(ABC);
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  
  //-- Output the gathered iteration
  sprintf(name, "iteration_Q%d_t%d", ABC->Q, t);
  outputIteration(name, ABC);
  
  printf("\n");
  printf("Number of attempts = %d\n", ABC->newIter->nbAttempts);
  printf("Success rate = %.5f\n", (double)ABC->newIter->nbParticles / (double)ABC->newIter->nbAttempts);
  free_PMC_ABC_t(ABC);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

