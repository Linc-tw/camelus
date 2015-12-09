

  /***********************************************
   **  ABC.c					**
   **  Chieh-An Lin				**
   **  Version 2015.12.09			**
   **						**
   **  References:				**
   **  - Marin et al. (2011)			**
   **  - Weyant et al. (2013) - ApJ, 764, 116	**
   ***********************************************/


#include "ABC.h"


//----------------------------------------------------------------------
//-- Functions related to iteration_t

iteration_t *initialize_iteration_t(int f, int Q, int MPISize, error **err)
{
  iteration_t *iter = (iteration_t*)malloc_err(sizeof(iteration_t), err); forwardError(*err, __LINE__,);
  iter->f           = f;
  iter->Q           = Q;
  iter->Q_MPI       = (int)ceil((double)Q / (double)MPISize);
  iter->nbAttempts  = 0;
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

void print_iteration_t(iteration_t *iter)
{
  //-- WARNING: f-dependent
  
  printf("# iteration_t (first 20 elements)\n");
  int f2 = iter->f + 2;
  int L  = MIN(iter->Q, 20);
  
  double *begin = iter->matrix->matrix;
  double *end   = begin + f2 * L;
  double *part;
  for (part=begin; part<end; part+=f2) printf("(Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), weight = %.3f, delta = %9.5f\n", part[0], part[1], part[2], part[f2-2], part[f2-1]);
  return;
}

void output_gsl_matrix(FILE *file, gsl_matrix *mat)
{
  int f = mat->size1;
  int i, j;
  
  fprintf(file, "# [[% 9.5f", mat->data[0+0*f]);
  for (i=1; i<f; i++) fprintf(file, ", % 9.5f", mat->data[i+0*f]);
  fprintf(file, "]");
  
  for (j=1; j<f; j++) {
    fprintf(file, ",\n#  [% 9.5f", mat->data[0+j*f]);
    for (i=1; i<f; i++) fprintf(file, ", % 9.5f", mat->data[i+j*f]);
    fprintf(file, "]");
  }
  
  fprintf(file, "]\n");
  return;
}

void output1_iteration_t(FILE *file, iteration_t *iter)
{
  //-- WARNING: f-dependent
  
  int Q = iter->Q;
  fprintf(file, "# iteration_t\n");
  fprintf(file, "# f                  = %d\n", iter->f);
  fprintf(file, "# Q                  = %d\n", Q);
  fprintf(file, "# Number of attempts = %d\n", iter->nbAttempts);
  fprintf(file, "# Success rate       = %.5f\n", Q / (double)iter->nbAttempts);
  fprintf(file, "#\n");
  
  int f2        = iter->f + 2;
  double *begin = iter->matrix->matrix;
  double *end   = begin + f2 * Q;
  double *part;
  
  fprintf(file, "Omega_M = [");
  for (part=begin;   part<end; part+=f2) fprintf(file, ", %.5f", part[0]);
  fprintf(file, "]\n");
  
  fprintf(file, "sigma_8 = [");
  for (part=begin+1; part<end; part+=f2) fprintf(file, ", %.5f", part[0]);
  fprintf(file, "]\n");
  
  fprintf(file, "w0_de   = [");
  for (part=begin+2; part<end; part+=f2) fprintf(file, ", %.5f", part[0]);
  fprintf(file, "]\n");
  
  fprintf(file, "delta   = [");
  for (part=begin+3; part<end; part+=f2) fprintf(file, ", %.5f", part[f2-2]);
  fprintf(file, "]\n");
  
  fprintf(file, "weight  = [");
  for (part=begin+4; part<end; part+=f2) fprintf(file, ", %.5f", part[f2-1]);
  fprintf(file, "]\n");
  fprintf(file, "\n");
  
  fprintf(file, "# cov    =\n");
  output_gsl_matrix(file, iter->cov);
  fprintf(file, "# invCov =\n");
  output_gsl_matrix(file, iter->invCov);
  return;
}

void output2_iteration_t(FILE *file, iteration_t *iter)
{
  //-- WARNING: f-dependent
  
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
  fprintf(file, "# Omega_M  sigma_8     w0_de     delta    weight\n");
  
  int f2        = iter->f + 2;
  double *begin = iter->matrix->matrix;
  double *end   = begin + f2 * Q;
  double *part;
  
  for (part=begin; part<end; part+=f2) fprintf(file, "  %.5f  %.5f  % .5f  %9.5f  %.5f\n", part[0], part[1], part[2], part[f2-2], part[f2-1]);
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

void read_iteration_t(char name[], iteration_t *iter, double *deltaArr, error **err)
{
  //-- WARNING: f-dependent
  
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  int f2 = iter->f + 2;
  double *part = iter->matrix->matrix;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      testErrorRet(count==iter->Q, peak_overflow, "Too many particles", *err, __LINE__,);
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf %lf %lf %lf %lf\n", &part[0], &part[1], &part[2], &part[f2-2], &part[f2-1]);
      count++;
      part += f2;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  updateMean_iteration_t(iter);
  updateCovariance_iteration_t(iter);
  updateCholesky_iteration_t(iter);
  
  printf("\"%s\" read\n", name);
  printf("%d particles generated\n", count);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to PMC_ABC_t

PMC_ABC_t *initialize_PMC_ABC_t(peak_param *peak, error **err)
{
  PMC_ABC_t *ABC = (PMC_ABC_t*)malloc_err(sizeof(PMC_ABC_t), err); forwardError(*err, __LINE__,);
  ABC->Q         = peak->ABC_Q;
  ABC->r_stop    = peak->ABC_r_stop;
  int i;
  STRING2ENUM(ABC->summ, peak->ABC_summ, summary_t, STR_SUMMARY_T, i, NB_SUMMARY_T, err); forwardError(*err, __LINE__,);
  
  ABC->f         = 3;
  ABC->priorFct  = prior_rectangle;
  
  ABC->Q_MPI     = (int)ceil((double)ABC->Q / (double)peak->MPISize);
  
  ABC->t         = 0;
  ABC->oldIter   = initialize_iteration_t(ABC->f, ABC->Q, peak->MPISize, err); forwardError(*err, __LINE__,);
  ABC->newIter   = initialize_iteration_t(ABC->f, ABC->Q, peak->MPISize, err); forwardError(*err, __LINE__,);
  ABC->deltaList = initialize_double_arr(ABC->Q);
  ABC->newIter->epsilon = 1e+15;
  
  int N_bin      = peak->doSmoothing < 3 ? peak->N_nu : (peak->doSmoothing == 3 ? peak->N_kappa : MAX(peak->N_nu, peak->N_kappa));
  ABC->x_obs     = initialize_double_arr(N_bin * peak->nbFilters);
  ABC->x_mod     = initialize_double_arr(N_bin * peak->nbFilters);
  ABC->summFct   = summary_multiscale;
  
  //-- WARNING: For Paper III, the distance is defined somewhere else
  if (ABC->summ == summ_gauss)     ABC->distFct = dist_gauss;
  else if (ABC->summ == summ_star) ABC->distFct = dist_star;
  else if (ABC->summ == summ_mrlens);
  else {*err = addError(peak_unknown, "Unknown summary type", *err, __LINE__); forwardError(*err, __LINE__,);}
  
  int N1_mask    = (int)(peak->Omega[0] * peak->theta_CCD_inv);
  int N2_mask    = (int)(peak->Omega[1] * peak->theta_CCD_inv);
  int length     = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  cosmo_hm *cmhm = initialize_cosmo_hm_default(err);                                           forwardError(*err, __LINE__,);
  
  //-- Camelus pipeline
  ABC->sampArr    = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
  ABC->hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  ABC->galSamp    = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, ABC->galSamp, err);                                             forwardError(*err, __LINE__,);
  ABC->gMap       = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  ABC->CCDMask    = initialize_short_mat(N1_mask, N2_mask);
  if (peak->doMask == 1) fillMask_CFHTLenS_W1(peak, ABC->CCDMask);
  ABC->smoother   = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernel(peak, ABC->smoother);
  ABC->kMap       = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  ABC->variance   = initialize_FFT_arr(peak->smootherSize, peak->FFTSize);
  if (peak->doSmoothing == 1 || peak->doSmoothing == 4) makeKernelForVariance(peak, ABC->variance);
  ABC->peakList   = initialize_double_arr(length);
  ABC->hist       = initialize_hist_t(peak->N_nu);
  setHist(peak, ABC->hist);
  ABC->hist2      = initialize_hist_t(peak->N_kappa);
  setHist2(peak, ABC->hist2);
  ABC->multiscale = initialize_double_mat(N_bin, peak->nbFilters);
  
  free_parameters_hm(&cmhm);
  
  if (peak->MPIInd == 0) printf("ABC initialization done\n");
  return ABC;
}

void free_PMC_ABC_t(PMC_ABC_t *ABC)
{
  if (ABC->oldIter)    {free_iteration_t(ABC->oldIter);   ABC->oldIter = NULL;}
  if (ABC->newIter)    {free_iteration_t(ABC->newIter);   ABC->newIter = NULL;}
  if (ABC->deltaList)  {free_double_arr(ABC->deltaList);  ABC->deltaList = NULL;}
  
  if (ABC->x_obs)      {free_double_arr(ABC->x_obs);      ABC->x_obs = NULL;}
  if (ABC->x_mod)      {free_double_arr(ABC->x_mod);      ABC->x_mod = NULL;}
  
  if (ABC->sampArr)    {free_sampler_arr(ABC->sampArr);   ABC->sampArr = NULL;}
  if (ABC->hMap)       {free_halo_map(ABC->hMap);         ABC->hMap = NULL;}
  if (ABC->galSamp)    {free_sampler_t(ABC->galSamp);     ABC->galSamp = NULL;}
  if (ABC->gMap)       {free_gal_map(ABC->gMap);          ABC->gMap = NULL;}
  if (ABC->CCDMask)    {free_short_mat(ABC->CCDMask);     ABC->CCDMask = NULL;}
  if (ABC->smoother)   {free_FFT_arr(ABC->smoother);      ABC->smoother = NULL;}
  if (ABC->kMap)       {free_map_t(ABC->kMap);            ABC->kMap = NULL;}
  if (ABC->variance)   {free_FFT_arr(ABC->variance);      ABC->variance = NULL;}
  if (ABC->peakList)   {free_double_arr(ABC->peakList);   ABC->peakList = NULL;}
  if (ABC->hist)       {free_hist_t(ABC->hist);           ABC->hist = NULL;}
  if (ABC->hist2)      {free_hist_t(ABC->hist2);          ABC->hist2 = NULL;}
  if (ABC->multiscale) {free_double_mat(ABC->multiscale); ABC->multiscale = NULL;}
  
  free(ABC); ABC = NULL;
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
  else if (ABC->summ == summ_mrlens) {
    for (i=0; i<27; i++) ABC->x_obs->array[i] = buffer3[i+18];
  }
  
  if (peak->MPIInd == 0) printf("Observation data read\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to accept-reject algorithm

void generateParam(peak_param *peak, iteration_t *oldIter, double *newPart, prior_fct *prior)
{
  int f  = oldIter->f;
  int f2 = oldIter->f + 2;
  double *oldBegin     = oldIter->matrix->matrix;
  gsl_rng *generator   = peak->generator;
  gsl_vector *buffer   = oldIter->buffer;
  gsl_matrix *cholesky = oldIter->cholesky;
  
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
    for(i=0; i<f; i++) buffer->data[i] = gsl_ran_ugaussian(generator);
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, cholesky, buffer);
    for(i=0; i<f; i++) newPart[i] = target[i] + buffer->data[i];
  } while (!prior(newPart));
  
  return;
}

#define NZBIN 1
#define NNZ 5
cosmo_hm *initialize_cosmo_hm_ABC(double Omega_M, double sigma_8, double w0_de, error **err)
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
  int Nnz[NZBIN]           = {NNZ};
  double par_nz[NZBIN*NNZ] = {0.0, 3.0, 2.0, 1.0, 0.5};
  nofz_t nofz[NZBIN]       = {ludo};
  cosmo_hm *cmhm = init_parameters_hm(Omega_M, 1.0 - Omega_M, w0_de, 0.0, NULL, 0, 
				      0.78, 0.047, 0.0, 0.0, sigma_8, 0.95,
				      NZBIN, Nnz, nofz, par_nz, -1, -1, 
				      smith03_revised, eisenhu_osc, growth_de, linder, norm_s8,
				      11.0, 1.0, 0.13, j01, halo_bias_sc,
				      1.0e13, 1.0e11, 1.0e13, 0.3, 0.5,
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				      0.0, 0.0, -1, -1, 1.0,
				      0.0, 0.0,
				      berwein02_hexcl, 60.0, err);
  forwardError(*err, __LINE__,);
  return cmhm;
}
#undef NZBIN
#undef NNZ

void generateModel(peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err)
{
  cosmo_hm *cmhmABC = initialize_cosmo_hm_ABC(newPart[0], newPart[1], newPart[2], err); forwardError(*err, __LINE__,);
  setMassSamplers(cmhmABC, peak, ABC->sampArr, err);                                    forwardError(*err, __LINE__,);
  if (peak->doRandGalPos == 0) updateCosmo_gal_map(cmhmABC, peak, ABC->gMap, err);      forwardError(*err, __LINE__,);
  
  multiscaleFromMassFct(cmhmABC, peak, ABC->sampArr, ABC->hMap, ABC->galSamp, ABC->gMap, ABC->CCDMask,
			ABC->smoother, ABC->kMap, ABC->variance, ABC->peakList, ABC->hist, ABC->hist2, ABC->multiscale, err); forwardError(*err, __LINE__,);
  
  ABC->summFct(ABC->peakList, ABC->hist, ABC->multiscale, ABC->x_mod->array);
  int f2 = ABC->f + 2;
  newPart[f2-2] = ABC->distFct(ABC->x_obs->array, ABC->x_mod->array);
  free_parameters_hm(&cmhmABC);
  return;
}

void acceptParticleFromPrior(peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err)
{
  //-- WARNING: f-dependent
  
  int f2 = ABC->f + 2;
  int reject = 1;
  do {
    do priorGenerator(peak->generator, newPart);
    while (!ABC->priorFct(newPart));
    if (peak->printMode == 2) printf("(Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), ", newPart[0], newPart[1], newPart[2]);
    
    generateModel(peak, ABC, newPart, err);
  
    if (isError(*err)) {
      if      (peak->printMode == 2) printf("get error and resample\n");
      else if (peak->printMode == 3) printf("Proc %2d: (Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), get error and resample\n", peak->MPIInd, newPart[0], newPart[1], newPart[2]);
      purgeError(err);
    }
    
    else {
      ABC->newIter->nbAttempts += 1;
      reject = 0;
      if      (peak->printMode == 2) printf("delta = %7.3f, accepted\n", newPart[f2-2]);
      else if (peak->printMode == 3) printf("Proc %2d: (Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, accepted\n", 
					    peak->MPIInd, newPart[0], newPart[1], newPart[2], ABC->hMap->total, ABC->gMap->total, ABC->gMap->kappa_mean, ABC->gMap->fillingThreshold, newPart[f2-2]);
    }
  }
  while (reject);
  return;
}

void acceptParticle(peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err)
{
  //-- WARNING: f-dependent
  
  int f2 = ABC->f + 2;
  int reject = 1;
  do {
    generateParam(peak, ABC->oldIter, newPart, ABC->priorFct);
    if (peak->printMode == 2) printf("(Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), ", newPart[0], newPart[1], newPart[2]);
    
    generateModel(peak, ABC, newPart, err);
    
    if (isError(*err)) {
      if      (peak->printMode == 2) printf("get error and resample\n");
      else if (peak->printMode == 3) printf("Proc %2d: (Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), get error and resample\n", peak->MPIInd, newPart[0], newPart[1], newPart[2]);
      purgeError(err);
    }
    
    else {
      ABC->newIter->nbAttempts += 1;
      
      if (newPart[f2-2] <= ABC->newIter->epsilon) {
	reject = 0;
	if      (peak->printMode == 2) printf("delta = %6.3f, accepted\n", newPart[f2-2]);
	else if (peak->printMode == 3) printf("Proc %2d: (Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, accepted\n", 
					      peak->MPIInd, newPart[0], newPart[1], newPart[2], ABC->hMap->total, ABC->gMap->total, ABC->gMap->kappa_mean, ABC->gMap->fillingThreshold, newPart[f2-2]);
      }
      else {
	if      (peak->printMode == 2) printf("delta = %6.3f, rejected\n", newPart[f2-2]);
	else if (peak->printMode == 3) printf("Proc %2d: (Omega_M, sigma_8, w0_de) = (%.3f, %.3f, % .3f), %6d halos, %d galaxies, kappa mean = %.5f, threshold = %.2f, delta = %7.3f, rejected\n", 
					      peak->MPIInd, newPart[0], newPart[1], newPart[2], ABC->hMap->total, ABC->gMap->total, ABC->gMap->kappa_mean, ABC->gMap->fillingThreshold, newPart[f2-2]);
      }
    }
  }
  while (reject);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to loop

void swapOldAndNew(PMC_ABC_t *ABC)
{
  iteration_t *buffer = ABC->oldIter;
  ABC->oldIter = ABC->newIter;
  ABC->newIter = buffer;
  return;
}

double argOfExp(double *oldPa, double *newPart, gsl_matrix *invCov, gsl_vector *Delta_param, gsl_vector *intermediate)
{
  double value;
  int i;
  for (i=0; i<Delta_param->size; i++) Delta_param->data[i] = newPart[i] - oldPa[i];
  gsl_blas_dsymv(CblasUpper, 1.0, invCov, Delta_param, 0.0, intermediate); //-- intermediate = invCov * Delta_param
  gsl_blas_ddot(Delta_param, intermediate, &value);
  value *= -0.5;
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
  
  gsl_matrix *invCov       = ABC->oldIter->invCov;
  gsl_vector *Delta_param  = ABC->oldIter->buffer;
  gsl_vector *intermediate = ABC->newIter->buffer;
  double sum = 0.0;
  
  double *oldPa, *newPart;
  double arg, weight;
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    weight = 0.0;
    for (oldPa=oldBegin; oldPa<oldEnd; oldPa+=f2) {
      arg = argOfExp(oldPa, newPart, invCov, Delta_param, intermediate);
      weight += oldPa[f2-1] * exp(arg);
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
  
  int f2 = ABC->f + 2;
  int Q  = ABC->Q;
  double *oldBegin = ABC->oldIter->matrix->matrix + f2 - 2;
  double *deltaArr = ABC->deltaList->array;
  
  double *oldPa;
  int i;
  
  for (i=0, oldPa=oldBegin; i<Q; i++, oldPa+=f2) deltaArr[i] = oldPa[0];
  gsl_sort(deltaArr, 1, Q);
  ABC->newIter->epsilon = gsl_stats_median_from_sorted_data(deltaArr, 1, Q);
  return;
}

void loopZero(peak_param *peak, PMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  ABC->t = 0;
  ABC->newIter->nbAttempts = 0;
  
  if (peak->MPIInd == 0) {
    printf("-- Begin loop %d\n\n", ABC->t);
    printf("Tolerance = DBL_MAX\n\n");
  }
  MPI_Barrier(MPI_COMM_WORLD); //-- Wait for print
  
  int f2    = ABC->f + 2;
  int Q     = ABC->Q;
  int Q_MPI = ABC->Q_MPI;
  double weight    = 1.0 / (double)Q;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *newEnd   = newBegin + f2 * Q_MPI;
  double *newPart;
  
  //-- Size ajustment for MPI
  if (peak->MPIInd == peak->MPISize - 1) newEnd -= f2 * (peak->MPISize * Q_MPI - Q);
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    acceptParticleFromPrior(peak, ABC, newPart, err); forwardError(*err, __LINE__,);
    newPart[f2-1] = weight;
  }
  
  printf("Proc %2d finished\n", peak->MPIInd);
  
  MPI_Allgather(newBegin, f2 * Q_MPI, MPI_DOUBLE, newBegin, f2 * Q_MPI, MPI_DOUBLE, MPI_COMM_WORLD);        //-- Communicate new particles
  MPI_Allreduce(&ABC->newIter->nbAttempts, &ABC->newIter->nbAttempts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //-- Communicate the number of attempts
  
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
  
  if (peak->MPIInd == 0) {
    printf("\nSuccess rate = %.5f\n", Q / (double)ABC->newIter->nbAttempts);
    printf("\n-- End loop 0, "); routineTime(start, clock());
    printf("------------------------------------------------------------------------\n");
  }
  return;
}

void loop(peak_param *peak, PMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  swapOldAndNew(ABC);
  setEpsilon(ABC);
  ABC->t++;
  ABC->newIter->nbAttempts = 0;
  
  if (peak->MPIInd == 0) {
    printf("-- Begin loop %d\n\n", ABC->t);
    printf("Tolerance = %.5f\n\n", ABC->newIter->epsilon);
  }
  MPI_Barrier(MPI_COMM_WORLD); //-- Wait for print
  
  int f2    = ABC->f + 2;
  int Q     = ABC->Q;
  int Q_MPI = ABC->Q_MPI;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *newEnd   = newBegin + f2 * Q_MPI;
  double *newPart;
  
  //-- Size ajustment for MPI
  if (peak->MPIInd == peak->MPISize - 1) newEnd -= f2 * (peak->MPISize * Q_MPI - Q);
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    acceptParticle(peak, ABC, newPart, err);
    forwardError(*err, __LINE__,);
  }
  
  printf("Proc %2d finished\n", peak->MPIInd);
  
  MPI_Allgather(newBegin, f2 * Q_MPI, MPI_DOUBLE, newBegin, f2 * Q_MPI, MPI_DOUBLE, MPI_COMM_WORLD);        //-- Communicate new particles
  MPI_Allreduce(&ABC->newIter->nbAttempts, &ABC->newIter->nbAttempts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //-- Communicate the number of attempts
  
  setWeights(ABC);
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
  
  if (peak->MPIInd == 0) {
    printf("\nSuccess rate = %.5f\n", Q / (double)ABC->newIter->nbAttempts);
    printf("\n-- End loop %d, ", ABC->t); routineTime(start, clock());
    printf("------------------------------------------------------------------------\n");
  }
  return;
}

void readAndLoop(char name[], peak_param *peak, PMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  read_iteration_t(name, ABC->oldIter, ABC->deltaList->array, err); forwardError(*err, __LINE__,);
  setEpsilon(ABC);
  ABC->t++;
  
  if (peak->MPIInd == 0) {
    printf("-- Begin loop %d\n\n", ABC->t);
    printf("Tolerance = %.5f\n\n", ABC->newIter->epsilon);
  }
  MPI_Barrier(MPI_COMM_WORLD); //-- Wait for print
  
  int f2    = ABC->f + 2;
  int Q     = ABC->Q;
  int Q_MPI = ABC->Q_MPI;
  double *newBegin = ABC->newIter->matrix->matrix;
  double *newEnd   = newBegin + f2 * Q_MPI;
  double *newPart;
  
  //-- Size ajustment for MPI
  if (peak->MPIInd == peak->MPISize - 1) newEnd -= f2 * (peak->MPISize * Q_MPI - Q);
  
  for (newPart=newBegin; newPart<newEnd; newPart+=f2) {
    acceptParticle(peak, ABC, newPart, err);
    forwardError(*err, __LINE__,);
  }
  
  MPI_Allgather(newBegin, f2 * Q_MPI, MPI_DOUBLE, newBegin, f2 * Q_MPI, MPI_DOUBLE, MPI_COMM_WORLD);        //-- Communicate new particles
  MPI_Allreduce(&ABC->newIter->nbAttempts, &ABC->newIter->nbAttempts, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); //-- Communicate the number of attempts
  
  setWeights(ABC);
  updateMean_iteration_t(ABC->newIter);
  updateCovariance_iteration_t(ABC->newIter);
  updateCholesky_iteration_t(ABC->newIter);
    
  if (peak->MPIInd == 0) {
    printf("\nSuccess rate = %.5f\n", Q / (double)ABC->newIter->nbAttempts);
    printf("\n-- End loop %d, ", ABC->t); routineTime(start, clock());
    printf("------------------------------------------------------------------------\n");
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to prior

#define OMEGA_M_MIN 0.10
#define OMEGA_M_MAX 0.90
#define SIGMA_8_MIN 0.30
#define SIGMA_8_MAX 1.60
#define W0_DE_MIN -1.8
#define W0_DE_MAX 0.0
void priorGenerator(gsl_rng *generator, double *part)
{
  part[0] = gsl_ran_flat(generator, OMEGA_M_MIN, OMEGA_M_MAX);
  part[1] = gsl_ran_flat(generator, SIGMA_8_MIN, SIGMA_8_MAX);
  part[2] = gsl_ran_flat(generator, W0_DE_MIN, W0_DE_MAX);
  return;
}

int prior_rectangle(double *part)
{
  int boolean = (part[0] >= OMEGA_M_MIN) && (part[0] < OMEGA_M_MAX)
             && (part[1] >= SIGMA_8_MIN) && (part[1] < SIGMA_8_MAX)
             && (part[2] >= W0_DE_MIN)   && (part[2] < W0_DE_MAX);
  return boolean;
}

int prior_pentagon(double *part)
{
  int boolean = (part[0] + part[1] <= 2.0)
             && (part[0] >= OMEGA_M_MIN) && (part[0] < OMEGA_M_MAX)
             && (part[1] >= SIGMA_8_MIN) && (part[1] < SIGMA_8_MAX)
             && (part[2] >= W0_DE_MIN)   && (part[2] < W0_DE_MAX);
  return boolean;
}
#undef OMEGA_M_MIN
#undef OMEGA_M_MAX
#undef SIGMA_8_MIN
#undef SIGMA_8_MAX
#undef W0_DE_MIN
#undef W0_DE_MAX

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
  PMC_ABC_t *ABC = initialize_PMC_ABC_t(peak, err); forwardError(*err, __LINE__,);
  fillObservation("../demo/x_obs", peak, ABC, err); forwardError(*err, __LINE__,);
  if (peak->MPIInd == 0) printf("------------------------------------------------------------------------\n");
  
  char name[STRING_LENGTH_MAX];
  FILE *file;
  
  //-- Loop 0
  loopZero(peak, ABC, err); forwardError(*err, __LINE__,);
  if (peak->MPIInd == 0) {
    sprintf(name, "iteration_%s_p%d_t%d", STR_SUMMARY_T(ABC->summ), ABC->Q, ABC->t);
    file = fopen(name, "w");
    output2_iteration_t(file, ABC->newIter);
    fclose(file);
  }
  MPI_Barrier(MPI_COMM_WORLD); //-- Wait for output
  
  //-- Other loops
  while (ABC->newIter->nbAttempts * ABC->r_stop <= ABC->Q) {
    loop(peak, ABC, err); forwardError(*err, __LINE__,); //-- Switch old and new inside
    if (peak->MPIInd == 0) {
      sprintf(name, "iteration_%s_p%d_t%d", STR_SUMMARY_T(ABC->summ), ABC->Q, ABC->t);
      file = fopen(name, "w");
      output2_iteration_t(file, ABC->newIter);
      fclose(file);
    }
    MPI_Barrier(MPI_COMM_WORLD); //-- Wait for output
  }
  
  if (peak->MPIInd == 0) printf("%d iterations done\n", ABC->t+1);
  free_PMC_ABC_t(ABC);
  return;
}

//----------------------------------------------------------------------

