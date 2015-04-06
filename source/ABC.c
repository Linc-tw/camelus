

  /***********************************************
   **  ABC.c					**
   **  Chieh-An Lin				**
   **  Version 2015.04.06			**
   **						**
   **  References:				**
   **  - Marin et al. (2011)			**
   **  - Weyant et al. (2013) - ApJ, 764, 116	**
   ***********************************************/


#include "ABC.h"


//----------------------------------------------------------------------
//-- Functions related to particle_t

particle_t *initialize_particle_t(int d, error **err)
{
  particle_t *pa = (particle_t*)malloc_err(sizeof(particle_t), err); forwardError(*err, __LINE__,);
  pa->d          = d;
  pa->diff       = 0.0;
  pa->weight     = 0.0;
  pa->param      = (double*)malloc_err(d * sizeof(double), err);     forwardError(*err, __LINE__,);
  pa->param_gsl  = gsl_vector_alloc(d);
  return pa;
}

void free_particle_t(particle_t *pa)
{
  if (pa->param) {free(pa->param); pa->param = NULL;}
  gsl_vector_free(pa->param_gsl);
  free(pa);
  return;
}

void print_particle_t(particle_t *pa)
{
  //-- WARNING: d-dependent
  printf("(Omega_M, sigma_8) = (%.3f, %.3f), weight = %.3f, diff = %9.5f\n", pa->param[0], pa->param[1], pa->weight, pa->diff);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to particle_arr

particle_arr *initialize_particle_arr(int d, int p, error **err)
{
  particle_arr *part = (particle_arr*)malloc_err(sizeof(particle_arr), err);   forwardError(*err, __LINE__,);
  part->d            = d;
  part->p            = p;
  part->nbAttempts   = 0;
  part->array        = (particle_t**)malloc_err(p * sizeof(particle_t*), err); forwardError(*err, __LINE__,);
  part->mean         = initialize_double_arr(d);
  part->cov          = gsl_matrix_alloc(d, d);
  part->cov2         = gsl_matrix_alloc(d, d);
  part->invCov       = gsl_matrix_alloc(d, d);
  part->perm         = gsl_permutation_alloc(d);
  
  int i;
  for (i=0; i<p; i++) {
    part->array[i] = initialize_particle_t(d, err);
    forwardError(*err, __LINE__,);
  }
  return part;
}

void free_particle_arr(particle_arr *part)
{
  int i;
  if (part->array) {
    for (i=0; i<part->p; i++) {free_particle_t(part->array[i]); part->array[i] = NULL;}
  }
  if (part->mean) {free_double_arr(part->mean); part->mean = NULL;}
  gsl_matrix_free(part->cov);
  gsl_matrix_free(part->cov2);
  gsl_matrix_free(part->invCov);
  gsl_permutation_free(part->perm);
  free(part);
  return;
}

void print_particle_arr(particle_arr *part)
{
  printf("# particle_arr (first 20 elements)\n");
  int L = MIN(part->p, 20);
  
  int i;
  for (i=0; i<L; i++) print_particle_t(part->array[i]);
  return;
}

void output1_particle_arr(FILE *file, particle_arr *part)
{
  //-- WARNING: d-dependent
  int p = part->p;
  fprintf(file, "# particle_arr\n");
  fprintf(file, "# p                  = %d\n", p);
  fprintf(file, "# Number of attempts = %d\n", part->nbAttempts);
  fprintf(file, "# Success rate       = %.5f\n", p / (double)part->nbAttempts);
  fprintf(file, "#\n");
  
  particle_t **partArr = part->array;
  int i; 
  
  fprintf(file, "Omega_M = [%.5f", partArr[0]->param[0]);
  for (i=1; i<p; i++) fprintf(file, ", %.5f", partArr[i]->param[0]);
  fprintf(file, "]\n");
  
  fprintf(file, "sigma_8 = [%.5f", partArr[0]->param[1]);
  for (i=1; i<p; i++) fprintf(file, ", %.5f", partArr[i]->param[1]);
  fprintf(file, "]\n");
  
  fprintf(file, "weight  = [%.5f", partArr[0]->weight);
  for (i=1; i<p; i++) fprintf(file, ", %.5f", partArr[i]->weight);
  fprintf(file, "]\n");
  
  fprintf(file, "diff    = [%.5f", partArr[0]->diff);
  for (i=1; i<p; i++) fprintf(file, ", %.5f", partArr[i]->diff);
  fprintf(file, "]\n");
  fprintf(file, "\n");
  
  fprintf(file, "cov    = [[% .5f, % .5f], [% .5f, % .5f]]\n",
          gsl_matrix_get(part->cov, 0, 0), gsl_matrix_get(part->cov, 0, 1),
          gsl_matrix_get(part->cov, 1, 0), gsl_matrix_get(part->cov, 1, 1));
  
  fprintf(file, "invCov = [[% .5f, % .5f], [% .5f, % .5f]]\n",
          gsl_matrix_get(part->invCov, 0, 0), gsl_matrix_get(part->invCov, 0, 1),
          gsl_matrix_get(part->invCov, 1, 0), gsl_matrix_get(part->invCov, 1, 1));
  return;
}

void output2_particle_arr(FILE *file, particle_arr *part)
{
  //-- WARNING: d-dependent
  int p = part->p;
  fprintf(file, "# particle_arr\n");
  fprintf(file, "# p                  = %d\n", p);
  fprintf(file, "# Number of attempts = %d\n", part->nbAttempts);
  fprintf(file, "# Success rate       = %.5f\n", p / (double)part->nbAttempts);
  fprintf(file, "# epsilon            = %.5f\n", part->epsilon);
  fprintf(file, "# cov                = [[% .5f, % .5f], [% .5f, % .5f]]\n",
          gsl_matrix_get(part->cov, 0, 0), gsl_matrix_get(part->cov, 0, 1),
          gsl_matrix_get(part->cov, 1, 0), gsl_matrix_get(part->cov, 1, 1));
  fprintf(file, "# invCov             = [[% .5f, % .5f], [% .5f, % .5f]]\n",
          gsl_matrix_get(part->invCov, 0, 0), gsl_matrix_get(part->invCov, 0, 1),
          gsl_matrix_get(part->invCov, 1, 0), gsl_matrix_get(part->invCov, 1, 1));
  fprintf(file, "#\n");
  fprintf(file, "# Omega_M  sigma_8   weight       diff\n");
  
  particle_t *pa;
  int i;
  for (i=0; i<p; i++) {
    pa = part->array[i];
    fprintf(file, "  %.5f  %.5f  %.5f  %9.5f\n", pa->param[0], pa->param[1], pa->weight, pa->diff);
  }
  return;
}

void updateEpsilon_particle_arr(particle_arr *part, double *diffArr)
{
  int p = part->p;
  int i;
  for (i=0; i<p; i++) diffArr[i] = part->array[i]->diff;
  
  gsl_sort(diffArr, 1, p);
  part->epsilon = gsl_stats_median_from_sorted_data(diffArr, 1, p);
  printf("epsilon_median = %.5f\n", part->epsilon);
  printf("Success rate   = %.5f\n", p / (double)part->nbAttempts);
  return;
}

void updateMean_particle_arr(particle_arr *part)
{
  int d = part->d;
  int p = part->p;
  double sum;
  
  int i, j;
  for (j=0; j<d; j++) {
    sum = 0.0;
    for (i=0; i<p; i++) sum += part->array[i]->weight * part->array[i]->param[j];
    part->mean->array[j] = sum;
  }
  return;
}

void updateCovariance_particle_arr(particle_arr *part)
{
  int d = part->d;
  int p = part->p;
  gsl_matrix *cov = part->cov;
  double *meanArr = part->mean->array;
  double w_tot    = 0.0;
  double w_sq_tot = 0.0;
  double mean_i, mean_j, sum;
  particle_t *pa;
  
  int k;
  for (k=0; k<p; k++) {
    w_tot    += part->array[k]->weight;
    w_sq_tot += pow(part->array[k]->weight, 2);
  }
  
  int i, j;
  for (j=0; j<d; j++) {
    for (i=0; i<d; i++) {
      sum = 0.0;
      if (i < j) {
	gsl_matrix_set(cov, i, j, gsl_matrix_get(cov, j, i));
	continue;
      }
      mean_i = meanArr[i];
      mean_j = meanArr[j];
      for (k=0; k<p; k++) {
	pa   = part->array[k];
	sum += pa->weight * (pa->param[i] - mean_i) * (pa->param[j] - mean_j);
      }
      sum *= w_tot / (pow(w_tot, 2) - w_sq_tot);
      gsl_matrix_set(cov, i, j, sum);
    }
  }
  
  //-- Make a copy, because the invertion will destroy cov
  gsl_matrix_memcpy(part->cov2, cov);
  
  //-- Make LU decomposition and invert
  int s;
  gsl_linalg_LU_decomp(part->cov2, part->perm, &s);
  gsl_linalg_LU_invert(part->cov2, part->perm, part->invCov);
  
  //-- Debias the C_{ij}^{-1}
  //-- C^-1_unbiased = (p - d - 2) / (p - 1) * C^-1_biased
  double factor = (p - d -2) / (double)(p - 1);
  gsl_matrix_scale(part->invCov, factor);
  return;
}

void read_particle_arr(char name[], particle_arr *part, double *diffArr, error **err)
{
  //-- WARNING: d-dependent
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  particle_t *pa;
  double pos[2], w, z, M;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      testError(count==part->p, peak_overflow, "Too many particles", *err, __LINE__); forwardError(*err, __LINE__,);
      ungetc(c, file);
      pa = part->array[count];
      buffer2 = fscanf(file, "%lf %lf %lf %lf\n", &pa->param[0], &pa->param[1], &pa->weight, &pa->diff);
      count++;
    }
    c = fgetc(file);
  }
  
  fclose(file);
  printf("\"%s\" read\n", name);
  printf("%d particles generated\n", count);
  updateEpsilon_particle_arr(part, diffArr);
  updateMean_particle_arr(part);
  updateCovariance_particle_arr(part);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to SMC_ABC_t

SMC_ABC_t *initialize_SMC_ABC_t(peak_param *peak, error **err)
{
  SMC_ABC_t *ABC = (SMC_ABC_t*)malloc_err(sizeof(SMC_ABC_t), err); forwardError(*err, __LINE__,);
  ABC->d         = peak->ABC_d;
  ABC->p         = peak->ABC_p;
  ABC->r_stop    = peak->ABC_r_stop;
  int i;
  STRING2ENUM(ABC->summ, peak->ABC_summ, summary_t, STR_SUMMARY_T, i, NB_SUMMARY_T, err); forwardError(*err, __LINE__,);
  
  ABC->epsilon_0 = DBL_MAX;
  ABC->priorFct  = prior_pentagon;
  
  ABC->t         = 0;
  ABC->oldPart   = initialize_particle_arr(ABC->d, ABC->p, err); forwardError(*err, __LINE__,);
  ABC->newPart   = initialize_particle_arr(ABC->d, ABC->p, err); forwardError(*err, __LINE__,);
  ABC->diffList  = initialize_double_arr(ABC->p);
  
  if (ABC->summ == abd6) {
    ABC->obsSummary   = initialize_double_arr(6);
    ABC->simulSummary = initialize_double_arr(6);
    ABC->peakHist     = initialize_hist_t(6);
    set_hist_t(ABC->peakHist, 3.5, 6.5);
    ABC->peakHist->x_max = 1000.0; //-- So that bins are [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 1000.0]
    ABC->summaryFct   = summary_abd_all;
    ABC->distFct      = dist_abd6;
  }
  else if (ABC->summ == pct6) {
    ABC->obsSummary   = initialize_double_arr(6);
    ABC->simulSummary = initialize_double_arr(6);
    ABC->summaryFct   = summary_pct6;
    ABC->distFct      = dist_6D;
  }
  else if (ABC->summ == cut6) {
    ABC->obsSummary   = initialize_double_arr(6);
    ABC->simulSummary = initialize_double_arr(6);
    ABC->summaryFct   = summary_cut6;
    ABC->distFct      = dist_6D;
  }
  else if (ABC->summ == abd5) {
    ABC->obsSummary   = initialize_double_arr(5);
    ABC->simulSummary = initialize_double_arr(5);
    ABC->peakHist     = initialize_hist_t(5);
    ABC->peakHist->x_lower[0] = 3.0;
    ABC->peakHist->x_lower[1] = 3.8;
    ABC->peakHist->x_lower[2] = 4.5;
    ABC->peakHist->x_lower[3] = 5.3;
    ABC->peakHist->x_lower[4] = 6.2;
    ABC->peakHist->x_max      = 1000.0; //-- So that bins are [3.0, 3.8, 4.5, 5.3, 6.2, 1000.0]
    ABC->summaryFct   = summary_abd_all;
    ABC->distFct      = dist_abd5;
  }
  else if (ABC->summ == pct5) {
    ABC->obsSummary   = initialize_double_arr(5);
    ABC->simulSummary = initialize_double_arr(5);
    ABC->summaryFct   = summary_pct5;
    ABC->distFct      = dist_pct5;
  }
  else if (ABC->summ == cut5) {
    ABC->obsSummary   = initialize_double_arr(5);
    ABC->simulSummary = initialize_double_arr(5);
    ABC->summaryFct   = summary_cut5;
    ABC->distFct      = dist_cut5;
  }
  else if (ABC->summ == pct4) {
    ABC->obsSummary   = initialize_double_arr(4);
    ABC->simulSummary = initialize_double_arr(4);
    ABC->summaryFct   = summary_pct4;
    ABC->distFct      = dist_4D;
  }
  else if (ABC->summ == pct998) {
    ABC->obsSummary   = initialize_double_arr(1);
    ABC->simulSummary = initialize_double_arr(1);
    ABC->summaryFct   = summary_pct998;
    ABC->distFct      = dist_1D;
  }
  else if (ABC->summ == pct996) {
    ABC->obsSummary   = initialize_double_arr(1);
    ABC->simulSummary = initialize_double_arr(1);
    ABC->summaryFct   = summary_pct996;
    ABC->distFct      = dist_1D;
  }
  else if (ABC->summ == cut900) {
    ABC->obsSummary   = initialize_double_arr(1);
    ABC->simulSummary = initialize_double_arr(1);
    ABC->summaryFct   = summary_cut900;
    ABC->distFct      = dist_1D;
  }
  else {*err = addError(peak_unknown, "Unknown summary type", *err, __LINE__); forwardError(*err, __LINE__,);}
  
  int length       = (peak->resol[0] - 2 * peak->bufferSize) * (peak->resol[1] - 2 * peak->bufferSize);
  cosmo_hm *cmhm   = initialize_cosmo_hm_default(err);                                          forwardError(*err, __LINE__,);
  
  ABC->sampArr     = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
  ABC->hMap        = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  ABC->gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  makeGalaxies(cmhm, peak, ABC->gMap, err);                                                     forwardError(*err, __LINE__,);
  ABC->kMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  ABC->nMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  ABC->transformer = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
  fillGaussianKernel(ABC->transformer, peak->s);
  ABC->peakList    = initialize_double_arr(length);
  
  free_parameters_hm(&cmhm);
  printf("ABC initialization done\n");
  return ABC;
}

void free_SMC_ABC_t(SMC_ABC_t *ABC)
{
  if (ABC->oldPart)      {free_particle_arr(ABC->oldPart);    ABC->oldPart = NULL;}
  if (ABC->newPart)      {free_particle_arr(ABC->newPart);    ABC->newPart = NULL;}
  if (ABC->diffList)     {free_double_arr(ABC->diffList);     ABC->diffList = NULL;}
  
  if (ABC->obsSummary)   {free_double_arr(ABC->obsSummary);   ABC->obsSummary = NULL;}
  if (ABC->simulSummary) {free_double_arr(ABC->simulSummary); ABC->simulSummary = NULL;}
  if (ABC->peakHist)     {free_hist_t(ABC->peakHist);         ABC->peakHist = NULL;}
  
  if (ABC->sampArr)      {free_sampler_arr(ABC->sampArr);     ABC->sampArr = NULL;}
  if (ABC->hMap)         {free_halo_map(ABC->hMap);           ABC->hMap = NULL;}
  if (ABC->gMap)         {free_gal_map(ABC->gMap);            ABC->gMap = NULL;}
  if (ABC->kMap)         {free_map_t(ABC->kMap);              ABC->kMap = NULL;}
  if (ABC->nMap)         {free_map_t(ABC->nMap);              ABC->nMap = NULL;}
  if (ABC->transformer)  {free_FFT_t(ABC->transformer);       ABC->transformer = NULL;}
  if (ABC->peakList)     {free_double_arr(ABC->peakList);     ABC->peakList = NULL;}
  free(ABC); ABC = NULL;
  return;
}

void fillObservation_SMC_ABC_t(char name[], SMC_ABC_t *ABC, error **err)
{
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double *array = ABC->peakList->array;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf\n", &array[count]);
      count++;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  ABC->peakList->length = count;
  ABC->summaryFct(ABC->peakList, ABC->peakHist, ABC->obsSummary->array);
  printf("Observation data read\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to accept-reject algorithm

void generateParam(peak_param *peak, particle_arr *oldPart, particle_t *newPa, prior_fct *prior)
{
  //-- WARNING: d-dependent
  if (newPa->d != 2) {
    printf("Need to implement in ABC.c/generateParam\n");
    exit(1);
  }
  
  gsl_rng *generator      = peak->generator;
  gsl_matrix *cov         = oldPart->cov;
  particle_t **oldPartArr = oldPart->array;
  double p, x1, x2, var2;
  
  particle_t *target;
  int i;
  do {
    i = 0;
    target = oldPartArr[0];
    p = gsl_ran_flat(generator, 0.0, 1.0);
    while (p > target->weight) {
      p -= target->weight;
      i++;
      target = oldPartArr[i];
    }
    
    //-- Compute pdf from multivariate Gaussian pdf
    x1   = gsl_ran_gaussian(generator, sqrt(gsl_matrix_get(cov, 0, 0)));
    var2 = gsl_matrix_get(cov, 1, 1) - SQ(gsl_matrix_get(cov, 1, 0)) / gsl_matrix_get(cov, 0, 0);
    x2   = gsl_ran_gaussian(generator, sqrt(var2));
    x2  += target->param[1] + x1 * gsl_matrix_get(cov, 1, 0) / gsl_matrix_get(cov, 0, 0);
    x1  += target->param[0];
    newPa->param[0] = x1;
    newPa->param[1] = x2;
  } while (!prior(newPa));
  
  return;
}

#define NZBIN 1
#define NNZ 5
cosmo_hm *initialize_cosmo_hm_ABC(double Omega_M, double sigma_8, error **err)
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
  cosmo_hm *cmhm = init_parameters_hm(Omega_M, 1-Omega_M, -1.0, 0.0, NULL, 0, 
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

void generateObs(peak_param *peak, SMC_ABC_t *ABC, particle_t *newPa, error **err)
{
  cosmo_hm *cmhmABC = initialize_cosmo_hm_ABC(newPa->param[0], newPa->param[1], err); forwardError(*err, __LINE__,);
  updateCosmo_gal_map(cmhmABC, peak, ABC->gMap, err);                                 forwardError(*err, __LINE__,);
  setMassSamplers(cmhmABC, peak, ABC->sampArr, err);                                  forwardError(*err, __LINE__,);
  peakListFromMassFct(cmhmABC, peak, ABC->sampArr, ABC->hMap, ABC->gMap,
		      ABC->kMap, ABC->nMap, ABC->transformer, ABC->peakList, err);    forwardError(*err, __LINE__,);
  
  ABC->summaryFct(ABC->peakList, ABC->peakHist, ABC->simulSummary->array);
  newPa->diff = ABC->distFct(ABC->obsSummary->array, ABC->simulSummary->array);
  free_parameters_hm(&cmhmABC);
  return;
}

void acceptParticleFromPrior(peak_param *peak, SMC_ABC_t *ABC, particle_t *newPa, error **err)
{
  //-- WARNING: d-dependent
  int reject = 1;
  do {
    do priorGenerator(peak->generator, newPa);
    while (!ABC->priorFct(newPa));
    printf("(Omega_M, sigma_8) = (%.3f, %.3f), ", newPa->param[0], newPa->param[1]);

    generateObs(peak, ABC, newPa, err);
    if (isError(*err)) {
      printf(", get error and resample\n");
      purgeError(err);
    }
    else {
      ABC->newPart->nbAttempts += 1;
      printf("diff = %6.3f, ", newPa->diff);
      reject = 0;
      printf("accepted\n");
    }
  }
  while (reject);
  return;
}

void acceptParticle(peak_param *peak, SMC_ABC_t *ABC, particle_t *newPa, error **err)
{
  //-- WARNING: d-dependent
  int reject = 1;
  do {
    generateParam(peak, ABC->oldPart, newPa, ABC->priorFct);
    printf("(Omega_M, sigma_8) = (%.3f, %.3f), ", newPa->param[0], newPa->param[1]);
    
    generateObs(peak, ABC, newPa, err);
    if (isError(*err)) {
      printf(", get error and resample\n");
      purgeError(err);
    }
    else {
      ABC->newPart->nbAttempts += 1;
      printf("diff = %6.3f, ", newPa->diff);
      if (newPa->diff <= ABC->oldPart->epsilon) {
	reject = 0;
	printf("accepted\n");
      }
      else printf("rejected\n");
    }
  }
  while (reject);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to loop

void swapOldAndNew(SMC_ABC_t *ABC)
{
  particle_arr *buffer = ABC->oldPart;
  ABC->oldPart = ABC->newPart;
  ABC->newPart = buffer;
  return;
}

double argOfExp(particle_t *oldPa, particle_t *newPa, gsl_matrix *invCov)
{
  gsl_vector *Delta_param  = oldPa->param_gsl;
  gsl_vector *intermediate = newPa->param_gsl;
  double value;
  
  int i;
  for (i=0; i<newPa->d; i++) gsl_vector_set(Delta_param, i, newPa->param[i] - oldPa->param[i]);
  gsl_blas_dsymv(CblasUpper, 1.0, invCov, Delta_param, 0.0, intermediate); //-- intermediate = invCov * Delta_param
  gsl_blas_ddot(Delta_param, intermediate, &value);
  value *= -0.5;
  return value;
}

void setWeights(SMC_ABC_t *ABC)
{
  particle_t **oldArr = ABC->oldPart->array;
  particle_t **newArr = ABC->newPart->array;
  gsl_matrix *invCov  = ABC->oldPart->invCov;
  double arg, weight, sum = 0.0;
  int p = ABC->p;
  
  particle_t *oldPa, *newPa;
  int i, j;
  for (j=0; j<p; j++) {
    newPa = newArr[j];
    weight = 0.0;
    for (i=0; i<p; i++) {
      oldPa = oldArr[i];
      arg = argOfExp(oldPa, newPa, invCov);
      weight += oldPa->weight * exp(arg);
    }
    weight = 1.0 / weight;
    newPa->weight = weight;
    sum += weight;
  }
  
  for (i=0; i<p; i++) newArr[i]->weight /= sum;
  return;
}

void loopZero(peak_param *peak, SMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  printf("-- Begin loop 0\n\n");
  int p = ABC->p;
  particle_arr *newPart = ABC->newPart;
  
  int i;
  for (i=0; i<p; i++) {
    acceptParticleFromPrior(peak, ABC, newPart->array[i], err); forwardError(*err, __LINE__,);
    newPart->array[i]->weight = 1.0 / (double)p;
  }
  printf("\n");
  
  print_particle_arr(newPart);
  updateEpsilon_particle_arr(newPart, ABC->diffList->array);
  updateMean_particle_arr(newPart);
  updateCovariance_particle_arr(newPart);
  
  printf("\n-- End loop 0, "); routineTime(start, clock());
  printf("------------------------------------------------------------------------\n");
  return;
}

void loop(peak_param *peak, SMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  swapOldAndNew(ABC);
  ABC->t++;
  ABC->newPart->nbAttempts = 0;
  particle_arr *newPart    = ABC->newPart;
  particle_t **newPartArr  = newPart->array;
  printf("-- Begin loop %d\n\n", ABC->t);
  printf("Tolerance = %.5f\n", ABC->oldPart->epsilon);
  
  int i;
  for (i=0; i<ABC->p; i++) {
    acceptParticle(peak, ABC, newPartArr[i], err);
    forwardError(*err, __LINE__,);
  }
  setWeights(ABC);
  printf("\n");
  
  print_particle_arr(newPart);
  updateEpsilon_particle_arr(newPart, ABC->diffList->array);
  updateMean_particle_arr(newPart);
  updateCovariance_particle_arr(newPart);
  
  printf("\n-- End loop %d, ", ABC->t); routineTime(start, clock());
  printf("------------------------------------------------------------------------\n");
  return;
}

void readAndLoop(char name[], peak_param *peak, SMC_ABC_t *ABC, error **err)
{
  clock_t start = clock();
  ABC->t++;
  read_particle_arr(name, ABC->oldPart, ABC->diffList->array, err); forwardError(*err, __LINE__,);
  particle_arr *newPart   = ABC->newPart;
  particle_t **newPartArr = newPart->array;
  printf("-- Begin loop %d\n\n", ABC->t);
  printf("Tolerance = %.5f\n", ABC->oldPart->epsilon);
  
  int i;
  for (i=0; i<ABC->p; i++) {
    acceptParticle(peak, ABC, newPartArr[i], err);
    forwardError(*err, __LINE__,);
  }
  setWeights(ABC);
  printf("\n");
  
  print_particle_arr(newPart);
  updateEpsilon_particle_arr(newPart, ABC->diffList->array);
  updateMean_particle_arr(newPart);
  updateCovariance_particle_arr(newPart);
  
  printf("\n-- End loop %d, ", ABC->t); routineTime(start, clock());
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to prior

#define OMEGA_M_MIN 0.05
#define OMEGA_M_MAX 0.95
#define SIGMA_8_MIN 0.05
#define SIGMA_8_MAX 1.5
void priorGenerator(gsl_rng *generator, particle_t *pa)
{
  pa->param[0] = gsl_ran_flat(generator, OMEGA_M_MIN, OMEGA_M_MAX);
  pa->param[1] = gsl_ran_flat(generator, SIGMA_8_MIN, SIGMA_8_MAX);
  return;
}

int prior_rectangle(particle_t *pa)
{
  int boolean = (pa->param[0] >= OMEGA_M_MIN) && (pa->param[0] < OMEGA_M_MAX)
             && (pa->param[1] >= SIGMA_8_MIN) && (pa->param[1] < SIGMA_8_MAX);
  return boolean;
}

int prior_pentagon(particle_t *pa)
{
  double Omega_M = pa->param[0];
  double sigma_8 = pa->param[1];
  int boolean = (Omega_M + sigma_8 <= 2.0)
             && (Omega_M >= OMEGA_M_MIN) && (Omega_M < OMEGA_M_MAX)
             && (sigma_8 >= SIGMA_8_MIN) && (sigma_8 < SIGMA_8_MAX);
  return boolean;
}
#undef OMEGA_M_MIN
#undef OMEGA_M_MAX
#undef SIGMA_8_MIN
#undef SIGMA_8_MAX

//----------------------------------------------------------------------
//-- Functions related to summary statistics

void summary_abd_all(double_arr *peakList, hist_t *hist, double *summ)
{
  reset_hist_t(hist);
  int i;
  for (i=0; i<peakList->length; i++) silentPush_hist_t(hist, peakList->array[i]);
  for (i=0; i<hist->length; i++) summ[i] = (double)hist->n[i];
  return;
}

void summary_pct6(double_arr *peakList, hist_t *hist, double *summ)
{
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.978);
  summ[1] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.988);
  summ[2] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.993);
  summ[3] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.996);
  summ[4] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.997);
  summ[5] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.998);
  return;
}

void summary_cut6(double_arr *peakList, hist_t *hist, double *summ)
{
  cutSmallPeaks(peakList, 3.0);
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.649);
  summ[1] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.809);
  summ[2] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.888);
  summ[3] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.930);
  summ[4] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.955);
  summ[5] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.970);
  return;
}

void summary_pct5(double_arr *peakList, hist_t *hist, double *summ)
{
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.969);
  summ[1] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.986);
  summ[2] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.994);
  summ[3] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.997);
  summ[4] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.999);
  return;
}

void summary_cut5(double_arr *peakList, hist_t *hist, double *summ)
{
  cutSmallPeaks(peakList, 3.0);
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.500);
  summ[1] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.7763932023);
  summ[2] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.900);
  summ[3] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.9552786405);
  summ[4] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.980);
  return;
}

void summary_pct4(double_arr *peakList, hist_t *hist, double *summ)
{
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.750);
  summ[1] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.950);
  summ[2] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.980);
  summ[3] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.998);
  return;
}

void summary_pct998(double_arr *peakList, hist_t *hist, double *summ)
{
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.998);
  return;
}

void summary_pct996(double_arr *peakList, hist_t *hist, double *summ)
{
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.996);
  return;
}

void summary_cut900(double_arr *peakList, hist_t *hist, double *summ)
{
  cutSmallPeaks(peakList, 3.0);
  gsl_sort(peakList->array, 1, peakList->length);
  summ[0] = gsl_stats_quantile_from_sorted_data(peakList->array, 1, peakList->length, 0.900);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to distance

#define invCov11 0.00985569
#define invCov22 0.01776991
#define invCov33 0.03681831
#define invCov44 0.06374419
#define invCov55 0.11458256
#define invCov66 0.05798793
double dist_abd6(double *a, double *b)
{
  double dist_sq = invCov11*pow(a[0]-b[0], 2) + invCov22*pow(a[1]-b[1], 2) + invCov33*pow(a[2]-b[2], 2)
                 + invCov44*pow(a[3]-b[3], 2) + invCov55*pow(a[4]-b[4], 2) + invCov66*pow(a[5]-b[5], 2);
  return sqrt(dist_sq);
}
#undef invCov11
#undef invCov22
#undef invCov33
#undef invCov44 
#undef invCov55
#undef invCov66

#define invCov11 0.00310011
#define invCov22 0.0112532
#define invCov33 0.02762047
#define invCov44 0.06786392
#define invCov55 0.06672277
double dist_abd5(double *a, double *b)
{
  double dist_sq = invCov11*pow(a[0]-b[0], 2) + invCov22*pow(a[1]-b[1], 2) + invCov33*pow(a[2]-b[2], 2)
                 + invCov44*pow(a[3]-b[3], 2) + invCov55*pow(a[4]-b[4], 2);
  return sqrt(dist_sq);
}
#undef invCov11
#undef invCov22
#undef invCov33
#undef invCov44 
#undef invCov55

#define invCov11 572.48978014
#define invCov22 186.75057353
#define invCov33 55.19602293
#define invCov44 18.55569778
#define invCov55 5.16860924
double dist_pct5(double *a, double *b)
{
  double dist_sq = invCov11*pow(a[0]-b[0], 2) + invCov22*pow(a[1]-b[1], 2) + invCov33*pow(a[2]-b[2], 2)
                 + invCov44*pow(a[3]-b[3], 2) + invCov55*pow(a[4]-b[4], 2);
  return sqrt(dist_sq);
}
#undef invCov11
#undef invCov22
#undef invCov33
#undef invCov44 
#undef invCov55

#define invCov11 884.77956728
#define invCov22 200.24750917
#define invCov33 57.10480582
#define invCov44 16.16071887
#define invCov55 6.59617695
double dist_cut5(double *a, double *b)
{
  double dist_sq = invCov11*pow(a[0]-b[0], 2) + invCov22*pow(a[1]-b[1], 2) + invCov33*pow(a[2]-b[2], 2)
                 + invCov44*pow(a[3]-b[3], 2) + invCov55*pow(a[4]-b[4], 2);
  return sqrt(dist_sq);
}
#undef invCov11
#undef invCov22
#undef invCov33
#undef invCov44 
#undef invCov55

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
  peak->printMode = 1; //-- Make camelus message silent
  
  //-- Initialization
  SMC_ABC_t *ABC = initialize_SMC_ABC_t(peak, err); forwardError(*err, __LINE__,);
  
  //-- Fill observation data
  char name[STRING_LENGTH_MAX];
  fillObservation_SMC_ABC_t("../demo/peakList", ABC, err); forwardError(*err, __LINE__,);
  printf("------------------------------------------------------------------------\n");
  
  //-- Loop 0
  loopZero(peak, ABC, err); forwardError(*err, __LINE__,);
  sprintf(name, "particles_%s_p%d_t%d", STR_SUMMARY_T(ABC->summ), ABC->p, ABC->t);
  FILE *file = fopen(name, "w");
  output2_particle_arr(file, ABC->newPart);
  fclose(file);
  
  //-- Other loops
  while (ABC->newPart->nbAttempts * ABC->r_stop <= ABC->p) {
    loop(peak, ABC, err); forwardError(*err, __LINE__,); //-- Switch old and new inside
    sprintf(name, "particles_%s_p%d_t%d", STR_SUMMARY_T(ABC->summ), ABC->p, ABC->t);
    FILE *file = fopen(name, "w");
    output2_particle_arr(file, ABC->newPart);
    fclose(file);
  }
  
  printf("%d iterations done\n", ABC->t+1);
  free_SMC_ABC_t(ABC);
  return;
}

//----------------------------------------------------------------------

