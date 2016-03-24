

  /*******************************
   **  peakParameters.c		**
   **  Chieh-An Lin		**
   **  Version 2016.02.02	**
   *******************************/


#include "peakParameters.h"


//----------------------------------------------------------------------
//-- Functions related to cosmo_hm

#define NZBIN 1
#define NNZ 5
cosmo_hm *initialize_cosmo_hm_default(error **err)
{
  int Nnz[NZBIN]           = {NNZ};
  double par_nz[NZBIN*NNZ] = {0.0, 3.0, 2.0, 1.0, 0.5};
  nofz_t nofz[NZBIN]       = {ludo};
  cosmo_hm *cmhm = init_parameters_hm(0.25, 0.75, -1.0, 0.0, NULL, 0, 
				      0.78, 0.047, 0.0, 0.0, 0.80, 0.95,
				      NZBIN, Nnz, nofz, par_nz, -1, -1, 
				      smith03, eisenhu, growth_de, linder, norm_s8, 
				      11.0, 1.0, 0.13, j01, halo_bias_sc,
				      0.0, 0.0, 0.0, 0.0, 0.0,
				      1.0e11, 0.5, 0.6, 1.5, 1.5, 10.62,
				      -0.13, 0.9, 1.0e10, -1, 0.0,
				      0.0, 0.0,
				      hamana04, 60.0, err);
  forwardError(*err, __LINE__, NULL);
  return cmhm;
}

void read_cosmo_hm(char name[], cosmo_hm **cmhm, error **err)
{
  //-- This is a modified version of nicaea2.5/halomodel.c/read_cosmological_parameters_hm.
  
  FILE *file = fopen_err(name, "r", err);           forwardError(*err, __LINE__,);

  cosmo_hm *tmp;
  config_element c = {0, 0.0, ""};
  char buffer[STRING_LENGTH_MAX];
  struct {char cosmo_file[CSLENS], shod[CSLENS], smassfct[CSLENS], shalo_bias[CSLENS]; double par_nz[NNZ];} tmp2;
  int j;
  FILE *FD;

  tmp = set_cosmological_parameters_to_default_hm(err);                                          forwardError(*err, __LINE__,);

  /* Cosmology */
  CONFIG_READ_S(&tmp2, cosmo_file, s, file, c, err);                                             forwardError(*err, __LINE__,);
  FD = fopen_err(tmp2.cosmo_file, "r", err);                                                     forwardError(*err, __LINE__,);
  read_cosmological_parameters(&(tmp->cosmo), FD, err);                                          forwardError(*err, __LINE__,);
  fclose(FD);

  /* Redshift distribution */
  int Nnz[NZBIN] = {NNZ};
  nofz_t nofz[NZBIN] = {ludo};
  CONFIG_READ_ARR(&tmp2, par_nz, d, j, NNZ, buffer, file, c, err);                               forwardError(*err, __LINE__,);
  tmp->redshift = init_redshift(NZBIN, Nnz, nofz, tmp2.par_nz, NULL, err);                       forwardError(*err, __LINE__,);
  
  /* Halomodel parameters (dark matter) */
  CONFIG_READ(tmp, alpha_NFW, d, file, c, err);                                                  forwardError(*err, __LINE__,);
  CONFIG_READ(tmp, c0, d, file, c, err);                                                         forwardError(*err, __LINE__,);
  CONFIG_READ(tmp, beta_NFW, d, file, c, err);                                                   forwardError(*err, __LINE__,);
  CONFIG_READ_S(&tmp2, smassfct, s, file, c, err);                                               forwardError(*err, __LINE__,);
  STRING2ENUM(tmp->massfct, tmp2.smassfct, massfct_t, smassfct_t, j, Nmassfct_t, err);           forwardError(*err, __LINE__,);

  CONFIG_READ_S(&tmp2, shalo_bias, s, file, c, err);                                             forwardError(*err, __LINE__,);
  STRING2ENUM(tmp->halo_bias, tmp2.shalo_bias, halo_bias_t, shalo_bias_t, j, Nhalo_bias_t, err); forwardError(*err, __LINE__,);

  /* for wp(rp) */
  CONFIG_READ(tmp, pi_max, d, file, c, err);                                                     forwardError(*err, __LINE__,);
  
  /* HOD model */
  CONFIG_READ_S(&tmp2, shod, s, file, c, err);                                                   forwardError(*err, __LINE__,);
  STRING2ENUM(tmp->hod, tmp2.shod, hod_t, shod_t, j, Nhod_t, err);                               forwardError(*err, __LINE__,);
  testErrorRetVA(tmp->hod!=hod_none && tmp->hod!=hamana04 && tmp->hod!=berwein02 && tmp->hod!=berwein02_hexcl && tmp->hod!=leauthaud11,
		 hm_hodtype, "HOD type (%d) unknown", *err, __LINE__,, tmp->hod, hamana04);
  
  /* sample properties */
  CONFIG_READ(tmp, Mstellar_min, d, file, c, err); forwardError(*err, __LINE__,);
  CONFIG_READ(tmp, Mstellar_max, d, file, c, err); forwardError(*err, __LINE__,);

  /* HOD parameters */
  switch (tmp->hod) {
    case berwein02 : case berwein02_hexcl: case hamana04 :
      CONFIG_READ(tmp, M_min, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, M1, d, file, c, err);          forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, M0, d, file, c, err);          forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, sigma_log_M, d, file, c, err); forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, alpha, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, eta, d, file, c, err);         forwardError(*err, __LINE__,);
      break;

    case leauthaud11:
      CONFIG_READ(tmp, M1, d, file, c, err);          forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, Mstar0, d, file, c, err);      forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, beta, d, file, c, err);        forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, delta, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, gamma,d, file, c, err);        forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, sigma_log_M, d, file, c, err); forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, B_cut, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, B_sat, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, beta_cut, d, file, c, err);    forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, beta_sat, d, file, c, err);    forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, alpha, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, fcen1, d, file, c, err);       forwardError(*err, __LINE__,);
      CONFIG_READ(tmp, fcen2, d, file, c, err);       forwardError(*err, __LINE__,);
      break;

    default:
      break;
  }

  *cmhm = copy_parameters_hm_only(tmp, err); forwardError(*err, __LINE__,);
  return;
}
#undef NZBIN
#undef NNZ

void outputCosmoParam(FILE *file, cosmo_hm *cmhm, peak_param *peak)
{
  fprintf(file, "# Simulation = %s\n", peak->simulName);
  fprintf(file, "# Omega_m  Omega_de  Omega_b    n_s  h_100  sigma_8  w0_de  w1_de\n");
  fprintf(file, "#  %6.4f    %6.4f   %6.4f  %5.3f  %5.3f    %5.3f  %5.2f   %4.2f\n",
	  cmhm->cosmo->Omega_m, cmhm->cosmo->Omega_de, cmhm->cosmo->Omega_b, cmhm->cosmo->n_spec, cmhm->cosmo->h_100,
	  cmhm->cosmo->normalization, cmhm->cosmo->w0_de, cmhm->cosmo->w1_de);
  return;
}

cosmo_hm *updateCmhm(cosmo_hm *oldCmhm, double Omega_m, double sigma_8, double w0_de, error **err)
{
  cosmo_hm *buffer   = oldCmhm; //-- To free
  cosmo *oldCm       = oldCmhm->cosmo;
  redshift_t * oldRs = oldCmhm->redshift;
  cosmo_hm *newCmhm  = init_parameters_hm(Omega_m,           1-Omega_m,          w0_de,                 oldCm->w1_de,          oldCm->w_poly_de,         oldCm->N_poly_de,
					  oldCm->h_100,      oldCm->Omega_b,     oldCm->Omega_nu_mass,  oldCm->Neff_nu_mass,   sigma_8,                  oldCm->n_spec,
					  oldRs->Nzbin,      oldRs->Nnz,         oldRs->nofz,           oldRs->par_nz,         oldCmhm->zmin,            oldCmhm->zmax,
					  oldCm->nonlinear,  oldCm->transfer,    oldCm->growth,         oldCm->de_param,       oldCmhm->cosmo->normmode,  
					  oldCmhm->c0,       oldCmhm->alpha_NFW, oldCmhm->beta_NFW,     oldCmhm->massfct,      oldCmhm->halo_bias,   
					  oldCmhm->M_min,    oldCmhm->M1,        oldCmhm->M0,           oldCmhm->sigma_log_M,  oldCmhm->alpha, 
					  oldCmhm->Mstar0,   oldCmhm->beta,      oldCmhm->delta,        oldCmhm->gamma,        oldCmhm->B_cut,           oldCmhm->B_sat, 
					  oldCmhm->beta_cut, oldCmhm->beta_sat,  oldCmhm->Mstellar_min, oldCmhm->Mstellar_max, oldCmhm->eta,
					  oldCmhm->fcen1,    oldCmhm->fcen2,
					  oldCmhm->hod,      oldCmhm->pi_max,    err);
  forwardError(*err, __LINE__, NULL);
  free_parameters_hm(&buffer);
  return newCmhm;
}

//----------------------------------------------------------------------
//-- Functions related to peak_param

void free_peak_param(peak_param *peak)
{
  if (peak->FFT_filter)      {free(peak->FFT_filter);        peak->FFT_filter      = NULL;}
  if (peak->FFT_scale)       {free(peak->FFT_scale);         peak->FFT_scale       = NULL;}
  if (peak->DC_filter)       {free(peak->DC_filter);         peak->DC_filter       = NULL;}
  if (peak->DC_scale)        {free(peak->DC_scale);          peak->DC_scale        = NULL;}
  if (peak->nu_bin)          {free(peak->nu_bin);            peak->nu_bin          = NULL;}
  if (peak->kappa_bin)       {free(peak->kappa_bin);         peak->kappa_bin       = NULL;}
  if (peak->ABC_doParam)     {free(peak->ABC_doParam);       peak->ABC_doParam     = NULL;}
  if (peak->generator)       {gsl_rng_free(peak->generator); peak->generator       = NULL;}
  if (peak->filter)          {free(peak->filter);            peak->filter          = NULL;}
  if (peak->FFT_sigma_noise) {free(peak->FFT_sigma_noise);   peak->FFT_sigma_noise = NULL;}
  if (peak->DC_sigma_noise)  {free(peak->DC_sigma_noise);    peak->DC_sigma_noise  = NULL;}
  if (peak->DC_scale_invSq)  {free(peak->DC_scale_invSq);    peak->DC_scale_invSq  = NULL;}
  if (peak->DC_cut_sq)       {free(peak->DC_cut_sq);         peak->DC_cut_sq       = NULL;}
  if (peak->FFT_scaleInPix)  {free(peak->FFT_scaleInPix);    peak->FFT_scaleInPix  = NULL;}
  if (peak->DC_scaleInPix)   {free(peak->DC_scaleInPix);     peak->DC_scaleInPix   = NULL;}
  free(peak);
  return;
}

void read_peak_param(char name[], peak_param *peak, error **err)
{
  //-- peak should have been initialized.
  
  struct {char strField[STRING_LENGTH_MAX], FFT_strFilter[16][STRING_LENGTH_MAX], DC_strFilter[16][STRING_LENGTH_MAX];} string;
  config_element c = {0, 0.0, ""};
  char buffer[STRING_LENGTH_MAX];
  int i, j;
  
  FILE *file = fopen_err(name, "r", err);                                                    forwardError(*err, __LINE__,);
  
  //-- Halos
  CONFIG_READ(peak, z_halo_max, d, file, c, err);                                            forwardError(*err, __LINE__,);
  CONFIG_READ(peak, N_z_halo, i, file, c, err);                                              forwardError(*err, __LINE__,);
  CONFIG_READ(peak, M_min, d, file, c, err);                                                 forwardError(*err, __LINE__,);
  CONFIG_READ(peak, M_max, d, file, c, err);                                                 forwardError(*err, __LINE__,);
  
  //-- Galaxies
  CONFIG_READ(peak, z_s, d, file, c, err);                                                   forwardError(*err, __LINE__,);
  CONFIG_READ(peak, doRandGalPos, i, file, c, err);                                          forwardError(*err, __LINE__,);
  CONFIG_READ(peak, n_gal, d, file, c, err);                                                 forwardError(*err, __LINE__,);
  CONFIG_READ(peak, sigma_eps, d, file, c, err);                                             forwardError(*err, __LINE__,);
  CONFIG_READ(peak, doKappa, i, file, c, err);                                               forwardError(*err, __LINE__,);
  CONFIG_READ(peak, doMask, i, file, c, err);                                                forwardError(*err, __LINE__,);
  CONFIG_READ_S(peak, maskPath, s, file, c, err);                                            forwardError(*err, __LINE__,);
  
  //-- Field & map
  CONFIG_READ_S(&string, strField, s, file, c, err);                                         forwardError(*err, __LINE__,);
  STRING2ENUM(peak->field, string.strField, field_t, STR_FIELD_T, i, NB_FIELD_T, err);       forwardError(*err, __LINE__,);
  CONFIG_READ_ARR(peak, Omega, d, i, 2, buffer, file, c, err);                               forwardError(*err, __LINE__,);
  CONFIG_READ(peak, theta_pix, d, file, c, err);                                             forwardError(*err, __LINE__,);
  
  //-- Filter
  CONFIG_READ(peak, doSmoothing, i, file, c, err);                                           forwardError(*err, __LINE__,);
  
  CONFIG_READ(peak, FFT_nbFilters, i, file, c, err);                                         forwardError(*err, __LINE__,);
  if (peak->FFT_nbFilters != 0) {
    peak->FFT_filter = (filter_t*)malloc_err(peak->FFT_nbFilters * sizeof(filter_t), err);   forwardError(*err, __LINE__,);
    CONFIG_READ_S_ARR(&string, FFT_strFilter, buffer, i, peak->FFT_nbFilters, file, c, err); forwardError(*err, __LINE__,);
    for (j=0; j<peak->FFT_nbFilters; j++) {STRING2ENUM(peak->FFT_filter[j], string.FFT_strFilter[j], filter_t, STR_FILTER_T, i, NB_FILTER_T, err); forwardError(*err, __LINE__,);}
    peak->FFT_scale  = (double*)malloc_err(peak->FFT_nbFilters * sizeof(double), err);       forwardError(*err, __LINE__,);
    CONFIG_READ_ARR(peak, FFT_scale, d, i, peak->FFT_nbFilters, buffer, file, c, err);       forwardError(*err, __LINE__,);
  }
  else {
    peak->FFT_filter = NULL;
    read_element(file, "", c, c_s, 1, err);                                                  forwardError(*err, __LINE__,);
    peak->FFT_scale  = NULL;
    read_element(file, "", c, c_s, 1, err);                                                  forwardError(*err, __LINE__,);
  }
  
  CONFIG_READ(peak, DC_nbFilters, i, file, c, err);                                          forwardError(*err, __LINE__,);
  if (peak->DC_nbFilters != 0) {
    peak->DC_filter = (filter_t*)malloc_err(peak->DC_nbFilters * sizeof(filter_t), err);     forwardError(*err, __LINE__,);
    CONFIG_READ_S_ARR(&string, DC_strFilter, buffer, i, peak->DC_nbFilters, file, c, err);   forwardError(*err, __LINE__,);
    for (j=0; j<peak->DC_nbFilters; j++) {STRING2ENUM(peak->DC_filter[j], string.DC_strFilter[j], filter_t, STR_FILTER_T, i, NB_FILTER_T, err); forwardError(*err, __LINE__,);}
    peak->DC_scale  = (double*)malloc_err(peak->DC_nbFilters * sizeof(double), err);         forwardError(*err, __LINE__,);
    CONFIG_READ_ARR(peak, DC_scale, d, i, peak->DC_nbFilters, buffer, file, c, err);         forwardError(*err, __LINE__,);
  }
  else {
    peak->DC_filter = NULL;
    read_element(file, "", c, c_s, 1, err);                                                  forwardError(*err, __LINE__,);
    peak->DC_scale  = NULL;
    read_element(file, "", c, c_s, 1, err);                                                  forwardError(*err, __LINE__,);
  }
  
  CONFIG_READ(peak, MRLens_nbScales, i, file, c, err);                                       forwardError(*err, __LINE__,);
  CONFIG_READ(peak, MRLens_FDR, d, file, c, err);                                            forwardError(*err, __LINE__,);
  
  //-- Peak historgram
  CONFIG_READ(peak, N_nu, i, file, c, err);                                                  forwardError(*err, __LINE__,);
  peak->nu_bin = (double*)malloc_err((peak->N_nu+1) * sizeof(double), err);                  forwardError(*err, __LINE__,);
  CONFIG_READ_ARR(peak, nu_bin, d, i, peak->N_nu+1, buffer, file, c, err);                   forwardError(*err, __LINE__,);
  CONFIG_READ(peak, N_kappa, i, file, c, err);                                               forwardError(*err, __LINE__,);
  peak->kappa_bin = (double*)malloc_err((peak->N_kappa+1) * sizeof(double), err);            forwardError(*err, __LINE__,);
  CONFIG_READ_ARR(peak, kappa_bin, d, i, peak->N_kappa+1, buffer, file, c, err);             forwardError(*err, __LINE__,);
  
  //-- ABC
  CONFIG_READ_S(peak, ABC_obsPath, s, file, c, err);                                         forwardError(*err, __LINE__,);
  CONFIG_READ(peak, ABC_f, i, file, c, err);                                                 forwardError(*err, __LINE__,);
  peak->ABC_doParam = (int*)malloc_err(peak->ABC_f * sizeof(int), err);                      forwardError(*err, __LINE__,);
  CONFIG_READ_ARR(peak, ABC_doParam, i, i, peak->ABC_f, buffer, file, c, err);               forwardError(*err, __LINE__,);
  CONFIG_READ(peak, ABC_Q, i, file, c, err);                                                 forwardError(*err, __LINE__,);
  CONFIG_READ(peak, ABC_r_stop, d, file, c, err);                                            forwardError(*err, __LINE__,);
  CONFIG_READ_S(peak, ABC_summ, s, file, c, err);                                            forwardError(*err, __LINE__,);
  
  //-- Label
  sprintf(peak->simulName, "default");
  fclose(file);
  return;
}

void set_peak_param(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  //-- Default part
  peak->dlogM  = 0.001; //-- Bin width for mass function
  peak->dz_gal = 0.01;
  peak->theta_CCD_inv[0] = 60.0;
  peak->theta_CCD_inv[1] = 60.0;
  sprintf(peak->tempPath, "%s", "");
  
  //-- Precomputed part, area
  switch (peak->field) {
    case rectangle:          peak->area = peak->Omega[0] * peak->Omega[1];        break;
    case circle:             peak->area = PI * peak->Omega[0] * peak->Omega[1];   break;
    case aardvark_hPatch04:  peak->area = peak->Omega[0] * peak->Omega[1];        break;
    case aardvark_gPatch086: peak->area = peak->Omega[0] * peak->Omega[1];        break;
    default: *err = addError(peak_unknown, "Unknown field type", *err, __LINE__); break;
  }
  forwardError(*err, __LINE__,);
  
  //-- Precomputed part, distance
  int wOmegar = 0;
  if (peak->z_s > 0) {
    peak->w_s  = w(cmhm->cosmo, 1.0/(1.0 + peak->z_s), wOmegar, err); forwardError(*err, __LINE__,);
    peak->D_s  = f_K(cmhm->cosmo, peak->w_s, err);                    forwardError(*err, __LINE__,);
    peak->D_s /= 1.0 + peak->z_s;
  }
  else {
    peak->w_s = -1.0;       //-- This means that w_s needs to be calculated in set_gal_t
    peak->D_s = -1.0;       //-- This means that D_s needs to be calculated in set_gal_t
    peak->doRandGalPos = 1; //-- If z_s < 0, doRandGalPos is automatically set to 1 
  }
  
  //-- Precomputed part, others
  peak->dz_halo       = peak->z_halo_max / (double)(peak->N_z_halo);
  peak->N_M           = (int)ceil(log10(peak->M_max / peak->M_min) / peak->dlogM); //-- Number of mass bins
  peak->N_z_gal       = (int)ceil((get_zmax(cmhm->redshift, 0) - get_zmin(cmhm->redshift, 0)) / peak->dz_gal);
  peak->sigma_half    = peak->sigma_eps * SQRT_2_INV;
  peak->theta_pix_inv = 1.0 / peak->theta_pix;
  peak->sigma_pix     = peak->sigma_half / sqrt(peak->n_gal * SQ(peak->theta_pix));
  peak->generator     = initializeGenerator();
  
  //-- Filter-related
  if ((peak->doSmoothing / 1) % 2 == 0) peak->FFT_nbFilters = 0;
  if ((peak->doSmoothing / 2) % 2 == 0) peak->DC_nbFilters = 0;
  peak->doNonlinear  = (peak->doSmoothing / 4) % 2;
  peak->nbFilters    = peak->FFT_nbFilters + peak->DC_nbFilters + peak->doNonlinear;
  peak->smootherSize = MAX(peak->FFT_nbFilters, 1); //-- Keep at least one table for inversion
  peak->filter       = (filter_t*)malloc_err(peak->nbFilters * sizeof(filter_t), err); forwardError(*err, __LINE__,);
  int count = 0;
  int i;
  for (i=0; i<peak->FFT_nbFilters; i++) {
    peak->filter[count] = peak->FFT_filter[i]; 
    count++;
  }
  for (i=0; i<peak->DC_nbFilters; i++) {
    peak->filter[count] = peak->DC_filter[i];
    count++;
  }
  if (peak->doNonlinear) {
    peak->filter[count] = mrlens;
    count++;
  }
  
  //-- Filter-related, FFT
  if (peak->FFT_nbFilters) {
    peak->FFT_sigma_noise = (double*)malloc_err(peak->FFT_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
    peak->FFT_scaleInPix  = (double*)malloc_err(peak->FFT_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
    for (i=0; i<peak->FFT_nbFilters; i++) {
      //-- See DC_sigma_noise for notes on sigma_noise
      if (peak->FFT_filter[i] == gauss) {
	peak->FFT_sigma_noise[i] = peak->sigma_half / (sqrt(peak->n_gal) * peak->FFT_scale[i]) * (1.0 / sqrt(TWO_PI));
      }
      else if (peak->FFT_filter[i] == star) {
	peak->FFT_sigma_noise[i] = peak->sigma_half / (sqrt(peak->n_gal) * peak->FFT_scale[i]) * (NORM2_STARLET / NORM1_STARLET);
      }
      else {
	*err = addError(peak_unknown, "Invalid filter type for FFT", *err, __LINE__);
	forwardError(*err, __LINE__,);
      }
    }
  }
  else {
    //-- Prevent crashes in free functions
    peak->FFT_sigma_noise = NULL;
    peak->FFT_scaleInPix  = NULL;
  }
  
  if (peak->DC_nbFilters) {
    peak->DC_sigma_noise = (double*)malloc_err(peak->DC_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
    peak->DC_scale_invSq = (double*)malloc_err(peak->DC_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
    peak->DC_cut_sq      = (double*)malloc_err(peak->DC_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
    peak->DC_scaleInPix  = (double*)malloc_err(peak->DC_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
    for (i=0; i<peak->DC_nbFilters; i++) {
      //-- Let W(x) be a kernel with || W ||_L1 = NORM1 and || W ||_L2^2 = NORM2^2.
      //-- The scale-dependent kernel should be defined as W_s(i) = W(i/s) / s^2 with
      //--     || W_s ||_L1   = || W ||_L1         = NORM1
      //-- and || W_s ||_L2^2 = || W ||_L2^2 / s^2 = NORM2^2 / s^2.
      //--
      //-- Thus, sigma_noise^2 = (sigma_epsilon^2 / 2 n_gal) (|| W_s ||_L2 / || W_s ||_L1)^2
      //--                     = (sigma_half^2 / n_gal) (NORM2 / NORM1 s)^2
      //-- So, sigma_noise = (sigma_half / sqrt(n_gal) s) (NORM2 / NORM1)
      //--
      //-- Gaussian:  W(x) = exp(-x^2),    NORM1 = pi,        NORM2^2 = pi / 2, NORM2 / NORM1 = 1 / sqrt(2 pi)
      //-- Starlet:   W(x, y) = psi(x, y), NORM1 = numerical, NORM2^2 = 5 I2^2 - 2 I3^2
      //-- M_ap tanh: W(x) = tanh(x / 0.1) / [x (1 + exp(5 - 150 x) + exp(-47 + 50 x))], 
      //--                                 NORM1 = numerical, NORM2^2 = numerical
      if (peak->DC_filter[i] == gauss) {
	peak->DC_sigma_noise[i] = peak->sigma_half / (sqrt(peak->n_gal) * peak->DC_scale[i]) * (1.0 / sqrt(TWO_PI));
	peak->DC_scale_invSq[i] = 1.0 / (SQ(peak->DC_scale[i]));
	peak->DC_cut_sq[i]      = SQ(CUTOFF_FACTOR_GAUSSIAN * peak->DC_scale[i]);
      }
      else if (peak->DC_filter[i] == star) {
	peak->DC_sigma_noise[i] = peak->sigma_half / (sqrt(peak->n_gal) * peak->DC_scale[i]) * (NORM2_STARLET / NORM1_STARLET);
	peak->DC_scale_invSq[i] = 1.0 / peak->DC_scale[i]; //-- Not squared for the starlet
	peak->DC_cut_sq[i]      = 2.0 * peak->DC_scale[i]; //-- Not squared for the starlet
      }
      else if (peak->DC_filter[i] == M_ap_tanh) {
	peak->DC_sigma_noise[i] = peak->sigma_half / (sqrt(peak->n_gal) * peak->DC_scale[i]) * (NORM2_M_AP_TANH / NORM1_M_AP_TANH);
	peak->DC_scale_invSq[i] = 1.0 / (SQ(peak->DC_scale[i]));
	peak->DC_cut_sq[i]      = SQ(CUTOFF_FACTOR_M_AP_TANH * peak->DC_scale[i]);
      }
      else {
	*err = addError(peak_unknown, "Invalid filter type for DC", *err, __LINE__);
	forwardError(*err, __LINE__,);
      }
    }
  }
  else {
    //-- Prevent crashes in free functions
    peak->DC_sigma_noise = NULL;
    peak->DC_scale_invSq = NULL;
    peak->DC_cut_sq      = NULL;
    peak->DC_scaleInPix  = NULL;
  }
  
  //-- Only used if field = rectangle
  double bufferSize   = 0.0;
  if (peak->field == rectangle) {
    peak->resol[0]    = (int)round(peak->Omega[0] / peak->theta_pix);
    peak->resol[1]    = (int)round(peak->Omega[1] / peak->theta_pix);
    for (i=0; i<peak->FFT_nbFilters; i++) {
      peak->FFT_scaleInPix[i] = peak->FFT_scale[i] * peak->theta_pix_inv;
      if      (peak->FFT_filter[i] == gauss) bufferSize = MAX(bufferSize, CUTOFF_FACTOR_GAUSSIAN * peak->FFT_scaleInPix[i]);
      else if (peak->FFT_filter[i] == star)  bufferSize = MAX(bufferSize, 2.0 * peak->FFT_scaleInPix[i]);
    }
    for (i=0; i<peak->DC_nbFilters; i++) {
      peak->DC_scaleInPix[i] = peak->DC_scale[i] * peak->theta_pix_inv;
      if      (peak->DC_filter[i] == gauss)     bufferSize = MAX(bufferSize, CUTOFF_FACTOR_GAUSSIAN * peak->DC_scaleInPix[i]);
      else if (peak->DC_filter[i] == star)      bufferSize = MAX(bufferSize, 2.0 * peak->DC_scaleInPix[i]);
      else if (peak->DC_filter[i] == M_ap_tanh) bufferSize = MAX(bufferSize, CUTOFF_FACTOR_M_AP_TANH * peak->DC_scaleInPix[i]);
    }
    if (peak->doNonlinear) bufferSize = MAX(bufferSize, 8.0 * 2.0); //-- 8 pixels
    peak->bufferSize    = (int)ceil(bufferSize);                                  //-- Buffer size for filter
    peak->FFTSize       = MAX(peak->resol[0], peak->resol[1]) + peak->bufferSize; //-- FFT size with buffer area
    peak->FFTSize      += 256 - (peak->FFTSize % 256);                            //-- Round the size to a multiple of 256
    peak->FFTNormFactor = 1.0 / (double)(SQ(peak->FFTSize));
    peak->limit[0] = 0;
    peak->limit[1] = peak->Omega[0];
    peak->limit[2] = 0;
    peak->limit[3] = peak->Omega[1];
  }
  
  //-- Running part, mode
  peak->printMode   = 0; //-- 0 = detailed, 1 = no flush, 2 = line mode, 3 = MPI
  peak->doNoise     = 1; //-- 0 = noiseless, 1 = noisy
  
  //-- Running part, pipeline index
  peak->cosmoInd; //-- Should not be initialized here
  peak->realizationInd = 0;
  peak->scaleInd       = 0;
  
  //-- Running part, MPI
  peak->MPISize;  //-- Should not be initialized here
  peak->MPIInd;   //-- Should not be initialized here
  
  //-- aardvark_hPatch04 and aardvark_gPatch086 only parameters
  //-- are initialized in _paperI__initialize_peak_param_aardvark
  return;
}

//----------------------------------------------------------------------
//-- Functions related to print

void printIntArray(int *array, int length)
{
  if (length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  printf("%d", array[0]);
  for (i=1; i<length; i++) printf(", %d", array[i]);
  return;
}

void printDoubleArray(double *array, int length, int digit)
{
  if (length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  if (digit == 0) {
    printf("%f", array[0]);
    for (i=1; i<length; i++) printf(", %f", array[i]);
  }
  else if (digit == 1) {
    printf("%.1f", array[0]);
    for (i=1; i<length; i++) printf(", %.1f", array[i]);
  }
  else if (digit == 2) {
    printf("%.2f", array[0]);
    for (i=1; i<length; i++) printf(", %.2f", array[i]);
  }
  else if (digit == 3) {
    printf("%.3f", array[0]);
    for (i=1; i<length; i++) printf(", %.3f", array[i]);
  }
  else if (digit == 4) {
    printf("%.4f", array[0]);
    for (i=1; i<length; i++) printf(", %.4f", array[i]);
  }
  else if (digit == 5) {
    printf("%.5f", array[0]);
    for (i=1; i<length; i++) printf(", %.5f", array[i]);
  }
  return;
}

void printGalaxyInfo(peak_param *peak)
{
  if (peak->doRandGalPos == 0) printf("Galaxy redshift fixed at %.3f, on a regular grid (%dx%d)\n", peak->z_s, peak->resol[0], peak->resol[1]);
  else if (peak->z_s > 0)      printf("Galaxy redshift fixed at %.3f, random angular position\n", peak->z_s);
  else                         printf("Random galaxy redshift, random angular position\n");
  return;
}

void printFilterArray(filter_t *array, int length)
{
  if (length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  printf("%s", STR_FILTER_T(array[0]));
  for (i=1; i<length; i++) printf(", %s", STR_FILTER_T(array[i]));
  return;
}

void printParam(cosmo_hm *cmhm, peak_param *peak)
{
  printf("Cosmological parameters\n");
  printf("Omega_m = %g\n", cmhm->cosmo->Omega_m);
  printf("sigma_8 = %g\n", cmhm->cosmo->sigma_8);
  printf("w0_de   = %g\n", cmhm->cosmo->w0_de);
  printf("\n");
  printf("Peak parameters\n");
  printGalaxyInfo(peak);
  printf("Omega           = (%.1f, %.1f) [arcmin]\n", peak->Omega[0], peak->Omega[1]);
  printf("theta_pix       = %.3f [arcmin]\n", peak->theta_pix);
  printf("Resolution      = (%d, %d) [pix]\n", peak->resol[0], peak->resol[1]);
  printf("Buffer size     = %d [pix]\n", peak->bufferSize);
  printf("Peak field      = (%d, %d) [pix]\n", peak->resol[0] - 2*peak->bufferSize, peak->resol[1] - 2*peak->bufferSize);
  printf("doKappa         = %d\n", peak->doKappa);
  if (peak->FFT_nbFilters) {
    printf("FFT_filter      = "); printFilterArray(peak->FFT_filter, peak->FFT_nbFilters); printf("\n");
    printf("FFT_scale       = "); printDoubleArray(peak->FFT_scale, peak->FFT_nbFilters, 3); printf(" [arcmin]\n");
  }
  if (peak->DC_nbFilters) {
    printf("DC_filter       = "); printFilterArray(peak->DC_filter, peak->DC_nbFilters); printf("\n");
    printf("DC_scale        = "); printDoubleArray(peak->DC_scale, peak->DC_nbFilters, 3); printf(" [arcmin]\n");
  }
  if (peak->doNonlinear) {
    printf("MRLens_nbScales = %d\n", peak->MRLens_nbScales);
    printf("MRLens_FDR      = %f\n", peak->MRLens_FDR);
  }
  printf("------------------------------------------------------------------------\n");
  return;
}

void printParam_complete(cosmo_hm *cmhm, peak_param *peak)
{
  int i;
  printf("Cosmological parameters\n");
  printf("Omega_m         = %.2f [-]\n", cmhm->cosmo->Omega_m);
  printf("sigma_8         = %.2f [-]\n", cmhm->cosmo->sigma_8);
  printf("w0_de           = %.2f [-]\n", cmhm->cosmo->w0_de);
  printf("par_nz          = "); printDoubleArray(cmhm->redshift->par_nz, cmhm->redshift->Nnz_max, 1); printf("\n");
  printGalaxyInfo(peak);
  printf("\n");
  printf("Peak parameters, customized part\n");
  printf("------ Halos ------\n");
  printf("z_halo_max      = %.3f [-]\n", peak->z_halo_max);
  printf("N_z_halo        = %d\n", peak->N_z_halo);
  printf("M_min, M_max    = %.3e, %.3e [M_sol/h]\n", peak->M_min, peak->M_max);
  printf("------ Galaxies ------\n");
  printf("z_s             = %.3f [-]\n", peak->z_s);
  printf("doRandGalPos    = %d\n", peak->doRandGalPos);
  printf("n_gal           = %.3f [arcmin^-2]\n", peak->n_gal);
  printf("sigma_eps       = %.3f [-]\n", peak->sigma_eps);
  printf("doKappa         = %d\n", peak->doKappa);
  printf("doMask          = %d\n", peak->doMask);
  printf("maskPath        = \"%s\"\n", peak->maskPath);
  printf("------ Field & map ------\n");
  printf("field           = %s\n", STR_FIELD_T(peak->field));
  printf("Omega           = (%.1f, %.1f) [arcmin]\n", peak->Omega[0], peak->Omega[1]);
  printf("theta_pix       = %.3f [arcmin]\n", peak->theta_pix);
  printf("------ Filter ------\n");
  printf("doSmoothing     = %d\n", peak->doSmoothing);
  printf("FFT_nbFilters   = %d\n", peak->FFT_nbFilters);
  printf("FFT_filter      = "); printFilterArray(peak->FFT_filter, peak->FFT_nbFilters); printf("\n");
  printf("FFT_scale       = "); printDoubleArray(peak->FFT_scale, peak->FFT_nbFilters, 3); printf(" [arcmin]\n");
  printf("DC_nbFilters    = %d\n", peak->DC_nbFilters);
  printf("DC_filter       = "); printFilterArray(peak->DC_filter, peak->DC_nbFilters); printf("\n");
  printf("DC_scale        = "); printDoubleArray(peak->DC_scale, peak->DC_nbFilters, 3); printf(" [arcmin]\n");
  printf("MRLens_nbScales = %d\n", peak->MRLens_nbScales);
  printf("MRLens_FDR      = %f\n", peak->MRLens_FDR);
  printf("------ Peak historgram ------\n");
  printf("N_nu            = %d\n", peak->N_nu);
  printf("nu_bin          = "); printDoubleArray(peak->nu_bin, peak->N_nu+1, 3); printf("\n");
  printf("N_kappa         = %d\n", peak->N_kappa);
  printf("kappa_bin       = "); printDoubleArray(peak->kappa_bin, peak->N_kappa+1, 3); printf("\n");
  printf("------ ABC ------\n");
  printf("ABC_obsPath     = \"%s\"\n", peak->ABC_obsPath);
  printf("ABC_f           = %d\n", peak->ABC_f);
  printf("ABC_doParam     = "); printIntArray(peak->ABC_doParam, peak->ABC_f); printf("\n");
  printf("ABC_Q           = %d\n", peak->ABC_Q);
  printf("ABC_r_stop      = %f\n", peak->ABC_r_stop);
  printf("ABC_summ        = %s\n", peak->ABC_summ);
  printf("------ Label ------\n");
  printf("simulName       = %s\n", peak->simulName);
  printf("\n");
  printf("Peak parameters, default part\n");
  printf("dlogM           = %.3f [-]\n", peak->dlogM);
  printf("dz_gal          = %.3f [-]\n", peak->dz_gal);
  printf("theta_CCD_inv   = %f, %f [arcmin^-1]\n", peak->theta_CCD_inv[0], peak->theta_CCD_inv[1]);
  printf("tempPath        = \"%s\"\n", peak->tempPath);
  printf("\n");
  printf("Peak parameters, precomputed part\n");
  printf("area            = %.3f [arcmin^2]\n", peak->area);
  printf("w_s             = %.3f [Mpc/h]\n", peak->w_s);
  printf("D_s             = %.3f [Mpc/h]\n", peak->D_s);
  printf("dz_halo         = %.3f [-]\n", peak->dz_halo);
  printf("N_M             = %d\n", peak->N_M);
  printf("N_z_gal         = %d\n", peak->N_z_gal);
  printf("sigma_half      = %.3f [-]\n", peak->sigma_half);
  printf("theta_pix_inv   = %.3f [arcmin^-1]\n", peak->theta_pix_inv);
  printf("sigma_pix       = %f [-]\n", peak->sigma_pix);
  printf("------ Filter-related ------\n");
  printf("doNonlinear     = %d\n", peak->doNonlinear);
  printf("nbFilters       = %d\n", peak->nbFilters);
  printf("smootherSize    = %d\n", peak->smootherSize);
  printf("filter          = "); printFilterArray(peak->filter, peak->nbFilters); printf("\n");
  printf("FFT_sigma_noise = "); printDoubleArray(peak->FFT_sigma_noise, peak->FFT_nbFilters, 0); printf(" [-]\n");
  printf("DC_sigma_noise  = "); printDoubleArray(peak->DC_sigma_noise, peak->DC_nbFilters, 0); printf(" [-]\n");
  printf("DC_scale_invSq  = "); printDoubleArray(peak->DC_scale_invSq, peak->DC_nbFilters, 3); printf(" [arcmin^-2]\n");
  printf("DC_cut_sq       = "); printDoubleArray(peak->DC_cut_sq, peak->DC_nbFilters, 3); printf(" [arcmin^2]\n");
  printf("------ Only used if field = rectangle ------\n");
  printf("resol           = (%d, %d) [pix]\n", peak->resol[0], peak->resol[1]);
  printf("FFT_scaleInPix  = "); printDoubleArray(peak->FFT_scaleInPix, peak->FFT_nbFilters, 3); printf(" [pix]\n");
  printf("DC_scaleInPix   = "); printDoubleArray(peak->DC_scaleInPix, peak->DC_nbFilters, 3); printf(" [pix]\n");
  printf("bufferSize      = %d [pix]\n", peak->bufferSize);
  printf("FFTSize         = %d [pix]\n", peak->FFTSize);
  printf("FFTNormFactor   = %.9f [-]\n", peak->FFTNormFactor);
  printf("limit           = (%.1f, %.1f, %.1f, %.1f) [arcmin]\n", peak->limit[0], peak->limit[1], peak->limit[2], peak->limit[3]);
  printf("\n");
  printf("Peak parameters, running part\n");
  printf("------ Mode ------\n");
  printf("printMode       = %d\n", peak->printMode);
  printf("doNoise         = %d\n", peak->doNoise);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------

