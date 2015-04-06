

  /*******************************
   **  peakParameters.c		**
   **  Chieh-An Lin		**
   **  Version 2015.04.05	**
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
  forwardError(*err, __LINE__,);
  return cmhm;
}
#undef NZBIN
#undef NNZ

void read_cosmo_hm(char name[], cosmo_hm **cmhm, error **err)
{
  FILE *file = fopen_err(name, "r", err);           forwardError(*err, __LINE__,);
  read_cosmological_parameters_hm(cmhm, file, err); forwardError(*err, __LINE__,);
  return;
}

void outputCosmoParam(FILE *file, cosmo_hm *cmhm, peak_param *peak)
{
  fprintf(file, "# Simulation = %s\n", peak->simulName);
  fprintf(file, "# Omega_m  Omega_de  Omega_b    n_s  h_100  sigma_8  w0_de  w1_de\n");
  fprintf(file, "#  %6.4f    %6.4f   %6.4f  %5.3f  %5.3f    %5.3f  %5.2f   %4.2f\n",
	  cmhm->cosmo->Omega_m, cmhm->cosmo->Omega_de, cmhm->cosmo->Omega_b, cmhm->cosmo->n_spec, cmhm->cosmo->h_100,
	  cmhm->cosmo->normalization, cmhm->cosmo->w0_de, cmhm->cosmo->w1_de);
  return;
}

cosmo_hm *updateCmhm(cosmo_hm *old, double Omega_m, double sigma_8, error **err)
{
   cosmo_hm *new, *buffer = old;
   new = init_parameters_hm(Omega_m, 1-Omega_m, old->cosmo->w0_de,
			    old->cosmo->w1_de,  old->cosmo->w_poly_de, old->cosmo->N_poly_de,
			    old->cosmo->h_100,  old->cosmo->Omega_b,
			    old->cosmo->Omega_nu_mass, old->cosmo->Neff_nu_mass,
			    sigma_8, old->cosmo->n_spec,
			    old->redshift->Nzbin,  old->redshift->Nnz, old->redshift->nofz,
			    old->redshift->par_nz,
			    old->zmin,             old->zmax,
			    old->cosmo->nonlinear, old->cosmo->transfer,
			    old->cosmo->growth,    old->cosmo->de_param, 
			    old->cosmo->normmode,  old->c0,  old->alpha_NFW, old->beta_NFW,
			    old->massfct,  old->halo_bias,   old->M_min,     old->M1, 
			    old->M0,       old->sigma_log_M, old->alpha, 
			    old->Mstar0,   old->beta,     old->delta,        old->gamma,        old->B_cut, old->B_sat, 
			    old->beta_cut, old->beta_sat, old->Mstellar_min, old->Mstellar_max, old->eta,
			    old->fcen1, old->fcen2,
			    old->hod,   old->pi_max, err);
   forwardError(*err, __LINE__,);
   free_parameters_hm(&buffer);
   return new;
}

//----------------------------------------------------------------------
//-- Functions related to peak_param
/*
peak_param *initialize_peak_param_default(error **err)
{
  peak_param *peak   = (peak_param*)malloc_err(sizeof(peak_param), err); forwardError(*err, __LINE__,);
  peak->field        = FIELD;
  peak->Omega[0]     = OMEGA_0; //-- [arcmin]
  peak->Omega[1]     = OMEGA_1; //-- [arcmin]
  
  peak->z_halo_max   = Z_HALO_MAX;
  peak->N_z_halo     = N_Z_HALO;
  peak->M_min        = M_MIN;   //-- [M_sol/h]
  peak->M_max        = M_MAX;   //-- [M_sol/h]
  
  peak->z_s          = Z_S;
  peak->doRandGalPos = DO_RAND_GAL_POS;
  peak->n_gal        = N_GAL;     //-- [arcmin^-2]
  peak->theta_pix    = THETA_PIX; //-- [arcmin]
  peak->sigma_eps    = SIGMA_EPS;
  
  peak->filter       = FILTER;
  peak->theta_G      = THETA_G;   //-- [arcmin]
  
  peak->simulName    = simulName;
  peak->paperII_index = -1; //-- Don't remove this.
  return peak;
  
  //-- This function is not used by main.c
  peak_param *peak = initialize_peak_param(rectangle, 409.6, 409.6,
					   1.0, 10, 1e+12, 1e+17,
					   1.0, 0, 25.0, 1.0, 0.4, 
					   gauss, 1.0,
					   "default", err);
  forwardError(*err, __LINE__,);
  return peak;
}
*/
void free_peak_param(peak_param *peak)
{
  if (peak->generator) {gsl_rng_free(peak->generator); peak->generator = NULL;}
  free(peak);
  return;
}

void read_peak_param(char name[], peak_param *peak, error **err)
{
  //-- peak should have been initialized.
  
  FILE *file = fopen(name, "r");
  testErrorRet(file==NULL, peak_null, "Peak parameter file not found", *err, __LINE__,);
  
  struct {char strField[STRING_LENGTH_MAX], strFilter[STRING_LENGTH_MAX];} string;
  config_element c = {0, 0.0, ""};
  int i;
  char buffer[STRING_LENGTH_MAX];
  
  CONFIG_READ_S(&string, strField, s, file, c, err);
  STRING2ENUM(peak->field, string.strField, field_t, STR_FIELD_T, i, NB_FIELD_T, err);
  CONFIG_READ_ARR(peak, Omega, d, i, 2, buffer, file, c, err);
  
  CONFIG_READ(peak, z_halo_max, d, file, c, err);
  CONFIG_READ(peak, N_z_halo, i, file, c, err);
  CONFIG_READ(peak, M_min, d, file, c, err);
  CONFIG_READ(peak, M_max, d, file, c, err);
  
  CONFIG_READ(peak, doKappa, i, file, c, err);
  CONFIG_READ(peak, z_s, d, file, c, err);
  CONFIG_READ(peak, doRandGalPos, i, file, c, err);
  CONFIG_READ(peak, n_gal, d, file, c, err);
  CONFIG_READ(peak, sigma_eps, d, file, c, err);
  
  CONFIG_READ(peak, theta_pix, d, file, c, err);
  CONFIG_READ_S(&string, strFilter, s, file, c, err);
  STRING2ENUM(peak->filter, string.strFilter, filter_t, STR_FILTER_T, i, NB_FILTER_T, err);
  CONFIG_READ(peak, theta_G, d, file, c, err);
  
  CONFIG_READ(peak, ABC_d, i, file, c, err);
  CONFIG_READ(peak, ABC_p, i, file, c, err);
  CONFIG_READ(peak, ABC_r_stop, d, file, c, err);
  CONFIG_READ_S(peak, ABC_summ, s, file, c, err);
  
  forwardError(*err, __LINE__,);
  peak->simulName = "default";
  fclose(file);
  return;
}

void set_peak_param(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  peak->printMode  = 0;
  peak->doNewNoise = 1;
  
  peak->generator  = initializeGenerator();
  peak->dlogM      = 0.001; //-- Bin width for mass function
  peak->nbMassBins = (int)ceil(log10(peak->M_max / peak->M_min) / peak->dlogM); //-- Number of mass bins
  
  switch (peak->field) {
    case rectangle:          peak->area = peak->Omega[0] * peak->Omega[1];       break;
    case circle:             peak->area = PI * peak->Omega[0] * peak->Omega[1];  break;
    case aardvark_hPatch04:  peak->area = peak->Omega[0] * peak->Omega[1];       break;
    case aardvark_gPatch086: peak->area = peak->Omega[0] * peak->Omega[1];       break;
    default: *err = addError(peak_unknown, "Unknown field type", *err, __LINE__); break;
  }
  forwardError(*err, __LINE__,);
  
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
  
  peak->sigma_noise = peak->sigma_eps / sqrt(FOUR_PI * peak->n_gal * SQ(peak->theta_G));
  peak->sigma_pix   = peak->sigma_eps / sqrt(2 * peak->n_gal * SQ(peak->theta_pix));
  peak->theta_G_sq  = SQ(peak->theta_G);
  
  if (peak->field == rectangle) {
    peak->s          = peak->theta_G / peak->theta_pix;
    peak->resol[0]   = (int)round(peak->Omega[0] / peak->theta_pix);
    peak->resol[1]   = (int)round(peak->Omega[1] / peak->theta_pix);
    peak->bufferSize = (int)ceil(CUTOFF_FACTOR_FILTER * peak->s) + 1;
    peak->FFTSize    = MAX(peak->resol[0], peak->resol[1]) + peak->bufferSize + 1;
  }
  else if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) {
    peak->resol[0]   = (int)round(peak->Omega[0] / peak->theta_pix);
    peak->resol[1]   = (int)round(peak->Omega[1] / peak->theta_pix);
    peak->bufferSize = (int)ceil(CUTOFF_FACTOR_FILTER * peak->s) + 1;
    //peak->cutoff      = CUTOFF_FACTOR_FILTER * peak->theta_G;
    //peak->cutoff_sq   = SQ(peak->cutoff);
    if (peak->field == aardvark_hPatch04)  peak->randPosFct = randPos_aardvark_hPatch04;
    if (peak->field == aardvark_gPatch086) peak->randPosFct = randPos_aardvark_gPatch086;
  }
    
  return;
}

void printParam(cosmo_hm *cmhm, peak_param *peak)
{
  printf("Cosmological parameters\n");
  printf("Omega_m = %.2f\n", cmhm->cosmo->Omega_m);
  printf("sigma_8 = %.2f\n", cmhm->cosmo->sigma_8);
  printf("\n");
  printf("Peak parameters\n");
  printf("Field       = %s\n", STR_FIELD_T(peak->field));
  printf("Omega       = (%.1f, %.1f) [arcmin]\n", peak->Omega[0], peak->Omega[1]);
  printf("sigma_noise = %.5f\n", peak->sigma_noise);
  printf("sigma_pix   = %.5f\n", peak->sigma_pix);
  printf("theta_G_sq  = %.3f\n", peak->theta_G_sq);
  printf("------------------------------------------------------------------------\n");
  return;
}

void printParam_ABC(cosmo_hm *cmhm, peak_param *peak)
{
  printf("ABC parameters\n");
  printf("Parameter dimension  = %d\n", peak->ABC_d);
  printf("Number of particles  = %d\n", peak->ABC_p);
  printf("Shutoff success rate = %f\n", peak->ABC_r_stop);
  printf("Summary type         = %s\n", peak->ABC_summ);
  printf("------------------------------------------------------------------------\n");
  return;
}

void printParam_complete(cosmo_hm *cmhm, peak_param *peak)
{
  printf("Cosmological parameters\n");
  printf("Omega_m = %.2f\n", cmhm->cosmo->Omega_m);
  printf("sigma_8 = %.2f\n", cmhm->cosmo->sigma_8);
  printf("\n");
  printf("Peak parameters\n");
  printf("------ ------ Field\n");
  printf("Field = %s\n", STR_FIELD_T(peak->field));
  printf("Omega = (%.1f, %.1f) [arcmin]\n", peak->Omega[0], peak->Omega[1]);
  printf("------ ------ Halos\n");
  printf("z_halo_max   = %.3f\n", peak->z_halo_max);
  printf("N_z_halo     = %d\n", peak->N_z_halo);
  printf("M_min, M_max = %.3e, %.3e\n", peak->M_min, peak->M_max);
  printf("------ ------ Galaxies\n");
  printf("z_s          = %.3f\n", peak->z_s);
  printf("doRandGalPos = %d\n", peak->doRandGalPos);
  printf("n_gal        = %.3f [arcmin^-2]\n", peak->n_gal);
  printf("sigma_eps    = %.3f\n", peak->sigma_eps);
  /*
  printf("------ ------ Galaxy survey\n");
  printf("Survey_units  = %s\n", STR_UNITS_T(peak->survey_units));
  printf("Survey_order  = (%d, %d, %d, %d, %d)\n", peak->survey_order[0], peak->survey_order[1], peak->survey_order[2], peak->survey_order[3], peak->survey_order[4]);
  printf("Survey_center = (%.3f, %.3f)\n", peak->survey_center[0], peak->survey_center[1]);
  */
  printf("------ ------ Map, filter & peaks selection\n");
  printf("theta_pix = %.3f [arcmin]\n", peak->theta_pix);
  printf("Filter    = %s\n", STR_FILTER_T(peak->filter));
  printf("theta_G   = %.3f [arcmin]\n", peak->theta_G);
  //printf("filling_factor = %.3f\n", peak->filling_factor);
  printf("------ ------ ABC\n");
  printf("Parameter dimension  = %d\n", peak->ABC_d);
  printf("Number of particles  = %d\n", peak->ABC_p);
  printf("Shutoff success rate = %f\n", peak->ABC_r_stop);
  printf("Summary type         = %s\n", peak->ABC_summ);
  printf("------ ------ Label\n");
  printf("simulName = %s\n", peak->simulName);
  printf("\n");
  printf("Precomputed part\n");
  printf("------ ------ Modes\n");
  printf("printMode   = %d\n", peak->printMode);
  printf("doKappa     = %d\n", peak->doKappa);
  printf("doNewNoise  = %d\n", peak->doNewNoise);
  printf("------ ------ Others\n");
  printf("dlogM       = %.3f\n", peak->dlogM);
  printf("nbMassBins  = %d\n", peak->nbMassBins);
  printf("area        = %.3f [arcmin^2]\n", peak->area);
  printf("w_s         = %.3f [Mpc/h]\n", peak->w_s);
  printf("D_s         = %.3f [Mpc/h]\n", peak->D_s);
  printf("sigma_noise = %.5f\n", peak->sigma_noise);
  printf("sigma_pix   = %.5f\n", peak->sigma_pix);
  printf("theta_G_sq  = %.3f [arcmin^2]\n", peak->theta_G_sq);
  printf("s           = %.3f [pix]\n", peak->s);
  printf("resol       = (%d, %d) [pix]\n", peak->resol[0], peak->resol[1]);
  printf("bufferSize  = %d [pix]\n", peak->bufferSize);
  printf("FFTSize     = %d [pix]\n", peak->FFTSize);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to position randomization

void randPos_aardvark_hPatch04(gsl_rng *generator, double pos[2])
{
  //-- Sample z and phi, z = cos(theta)
  double z_min = 1.0/3.0;
  double z_max = 11.0/12.0;
  double z_middle = 2.0/3.0;
  double z, phi, phi_min, phi_max;
  int badValues = 1;
  
  while (badValues) {
    z = gsl_ran_flat(generator, z_min, z_max);
    phi = gsl_ran_flat(generator, 0.0, 0.5);
    if (z <= z_middle) {
      phi_min = 0.5 - 0.75 * z;
      phi_max = 0.75 * z;
    }
    else {
      phi_min = 0;
      phi_max = 1 - 1.0/sqrt(12 * (1-z));
    }
    if (phi >= phi_min && phi <= phi_max) badValues = 0;
  }
  phi *= HALF_PI;
  
  //-- z, phi => x, y, z
  double sin_theta = sqrt(1 - z*z);
  double x = cos(phi) * sin_theta;
  double y = sin(phi) * sin_theta;
  
  //-- Rotation
  double x2 =  0.38197056*x + 0.05467823*y + 0.92255557*z;
  double y2 = -0.92417079*x + 0.02542496*y + 0.38113242*z;
  double z2 = -0.00261629*x - 0.99818028*y + 0.06024361*z;
  
  //-- x, y, z => RA, DEC [deg]
  double RA = atan2(y2, x2);
  double DEC = HALF_PI - acos(z2);
  pos[0] = RA  * RADIAN_TO_DEGREE;
  pos[1] = DEC * RADIAN_TO_DEGREE;
  return;
}

void randPos_aardvark_gPatch086(gsl_rng *generator, double pos[2])
{
  //-- Sample z and phi, z = cos(theta)
  double z_min = 2.0/3.0;
  double z_max = 13.0/16.0;
  double z_middle = 143.0/192.0;
  double r, z, phi, phi_min, phi_max, denom;
  int badValues = 1;
  
  while (badValues) {
    z = gsl_ran_flat(generator, z_min, z_max);
    denom = sqrt(192 * (1-z));
    r = gsl_ran_flat(generator, 0.0, 1.0);
    if (z <= z_middle) {
      phi_min = 1 - 5.0/denom;
      phi_max = 3.0/denom;
    }
    else {
      phi_min = 2.0/denom;
      phi_max = 1 - 4.0/denom;
    }
    phi = 2.0/7.0 * (1-r) + 3.0/7.0 * r;
    if (phi >= phi_min && phi <= phi_max) badValues = 0;
  }
  phi *= HALF_PI;
  
  //-- z, phi => x, y, z
  double sin_theta = sqrt(1 - z*z);
  double x = cos(phi) * sin_theta;
  double y = sin(phi) * sin_theta;
  
  //-- Rotation
  double x2 =  0.38197056*x + 0.05467823*y + 0.92255557*z;
  double y2 = -0.92417079*x + 0.02542496*y + 0.38113242*z;
  double z2 = -0.00261629*x - 0.99818028*y + 0.06024361*z;
  
  //-- x, y, z => RA, DEC [deg]
  double RA = atan2(y2, x2);
  double DEC = HALF_PI - acos(z2);
  pos[0] = RA  * RADIAN_TO_DEGREE;
  pos[1] = DEC * RADIAN_TO_DEGREE;
  return;
}

//----------------------------------------------------------------------


