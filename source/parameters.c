

  /*******************************************************
   **  parameters.c					**
   **  Version 2018.03.14				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "parameters.h"


//----------------------------------------------------------------------
//-- Functions related to parameter reading

void ignoreComments(char line[])
{
  char *end = strchr(line, '#');
  if (end == line) sprintf(end, "%s", "");
  else if (end != NULL) sprintf(end, "\n");
  return;
}

int getKeyAndValues(char line[], char kv[][256])
{
  if (!strcmp(line, "-h") || !strcmp(line, "-H") || !strcmp(line, "--help")) return -1; //-- This is for parameters updated from the command line.
  
  char buffer[STRING_LENGTH_MAX];
  sprintf(buffer, "%s", line);
  
  int count = 0;
  char *token = strtok(buffer, " ,=\"\t\n");
  while (token != NULL) {
    sprintf(kv[count], "%s", token);
    count++;
    token = strtok(NULL, " ,=\"\t\n");
  }
  return count;
}

int *makeIntArray(char kv[][256], int count, error **err)
{
  int *array = (count == 1) ? NULL : (int*)malloc_err((count-1) * sizeof(int), err); forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<count-1; i++) array[i] = atoi(kv[i+1]);
  return array;
}

double *makeDoubleArray(char kv[][256], int count, error **err)
{
  double *array = (count == 1) ? NULL : (double*)malloc_err((count-1) * sizeof(double), err); forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<count-1; i++) array[i] = atof(kv[i+1]);
  return array;
}

void setPathWhichCanBeBlank(char path[], char kv[][256], int count)
{
  if (count == 1) sprintf(path, "%s", "");
  else            sprintf(path, "%s", kv[1]);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to cosmo_hm

#define NZBIN 1
#define NNZ 5
cosmo_hm *initialize_cosmo_hm_default(error **err)
{
  int Nnz[NZBIN]           = {NNZ};
  double par_nz[NZBIN*NNZ] = {0.0, 3.0, 2.0, 1.0, 0.5};
  nofz_t nofz[NZBIN]       = {ludo};
  photz_t photz[NZBIN]     = {photz_no};
  cosmo_hm *chPar = init_parameters_hm(0.25, 0.75, -1.0, 0.0, NULL, 0, 
				       0.78, 0.047, 0.0, 0.0, 0.80, 0.95,
				       NZBIN, Nnz, nofz, photz, par_nz, -1, -1, 
				       smith03_revised, eisenhu_osc, growth_de, linder, norm_s8, 
				       11.0, 1.0, 0.13, j01, halo_bias_sc,
				       11.0, 12.0, 11.0, 0.3, 0.75,
				       11.0, 0.5, 0.6, 1.5, 1.5, 10.62,
				       -0.13, 0.9, 10.0, -1, 1.0,
				       0.15, 0.5, hamana04, 60.0, err);
  forwardError(*err, __LINE__, NULL);
  return chPar;
}
#undef NZBIN
#undef NNZ

int findKey_cosmo(cosmo *cmPar, peak_param *pkPar, char kv[][256], int count, error **err)
{
  int count2 = 0;
  char buffer[STRING_LENGTH_MAX];
  int j;
  
  if (count == 0) return 0;
  
  //-- Cosmology
  if (!strcmp(kv[0], "Omega_m"))       {pkPar->cmParTable[count2] = 1; cmPar->Omega_m       = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "Omega_de"))      {pkPar->cmParTable[count2] = 1; cmPar->Omega_de      = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "w0_de"))         {pkPar->cmParTable[count2] = 1; cmPar->w0_de         = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "w1_de"))         {pkPar->cmParTable[count2] = 1; cmPar->w1_de         = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "h_100"))         {pkPar->cmParTable[count2] = 1; cmPar->h_100         = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "Omega_b"))       {pkPar->cmParTable[count2] = 1; cmPar->Omega_b       = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "Omega_nu_mass")) {pkPar->cmParTable[count2] = 1; cmPar->Omega_nu_mass = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "Neff_nu_mass"))  {pkPar->cmParTable[count2] = 1; cmPar->Neff_nu_mass  = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "normalization") || !strcmp(kv[0], "sigma_8")) {pkPar->cmParTable[count2] = 1; cmPar->normalization = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "n_spec"))        {pkPar->cmParTable[count2] = 1; cmPar->n_spec        = atof(kv[1]); return 0;} else count2++;
  
  if (!strcmp(kv[0], "snonlinear"))    {
                                        pkPar->cmParTable[count2] = 1;
                                        sprintf(buffer, "%s", kv[1]);
                                        STRING_TO_ENUM(cmPar->nonlinear, buffer, nonlinear_t, snonlinear_t, j, Nnonlinear_t, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "stransfer"))     {
                                        pkPar->cmParTable[count2] = 1;
                                        sprintf(buffer, "%s", kv[1]);
                                        STRING_TO_ENUM(cmPar->transfer, buffer, transfer_t, stransfer_t, j, Ntransfer_t, err);     forwardError(*err, __LINE__, -1);
                                        return 0;
  } else count2++;
  if (!strcmp(kv[0], "sgrowth"))       {
                                        pkPar->cmParTable[count2] = 1;
                                        sprintf(buffer, "%s", kv[1]);
                                        STRING_TO_ENUM(cmPar->growth, buffer, growth_t, sgrowth_t, j, Ngrowth_t, err);             forwardError(*err, __LINE__, -1);
                                        return 0;
  } else count2++;
  if (!strcmp(kv[0], "sde_param"))     {
                                        pkPar->cmParTable[count2] = 1;
                                        sprintf(buffer, "%s", kv[1]);
                                        STRING_TO_ENUM(cmPar->de_param, buffer, de_param_t, sde_param_t, j, Nde_param_t, err);     forwardError(*err, __LINE__, -1);
                                        return 0;
  } else count2++;
  if (!strcmp(kv[0], "normmode"))      {pkPar->cmParTable[count2] = 1; cmPar->normmode      = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "a_min"))         {pkPar->cmParTable[count2] = 1; cmPar->a_min         = atof(kv[1]); return 0;} else count2++;
  
  return 1;
}

void read_cosmo(char name[], cosmo *cmPar, peak_param *pkPar, error **err)
{
  //-- cmPar should have been initialized.
  
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char line[STRING_LENGTH_MAX], kv[STRING_LENGTH_MAX][256], *buffer;
  int count;
  
  //-- Read
  do {
    buffer = fgets(line, STRING_LENGTH_MAX, file);
    ignoreComments(line);
    count = getKeyAndValues(line, kv);
    findKey_cosmo(cmPar, pkPar, kv, count, err); forwardError(*err, __LINE__,);
  }
  while (buffer != NULL);
  
  fclose(file);
  return;
}

#define NZBIN 1
#define NNZ 5
int findKey_cosmo_hm(cosmo_hm *chPar, peak_param *pkPar, char kv[][256], int count, error **err)
{
  int count2 = 0;
  char buffer[STRING_LENGTH_MAX];
  double buffer2[NNZ];
  int j, unknown;
  
  if (count == 0) return 0;
  
  //-- Cosmology
  unknown = findKey_cosmo(chPar->cosmo, pkPar, kv, count, err); forwardError(*err, __LINE__, -1);
  
  //-- Parameter file
  if (!strcmp(kv[0], "cosmo_file")) {
                                     pkPar->hmParTable[count2] = 1;
                                     sprintf(pkPar->cmParPath, "%s", kv[1]);
                                     read_cosmo(pkPar->cmParPath, chPar->cosmo, pkPar, err); forwardError(*err, __LINE__, -1);
                                     return 0;
  } else count2++;
  
  //-- Galaxy redshift distribution
  if (!strcmp(kv[0], "par_nz"))     {
                                     pkPar->hmParTable[count2] = 1;
                                     if (chPar->redshift) free_redshift(&chPar->redshift);
                                     for (j=0; j<NNZ; j++) buffer2[j] = atof(kv[j+1]);
                                     int Nnz[NZBIN]       = {NNZ};
                                     nofz_t nofz[NZBIN]   = {ludo};
                                     photz_t photz[NZBIN] = {photz_no};
                                     chPar->redshift      = init_redshift(NZBIN, Nnz, nofz, photz, buffer2, NULL, err);         forwardError(*err, __LINE__, -1);
                                     return 0;
  } else count2++;
  
  //-- Halo model
  if (!strcmp(kv[0], "alpha_NFW"))  {pkPar->hmParTable[count2] = 1; chPar->alpha_NFW = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "c0"))         {pkPar->hmParTable[count2] = 1; chPar->c0        = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "beta_NFW"))   {pkPar->hmParTable[count2] = 1; chPar->beta_NFW  = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "smassfct"))   {
                                     pkPar->hmParTable[count2] = 1;
                                     sprintf(buffer, "%s", kv[1]);
                                     STRING_TO_ENUM(chPar->massfct, buffer, massfct_t, smassfct_t, j, Nmassfct_t, err);         forwardError(*err, __LINE__, -1);
                                     return 0;
  } else count2++;
  if (!strcmp(kv[0], "shalo_bias")) {
                                     pkPar->hmParTable[count2] = 1;
                                     sprintf(buffer, "%s", kv[1]);
                                     STRING_TO_ENUM(chPar->halo_bias, buffer, halo_bias_t, shalo_bias_t, j, Nhalo_bias_t, err); forwardError(*err, __LINE__, -1);
                                     return 0;
  } else count2++;
  
  //-- HOD
  if (!strcmp(kv[0], "shod")) {
                                     sprintf(buffer, "%s", kv[1]);
                                     STRING_TO_ENUM(chPar->hod, buffer, hod_t, shod_t, j, Nhod_t, err);                         forwardError(*err, __LINE__, -1);
  }
  else if (!strcmp(kv[0], "log10M_min"))     chPar->log10M_min     = atof(kv[1]);
  else if (!strcmp(kv[0], "log10M1"))        chPar->log10M1        = atof(kv[1]);
  else if (!strcmp(kv[0], "log10M0"))        chPar->log10M0        = atof(kv[1]);
  else if (!strcmp(kv[0], "sigma_log_M"))    chPar->sigma_log_M    = atof(kv[1]);
  else if (!strcmp(kv[0], "alpha"))          chPar->alpha          = atof(kv[1]);
  else if (!strcmp(kv[0], "eta"))            chPar->eta            = atof(kv[1]);
  else if (!strcmp(kv[0], "pi_max"))         chPar->pi_max         = atof(kv[1]);
  else if (!strcmp(kv[0], "log10Mstar_min")) chPar->log10Mstar_min = atof(kv[1]);
  else if (!strcmp(kv[0], "log10Mstar_max")) chPar->log10Mstar_max = atof(kv[1]);
  else return (int)(unknown > 0);
  
  return 0;
}
#undef NZBIN
#undef NNZ

void read_cosmo_hm(char name[], cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  //-- chPar should have been initialized.
  
  FILE *file = fopen_err(name, "r", err);    forwardError(*err, __LINE__,);
  char line[STRING_LENGTH_MAX], kv[STRING_LENGTH_MAX][256], *buffer;
  int count, unknown;
  
  //-- Read
  do {
    buffer = fgets(line, STRING_LENGTH_MAX, file);
    ignoreComments(line);
    count = getKeyAndValues(line, kv);
    unknown = findKey_cosmo_hm(chPar, pkPar, kv, count, err); forwardError(*err, __LINE__,);
    if (unknown == 1) printf("Detected the keyword \"%s\" but unable to interpret, skipped\n", kv[0]);
  }
  while (buffer != NULL);
  
  fclose(file);
  return;
}

cosmo_hm *reinitialize_cosmo_hm(cosmo_hm *chPar, error **err)
{
  cosmo_hm *newChPar = copy_parameters_hm_only(chPar, err); forwardError(*err, __LINE__, NULL);
  free_parameters_hm(&chPar);
  return newChPar;
}

void outAsciiCosmoParam(FILE *file, cosmo_hm *chPar, peak_param *pkPar)
{
  fprintf(file, "# Cosmology = %s\n", pkPar->cmParPath);
  fprintf(file, "# Omega_m  Omega_de  Omega_b    n_s  h_100  sigma_8  w0_de  w1_de\n");
  fprintf(file, "#  %6.4f    %6.4f   %6.4f  %5.3f  %5.3f    %5.3f  %5.2f   %4.2f\n",
	  chPar->cosmo->Omega_m, chPar->cosmo->Omega_de, chPar->cosmo->Omega_b, chPar->cosmo->n_spec, chPar->cosmo->h_100,
	  chPar->cosmo->normalization, chPar->cosmo->w0_de, chPar->cosmo->w1_de);
  return;
}

#ifdef __CAMELUS_USE_FITS__
void outFitsCosmoParam(FITS_t *fits, cosmo_hm *chPar, peak_param *pkPar)
{
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "COSMO",   pkPar->cmParPath,             "[-] Path of the cosmological parameter file");
  addKeyword(fits, TDOUBLE, "OMEGAM",  &chPar->cosmo->Omega_m,       "[-]");
  addKeyword(fits, TDOUBLE, "OMEGADE", &chPar->cosmo->Omega_de,      "[-]");
  addKeyword(fits, TDOUBLE, "OMEGAB",  &chPar->cosmo->Omega_b,       "[-]");
  addKeyword(fits, TDOUBLE, "NS",      &chPar->cosmo->n_spec,        "[-]");
  addKeyword(fits, TDOUBLE, "H100",    &chPar->cosmo->h_100,         "[-]");
  addKeyword(fits, TDOUBLE, "SIGMA8",  &chPar->cosmo->normalization, "[-]");
  addKeyword(fits, TDOUBLE, "W0DE",    &chPar->cosmo->w0_de,         "[-]");
  addKeyword(fits, TDOUBLE, "W1DE",    &chPar->cosmo->w1_de,         "[-]");
  return;
}
#endif

//----------------------------------------------------------------------
//-- Functions related to peak_param

peak_param *initialize_peak_param(error **err)
{
  peak_param *pkPar = (peak_param*)malloc_err(sizeof(peak_param), err); forwardError(*err, __LINE__, NULL);
  return pkPar;
}

void free_peak_param(peak_param *pkPar)
{
  if (pkPar) {
    if (pkPar->FFT_filter)      {free(pkPar->FFT_filter);        pkPar->FFT_filter      = NULL;}
    if (pkPar->FFT_scale)       {free(pkPar->FFT_scale);         pkPar->FFT_scale       = NULL;}
    if (pkPar->DC_filter)       {free(pkPar->DC_filter);         pkPar->DC_filter       = NULL;}
    if (pkPar->DC_scale)        {free(pkPar->DC_scale);          pkPar->DC_scale        = NULL;}
    if (pkPar->bin_nu)          {free(pkPar->bin_nu);            pkPar->bin_nu          = NULL;}
    if (pkPar->ABC_doParam)     {free(pkPar->ABC_doParam);       pkPar->ABC_doParam     = NULL;}
    if (pkPar->generator)       {gsl_rng_free(pkPar->generator); pkPar->generator       = NULL;}
    if (pkPar->filter)          {free(pkPar->filter);            pkPar->filter          = NULL;}
    if (pkPar->FFT_scaleInPix)  {free(pkPar->FFT_scaleInPix);    pkPar->FFT_scaleInPix  = NULL;}
    if (pkPar->FFT_sigma_noise) {free(pkPar->FFT_sigma_noise);   pkPar->FFT_sigma_noise = NULL;}
    if (pkPar->DC_scaleInPix)   {free(pkPar->DC_scaleInPix);     pkPar->DC_scaleInPix   = NULL;}
    if (pkPar->DC_scale_inv)    {free(pkPar->DC_scale_inv);      pkPar->DC_scale_inv    = NULL;}
    if (pkPar->DC_cut)          {free(pkPar->DC_cut);            pkPar->DC_cut          = NULL;}
    if (pkPar->DC_sigma_noise)  {free(pkPar->DC_sigma_noise);    pkPar->DC_sigma_noise  = NULL;}
    free(pkPar); pkPar = NULL;
  }
  return;
}

int findKey_peak_param(peak_param *pkPar, char kv[][256], int count, error **err)
{
  int count2 = 0;
  int i;
  
  if (count == 0) return 0;
  
  //-- Parameter files
  if (!strcmp(kv[0], "hmParPath"))      {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->hmParPath, kv, count); return 0;} else count2++;
  if (!strcmp(kv[0], "seed"))           {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->seed, kv, count);      return 0;} else count2++;
  if (!strcmp(kv[0], "verbose"))        {pkPar->pkParTable[count2] = 1; pkPar->verbose = atoi(kv[1]);                        return 0;} else count2++;
  
  //-- Field
  if (!strcmp(kv[0], "field"))          {pkPar->pkParTable[count2] = 1; pkPar->field = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "Omega"))          {
                                         pkPar->pkParTable[count2] = 1;
                                         pkPar->Omega[0]  = atof(kv[1]) * ARCMIN_TO_RADIAN;
                                         pkPar->Omega[1]  = atof(kv[2]) * ARCMIN_TO_RADIAN;
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "theta_pix"))      {pkPar->pkParTable[count2] = 1; pkPar->theta_pix = atof(kv[1]) * ARCMIN_TO_RADIAN; return 0;} else count2++;
  if (!strcmp(kv[0], "nside"))          {pkPar->pkParTable[count2] = 1; pkPar->nside     = strtol(kv[1], NULL, 10);        return 0;} else count2++;
  if (!strcmp(kv[0], "patch"))          {pkPar->pkParTable[count2] = 1; pkPar->patch     = strtol(kv[1], NULL, 10);        return 0;} else count2++;
  if (!strcmp(kv[0], "rotAng"))         {pkPar->pkParTable[count2] = 1; pkPar->rotAng    = atof(kv[1]) * DEGREE_TO_RADIAN; return 0;} else count2++;
  if (!strcmp(kv[0], "HP_resol"))       {pkPar->pkParTable[count2] = 1; pkPar->HP_resol  = atof(kv[1]);                    return 0;} else count2++;
  
  //-- Halos
  if (!strcmp(kv[0], "inHaloCatPath"))  {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->inHaloCatPath, kv, count); return 0;} else count2++;
  if (!strcmp(kv[0], "inHaloCatCol"))   {
                                         pkPar->pkParTable[count2] = 1; 
                                         pkPar->inHaloCatCol[0] = atoi(kv[1]);
                                         pkPar->inHaloCatCol[1] = atoi(kv[2]);
                                         pkPar->inHaloCatCol[2] = atoi(kv[3]);
                                         pkPar->inHaloCatCol[3] = atoi(kv[4]);
                                         pkPar->inHaloCatCol[4] = atoi(kv[5]);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "z_halo_min"))     {pkPar->pkParTable[count2] = 1; pkPar->z_halo_min = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "z_halo_max"))     {pkPar->pkParTable[count2] = 1; pkPar->z_halo_max = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "N_z_halo"))       {pkPar->pkParTable[count2] = 1; pkPar->N_z_halo   = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "M_min"))          {pkPar->pkParTable[count2] = 1; pkPar->M_min      = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "M_max"))          {pkPar->pkParTable[count2] = 1; pkPar->M_max      = atof(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "dlogM"))          {pkPar->pkParTable[count2] = 1; pkPar->dlogM      = atof(kv[1]); return 0;} else count2++;
  
  //-- Galaxies
  if (!strcmp(kv[0], "inGalCatPath"))   {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->inGalCatPath, kv, count); return 0;} else count2++;
  if (!strcmp(kv[0], "inGalCatCol"))    {
                                         pkPar->pkParTable[count2] = 1; 
                                         pkPar->inGalCatCol[0] = atoi(kv[1]);
                                         pkPar->inGalCatCol[1] = atoi(kv[2]);
                                         pkPar->inGalCatCol[2] = atoi(kv[3]);
                                         pkPar->inGalCatCol[3] = atoi(kv[4]);
                                         pkPar->inGalCatCol[4] = atoi(kv[5]);
                                         pkPar->inGalCatCol[5] = atoi(kv[6]);
                                         pkPar->inGalCatCol[6] = atoi(kv[7]);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "z_s"))            {pkPar->pkParTable[count2] = 1; pkPar->z_s          = atof(kv[1]);                          return 0;} else count2++;
  if (!strcmp(kv[0], "dz_gal"))         {pkPar->pkParTable[count2] = 1; pkPar->dz_gal       = atof(kv[1]);                          return 0;} else count2++;
  if (!strcmp(kv[0], "doRandGalPos"))   {pkPar->pkParTable[count2] = 1; pkPar->doRandGalPos = atoi(kv[1]);                          return 0;} else count2++;
  if (!strcmp(kv[0], "n_gal"))          {pkPar->pkParTable[count2] = 1; pkPar->n_gal        = atof(kv[1]) / ARCMIN_SQ_TO_RADIAN_SQ; return 0;} else count2++;
  if (!strcmp(kv[0], "sigma_eps"))      {pkPar->pkParTable[count2] = 1; pkPar->sigma_eps    = atof(kv[1]);                          return 0;} else count2++;
  
  //-- Masks
  if (!strcmp(kv[0], "doMask"))         {pkPar->pkParTable[count2] = 1; pkPar->doMask      = atoi(kv[1]);                   return 0;} else count2++;
  if (!strcmp(kv[0], "maskPath"))       {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->maskPath, kv, count); return 0;} else count2++;
  if (!strcmp(kv[0], "nbHoleTypes"))    {pkPar->pkParTable[count2] = 1; pkPar->nbHoleTypes = atoi(kv[1]);                   return 0;} else count2++;
  if (!strcmp(kv[0], "holeRadius"))     {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->holeRadius) free(pkPar->holeRadius); 
                                         pkPar->holeRadius = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         for (i=0; i<count-1; i++) pkPar->holeRadius[i] *= ARCMIN_TO_RADIAN;
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "holeDensity"))    {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->holeDensity) free(pkPar->holeDensity);
                                         pkPar->holeDensity = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         for (i=0; i<count-1; i++) pkPar->holeDensity[i] /= DEGREE_SQ_TO_RADIAN_SQ;
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "stripeLRatio"))   {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->stripeLRatio) free(pkPar->stripeLRatio);
                                         pkPar->stripeLRatio = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "stripeWRatio"))   {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->stripeWRatio) free(pkPar->stripeWRatio);
                                         pkPar->stripeWRatio = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  
  //-- Lensing
  if (!strcmp(kv[0], "doLensing"))      {pkPar->pkParTable[count2] = 1; pkPar->doLensing     = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "doKappa"))        {pkPar->pkParTable[count2] = 1; pkPar->doKappa       = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "doSubtraction"))  {pkPar->pkParTable[count2] = 1; pkPar->doSubtraction = atoi(kv[1]); return 0;} else count2++;
  
  //-- Filters
  if (!strcmp(kv[0], "doSmoothing"))    {pkPar->pkParTable[count2] = 1; pkPar->doSmoothing   = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "FFT_nbFilters"))  {pkPar->pkParTable[count2] = 1; pkPar->FFT_nbFilters = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "FFT_filter"))     {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->FFT_filter) free(pkPar->FFT_filter); 
                                         pkPar->FFT_filter = makeIntArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "FFT_scale"))      {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->FFT_scale) free(pkPar->FFT_scale);
                                         pkPar->FFT_scale = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         for (i=0; i<count-1; i++) pkPar->FFT_scale[i] *= ARCMIN_TO_RADIAN;
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "DC_nbFilters"))   {pkPar->pkParTable[count2] = 1; pkPar->DC_nbFilters = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "DC_filter"))      {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->DC_filter) free(pkPar->DC_filter); 
                                         pkPar->DC_filter = makeIntArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "DC_scale"))       {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->DC_scale) free(pkPar->DC_scale); 
                                         pkPar->DC_scale = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         for (i=0; i<count-1; i++) pkPar->DC_scale[i] *= ARCMIN_TO_RADIAN;
                                         return 0;
  } else count2++;
  
  //-- Histograms
  if (!strcmp(kv[0], "doLocalNoise"))   {pkPar->pkParTable[count2] = 1; pkPar->doLocalNoise = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "N_nu"))           {pkPar->pkParTable[count2] = 1; pkPar->N_nu         = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "bin_nu")) {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->bin_nu) free(pkPar->bin_nu); 
                                         pkPar->bin_nu = makeDoubleArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  
  //-- Outputs
  if (!strcmp(kv[0], "prefix"))         {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->prefix, kv, count); return 0;} else count2++;
  if (!strcmp(kv[0], "doFITS"))         {pkPar->pkParTable[count2] = 1; pkPar->doFITS        = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outHaloCat"))     {pkPar->pkParTable[count2] = 1; pkPar->outHaloCat    = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outGalCat"))      {pkPar->pkParTable[count2] = 1; pkPar->outGalCat     = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outMaps"))        {pkPar->pkParTable[count2] = 1; pkPar->outMaps       = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outTruth"))       {pkPar->pkParTable[count2] = 1; pkPar->outTruth      = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outMask"))        {pkPar->pkParTable[count2] = 1; pkPar->outMask       = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outPeakList"))    {pkPar->pkParTable[count2] = 1; pkPar->outPeakList   = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outHist"))        {pkPar->pkParTable[count2] = 1; pkPar->outHist       = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "outMultiscale"))  {pkPar->pkParTable[count2] = 1; pkPar->outMultiscale = atoi(kv[1]); return 0;} else count2++;
  
  //-- ABC
  if (!strcmp(kv[0], "ABC_f"))          {pkPar->pkParTable[count2] = 1; pkPar->ABC_f = atoi(kv[1]); return 0;} else count2++;
  if (!strcmp(kv[0], "ABC_doParam")) {
                                         pkPar->pkParTable[count2] = 1; 
                                         if (pkPar->ABC_doParam) free(pkPar->ABC_doParam); 
                                         pkPar->ABC_doParam = makeIntArray(kv, count, err); forwardError(*err, __LINE__, -1);
                                         return 0;
  } else count2++;
  if (!strcmp(kv[0], "ABC_Q"))          {pkPar->pkParTable[count2] = 1; pkPar->ABC_Q       = atoi(kv[1]);                         return 0;} else count2++;
  if (!strcmp(kv[0], "ABC_r_stop"))     {pkPar->pkParTable[count2] = 1; pkPar->ABC_r_stop  = atof(kv[1]);                         return 0;} else count2++;
  if (!strcmp(kv[0], "ABC_obsPath"))    {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->ABC_obsPath, kv, count);    return 0;} else count2++;
  if (!strcmp(kv[0], "ABC_doCorr"))     {pkPar->pkParTable[count2] = 1; pkPar->ABC_doCorr  = atoi(kv[1]);                         return 0;} else count2++;
  if (!strcmp(kv[0], "ABC_invCovPath")) {pkPar->pkParTable[count2] = 1; setPathWhichCanBeBlank(pkPar->ABC_invCovPath, kv, count); return 0;} else count2++;
  
  //-- Running part
  if (!strcmp(kv[0], "doSurvey"))            pkPar->doSurvey       = atoi(kv[1]);
  else if (!strcmp(kv[0], "realizationInd")) pkPar->realizationInd = atoi(kv[1]);
  else if (!strcmp(kv[0], "MPISize"))        pkPar->MPISize        = atoi(kv[1]);
  else if (!strcmp(kv[0], "MPIInd"))         pkPar->MPIInd         = atoi(kv[1]);
  else return 1;
  
  return 0;
}

void read_peak_param(char name[], peak_param *pkPar, error **err)
{
  //-- pkPar should have been initialized.
  
  sprintf(pkPar->pkParPath, "%s", name);
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  char line[STRING_LENGTH_MAX], kv[STRING_LENGTH_MAX][256], *buffer;
  int count, unknown;
  
  //-- Read
  do {
    buffer = fgets(line, STRING_LENGTH_MAX, file);
    ignoreComments(line);
    count = getKeyAndValues(line, kv);
    unknown = findKey_peak_param(pkPar, kv, count, err); forwardError(*err, __LINE__,);
    if (unknown == 1) printf("Detected the keyword \"%s\" but unable to interpret, skipped\n", kv[0]);
  }
  while (buffer != NULL);
  
  fclose(file);
  return;
}

int updateFromCommandLine(int argc, char *argv[], cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  //-- Update from the command line
  
  char kv[STRING_LENGTH_MAX][256];
  int i, count, unknown;
  
  for (i=3; i<argc; i++) {
    count = getKeyAndValues(argv[i], kv);
    if (count == -1) break; //-- Help detected
    unknown  = findKey_cosmo_hm(chPar, pkPar, kv, count, err); forwardError(*err, __LINE__, -1);
    unknown *= findKey_peak_param(pkPar, kv, count, err);      forwardError(*err, __LINE__, -1);
    if (unknown == 1 && count > 1) printf("Detected the keyword \"%s\" but unable to interpret, skipped\n", kv[0]);
  }
  
  return (int)(count == -1);
}

#define nbCmPar 16
#define nbHmPar 7
#define nbPkPar 72
void checkParam(peak_param *pkPar, error **err)
{
  
  char *cmParStr[nbCmPar] = {
    "Omega_m", "Omega_de", "w0_de", "w1_de", "h_100", "Omega_b", "Omega_nu_mass", "Neff_nu_mass", "normalization", "n_spec", 
    "snonlinear", "stransfer", "sgrowth", "sde_param", "normmode", "a_min"
  };
  
  char *hmParStr[nbHmPar] = {
    "cosmo_file", "par_nz", "alpha_NFW", "c0", "beta_NFW", "smassfct", "shalo_bias"
  };
  
  char *pkParStr[nbPkPar] = {
    "hmParPath", "seed", "verbose", 
    "field", "Omega", "theta_pix", "nside", "patch", "rotAng", "HP_resol", 
    "inHaloCatPath", "inHaloCatCol", "z_halo_min", "z_halo_max", "N_z_halo", "M_min", "M_max", "dlogM", 
    "inGalCatPath", "inGalCatCol", "z_s", "dz_gal", "doRandGalPos", "n_gal", "sigma_eps", 
    "doMask", "maskPath", "nbHoleTypes", "holeRadius", "holeDensity", "stripeLRatio", "stripeWRatio", 
    "doLensing", "doKappa", "doSubtraction", 
    "doSmoothing", "FFT_nbFilters", "FFT_filter", "FFT_scale", "DC_nbFilters", "DC_filter", "DC_scale",
    "doLocalNoise", "N_nu", "bin_nu", 
    "prefix", "doFITS", "outHaloCat", "outGalCat", "outMaps", "outTruth", "outMask", "outPeakList", "outHist", "outMultiscale",
    "ABC_f", "ABC_doParam", "ABC_Q", "ABC_r_stop", "ABC_obsPath", "ABC_doCorr", "ABC_invCovPath",
    "", "", "", "", "", "", "", "", "", ""
  };

  int hasError = 0;
  int i;
  
  for (i=0; i<nbCmPar; i++) {
    if (strcmp(cmParStr[i], "") && pkPar->cmParTable[i] < 1) {
      printf("Cosmological parameter \"%s\" in not defined!\n", cmParStr[i]);
      hasError++;
    }
  }
  
  for (i=0; i<nbHmPar; i++) {
    if (strcmp(hmParStr[i], "") && pkPar->hmParTable[i] < 1) {
      printf("Halo model parameter \"%s\" not defined!\n", hmParStr[i]);
      hasError++;
    }
  }
  
  for (i=0; i<nbPkPar; i++) {
    if (strcmp(pkParStr[i], "") && pkPar->pkParTable[i] < 1) {
      printf("Peak parameter \"%s\" not defined!\n", pkParStr[i]);
      hasError++;
    }
  }
  
  testErrorRet(hasError>0, peak_unknown, "Some parameters not defined", *err, __LINE__,);
  return;
}
#undef nbCmPar
#undef nbHmPar
#undef nbPkPar

void set_peak_param(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  //-- Precomputed part - Parameter files
  u_int32_t seed   = (strchr(pkPar->seed, 'r') == NULL) ? strtoul(pkPar->seed, NULL, 10) : renewSeed();
  pkPar->generator = initializeGenerator(seed);
  
  //-- Precomputed part - Field
  if (pkPar->field == 2)      pkPar->theta_pix = 1.0 * ARCMIN_TO_RADIAN; //-- [rad] Forced to 1.0 if field = 2
  if (pkPar->field != 1)      pkPar->rotAng    = 0.0 * ARCMIN_TO_RADIAN; //-- [rad] Forced to 0.0 if field != 1
  if      (pkPar->field == 0) pkPar->area = pkPar->Omega[0] * pkPar->Omega[1];
  else if (pkPar->field <= 2) pkPar->area = FULL_SKY * DEGREE_SQ_TO_RADIAN_SQ / (12.0 * SQ((double)pkPar->nside)); //-- [rad^2] area is related to sampling. For field = 1, it is the HEALPix patch area.
  else {
    testErrorRet(1, peak_unknown, "Unknown field type", *err, __LINE__,);
  }
  pkPar->theta_pix_inv = 1.0 / pkPar->theta_pix;
  if (pkPar->field < 2) {
    pkPar->resol[0] = (int)round(pkPar->Omega[0] / pkPar->theta_pix);
    pkPar->resol[1] = (int)round(pkPar->Omega[1] / pkPar->theta_pix);
  }
  else {
    pkPar->resol[0] = 6400 / pkPar->nside;
    pkPar->resol[1] = 6400 / pkPar->nside;
  }
  
  //-- Precomputed part - Halos
  pkPar->doInHaloCat = (!strcmp(pkPar->inHaloCatPath, "")) ? 0 : 1;
  pkPar->dz_halo     = (pkPar->z_halo_max - pkPar->z_halo_min) / (double)(pkPar->N_z_halo);
  pkPar->N_M         = (int)round(log10(pkPar->M_max / pkPar->M_min) / pkPar->dlogM); //-- Number of mass bins
  
  //-- Precomputed part - Galaxies
  pkPar->doInGalCat  = (!strcmp(pkPar->inGalCatPath, "")) ? 0 : 1;
  
  if (pkPar->z_s > 0.0) {
    pkPar->w_s  = w(chPar->cosmo, 1.0/(1.0 + pkPar->z_s), 0, err); forwardError(*err, __LINE__,); //-- wOmegar = 0
    pkPar->D_s  = f_K(chPar->cosmo, pkPar->w_s, err);              forwardError(*err, __LINE__,);
    pkPar->D_s /= 1.0 + pkPar->z_s;
  }
  else {
    pkPar->w_s  = -1.0; //-- w_s to be calculated in set_gal_t
    pkPar->D_s  = -1.0; //-- D_s to be calculated in set_gal_t
    if (pkPar->doRandGalPos == 0 && pkPar->verbose < 99) printf("Found z_s <= 0, doRandGalPos forced to 1\n");
    pkPar->doRandGalPos = 1; //-- Forced to 1 if z_s < 0
  }
  pkPar->N_z_gal    = (int)round((chPar->redshift->par_nz[1] - chPar->redshift->par_nz[0]) / pkPar->dz_gal);
  pkPar->doNoise    = (pkPar->sigma_eps > 0.0) ? 1 : 0;
  pkPar->sigma_half = pkPar->sigma_eps * SQRT_2_INV;
  pkPar->sigma_pix  = pkPar->sigma_half / sqrt(pkPar->n_gal * SQ(pkPar->theta_pix));
  
  //-- Precomputed part - Lensing
  if (pkPar->doInGalCat == 1 && (pkPar->inGalCatCol[4] > -1 || pkPar->inGalCatCol[5] > -1 || pkPar->inGalCatCol[6] > -1)) {
    if (pkPar->doLensing == 1 && pkPar->verbose < 99) printf("Found inGalCatPath not blank and the signal not skipped, doLensing forced to 0\n");
    pkPar->doLensing = 0; //-- Forced to 0 if inGalCatCol does not skip any of kappa, e_1, or e_2
  }
  
  //-- Precomputed part - Filters
  if (pkPar->field == 2) {
    if (pkPar->FFT_nbFilters > 0 && pkPar->verbose < 99) printf("Found field = 2, FFT_nbFilters forced to 0\n");
    pkPar->FFT_nbFilters = 0;
  }
  if (pkPar->field == 1) {
    if (pkPar->DC_nbFilters > 0 && pkPar->verbose < 99) printf("Found field = 1, DC_nbFilters forced to 0\n");
    pkPar->DC_nbFilters = 0;
  }
  if ((pkPar->doSmoothing / 1) % 2 == 0) {
    if (pkPar->FFT_nbFilters > 0 && pkPar->verbose < 99) printf("Found doSmoothing without bit 1, FFT_nbFilters forced to 0\n");
    pkPar->FFT_nbFilters = 0;
  }
  if ((pkPar->doSmoothing / 2) % 2 == 0) {
    if (pkPar->DC_nbFilters > 0 && pkPar->verbose < 99) printf("Found doSmoothing without bit 2, DC_nbFilters forced to 0\n");
    pkPar->DC_nbFilters = 0;
  }
  pkPar->nbFilters    = MAX(pkPar->FFT_nbFilters + pkPar->DC_nbFilters, 1); //-- Keep at least one
  pkPar->smootherSize = MAX(pkPar->FFT_nbFilters, pkPar->DC_nbFilters);
  pkPar->smootherSize = MAX(pkPar->smootherSize, 1); //-- Keep at least one table for inversion
  pkPar->filter       = (int*)malloc_err(pkPar->nbFilters * sizeof(int), err); forwardError(*err, __LINE__,);
  int count = 0;
  int i;
  for (i=0; i<pkPar->FFT_nbFilters; i++, count++) pkPar->filter[count] = pkPar->FFT_filter[i];
  for (i=0; i<pkPar->DC_nbFilters; i++, count++)  pkPar->filter[count] = pkPar->DC_filter[i];
  
  //-- Precomputed part - Filters - FFT
  pkPar->FFT_firstToInvert   = pkPar->FFT_nbFilters;
  pkPar->FFT_firstToConserve = pkPar->FFT_nbFilters;
  pkPar->FFT_hasGammaT       = 0;
  for (i=0; i<pkPar->FFT_nbFilters; i++) {
    if (pkPar->FFT_filter[i] < 2)  pkPar->FFT_firstToInvert   = MIN(pkPar->FFT_firstToInvert, i);
    if (pkPar->FFT_filter[i] >= 2) pkPar->FFT_firstToConserve = MIN(pkPar->FFT_firstToConserve, i);
    if (pkPar->FFT_filter[i] == 3) pkPar->FFT_hasGammaT       = 1;
  }
  pkPar->FFT_scaleInPix  = (pkPar->FFT_nbFilters == 0) ? NULL : (double*)malloc_err((pkPar->FFT_nbFilters+pkPar->FFT_hasGammaT) * sizeof(double), err); forwardError(*err, __LINE__,);
  pkPar->FFT_sigma_noise = (pkPar->FFT_nbFilters == 0) ? NULL : (double*)malloc_err(pkPar->FFT_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
  double bufferSize = 0.0;
  for (i=0; i<pkPar->FFT_nbFilters+pkPar->FFT_hasGammaT; i++) pkPar->FFT_scaleInPix[i] = pkPar->FFT_scale[i] * pkPar->theta_pix_inv;
  for (i=0; i<pkPar->FFT_nbFilters; i++) {
    //-- Let W(x) be a kernel with || W ||_L1 = NORM1 and || W ||_L2^2 = NORM2^2.
    //-- The scale-dependent kernel should be defined as W_s(i) = W(i/s) / s^2 with
    //--     || W_s ||_L1   = || W ||_L1         = NORM1
    //-- and || W_s ||_L2^2 = || W ||_L2^2 / s^2 = NORM2^2 / s^2.
    //--
    //-- Thus, sigma_noise^2 = (sigma_epsilon^2 / 2 n_gal) (|| W_s ||_L2 / || W_s ||_L1)^2
    //--                     = (sigma_half^2 / n_gal) (NORM2 / NORM1 s)^2
    //-- So, sigma_noise = (sigma_half / sqrt(n_gal) s) (NORM2 / NORM1)
    //--
    //-- Gaussian:     W(x) = exp(-x^2),    NORM1 = pi,             NORM2^2 = pi / 2, NORM2 / NORM1 = 1 / sqrt(2 pi)
    //-- Starlet:      W(x, y) = psi(x, y), NORM1 = numerical,      NORM2^2 = 5 I2^2 - 2 I3^2
    //-- M_ap tanh:    W(x) = tanh(x / 0.1) / [x (1 + exp(5 - 150 x) + exp(-47 + 50 x))], 
    //--                                    NORM1 = numerical,      NORM2^2 = numerical
    //-- M_ap gamma_t: W(x) = 1_{a<=x<b},   NORM1 = pi (b^2 - a^2), NORM2^2 = pi (b^2 - a^2)
    if (pkPar->FFT_filter[i] == 0) {
      //-- Gaussian
      bufferSize = MAX(bufferSize, CUTOFF_FACTOR_GAUSSIAN * pkPar->FFT_scaleInPix[i]);
      pkPar->FFT_sigma_noise[i] = pkPar->sigma_half / (sqrt(pkPar->n_gal) * pkPar->FFT_scale[i]) * (1.0 / sqrt(TWO_PI));
    }
    else if (pkPar->FFT_filter[i] == 1) {
      //-- Starlet
      bufferSize = MAX(bufferSize, 2.0 * pkPar->FFT_scaleInPix[i]);
      pkPar->FFT_sigma_noise[i] = pkPar->sigma_half / (sqrt(pkPar->n_gal) * pkPar->FFT_scale[i]) * (NORM2_STARLET / NORM1_STARLET);
    }
    else if (pkPar->FFT_filter[i] == 2) {
      //-- M_ap tanh
      testErrorRet(pkPar->doKappa==1, peak_setting, "FFT_filter = 2 only allowed if doKappa != 1", *err, __LINE__,);
      bufferSize = MAX(bufferSize, CUTOFF_FACTOR_M_AP_TANH * pkPar->FFT_scaleInPix[i]);
      pkPar->FFT_sigma_noise[i] = pkPar->sigma_half / (sqrt(pkPar->n_gal) * pkPar->FFT_scale[i]) * (NORM2_M_AP_TANH / NORM1_M_AP_TANH);
    }
    else if (pkPar->FFT_filter[i] == 3) {
      //-- M_ap gamma_t
      testErrorRet(pkPar->doKappa==1, peak_setting, "FFT_filter = 3 only allowed if doKappa != 1", *err, __LINE__,);
      bufferSize = MAX(bufferSize, pkPar->FFT_scaleInPix[i]);
      bufferSize = MAX(bufferSize, pkPar->FFT_scaleInPix[i+1]);
      pkPar->FFT_sigma_noise[i] = pkPar->sigma_half / sqrt(pkPar->n_gal * PI * (SQ(pkPar->FFT_scale[i+1]) - SQ(pkPar->FFT_scale[i])));
    }
    else {
      testErrorRet(1, peak_unknown, "Invalid filter type for FFT", *err, __LINE__,);
    }
  }
  
  //-- Precomputed part - Filters - DC
  pkPar->DC_hasGammaT = 0;
  for (i=0; i<pkPar->DC_nbFilters; i++) {
    if (pkPar->DC_filter[i] == 3) {
      pkPar->DC_hasGammaT = 1;
      break;
    }
  }
  pkPar->DC_scaleInPix  = (pkPar->DC_nbFilters == 0) ? NULL : (double*)malloc_err((pkPar->DC_nbFilters+pkPar->DC_hasGammaT) * sizeof(double), err); forwardError(*err, __LINE__,);
  pkPar->DC_scale_inv   = (pkPar->DC_nbFilters == 0) ? NULL : (double*)malloc_err((pkPar->DC_nbFilters+pkPar->DC_hasGammaT) * sizeof(double), err); forwardError(*err, __LINE__,);
  pkPar->DC_cut         = (pkPar->DC_nbFilters == 0) ? NULL : (double*)malloc_err((pkPar->DC_nbFilters+pkPar->DC_hasGammaT) * sizeof(double), err); forwardError(*err, __LINE__,);
  pkPar->DC_sigma_noise = (pkPar->DC_nbFilters == 0) ? NULL : (double*)malloc_err(pkPar->DC_nbFilters * sizeof(double), err); forwardError(*err, __LINE__,);
  for (i=0; i<pkPar->DC_nbFilters+pkPar->DC_hasGammaT; i++) pkPar->DC_scaleInPix[i] = pkPar->DC_scale[i] * pkPar->theta_pix_inv;
  for (i=0; i<pkPar->DC_nbFilters; i++) {
    //-- See DC_sigma_noise for notes on sigma_noise
    if (pkPar->DC_filter[i] == 0) {
      //-- Gaussian
      testErrorRet(pkPar->doKappa!=1, peak_setting, "DC_filter = 0 only allowed if doKappa = 1", *err, __LINE__,);
      bufferSize = MAX(bufferSize, CUTOFF_FACTOR_GAUSSIAN * pkPar->DC_scaleInPix[i]);
      pkPar->DC_scale_inv[i]   = 1.0 / (SQ(pkPar->DC_scale[i]));                  //-- Squared for the Gaussian
      pkPar->DC_cut[i]         = SQ(CUTOFF_FACTOR_GAUSSIAN * pkPar->DC_scale[i]); //-- Squared for the Gaussian
      pkPar->DC_sigma_noise[i] = pkPar->sigma_half / (sqrt(pkPar->n_gal) * pkPar->DC_scale[i]) * (1.0 / sqrt(TWO_PI));
    }
    else if (pkPar->DC_filter[i] == 1) {
      //-- Starlet
      testErrorRet(pkPar->doKappa!=1, peak_setting, "DC_filter = 1 only allowed if doKappa = 1", *err, __LINE__,);
      bufferSize = MAX(bufferSize, 2.0 * pkPar->DC_scaleInPix[i]);
      pkPar->DC_scale_inv[i]   = 1.0 / pkPar->DC_scale[i];
      pkPar->DC_cut[i]         = 2.0 * pkPar->DC_scale[i];
      pkPar->DC_sigma_noise[i] = pkPar->sigma_half / (sqrt(pkPar->n_gal) * pkPar->DC_scale[i]) * (NORM2_STARLET / NORM1_STARLET);
    }
    else if (pkPar->DC_filter[i] == 2) {
      //-- M_ap tanh
      testErrorRet(pkPar->doKappa==1, peak_setting, "DC_filter = 2 only allowed if doKappa != 1", *err, __LINE__,);
      bufferSize = MAX(bufferSize, CUTOFF_FACTOR_M_AP_TANH * pkPar->DC_scaleInPix[i]);
      pkPar->DC_scale_inv[i]   = 1.0 / pkPar->DC_scale[i];
      pkPar->DC_cut[i]         = CUTOFF_FACTOR_M_AP_TANH * pkPar->DC_scale[i];
      pkPar->DC_sigma_noise[i] = pkPar->sigma_half / (sqrt(pkPar->n_gal) * pkPar->DC_scale[i]) * (NORM2_M_AP_TANH / NORM1_M_AP_TANH);
    }
    else if (pkPar->DC_filter[i] == 3) {
      //-- M_ap gamma_t
      testErrorRet(pkPar->doKappa==1, peak_setting, "DC_filter = 3 only allowed if doKappa != 1", *err, __LINE__,);
      bufferSize = MAX(bufferSize, pkPar->DC_scaleInPix[i]);
      bufferSize = MAX(bufferSize, pkPar->DC_scaleInPix[i+1]);
      pkPar->DC_scale_inv[i]   = 1.0 / pkPar->DC_scale[i];
      pkPar->DC_scale_inv[i+1] = 1.0 / pkPar->DC_scale[i+1];
      pkPar->DC_cut[i]         = pkPar->DC_scale[i];
      pkPar->DC_cut[i+1]       = pkPar->DC_scale[i+1];
      pkPar->DC_sigma_noise[i] = pkPar->sigma_half / sqrt(pkPar->n_gal * PI * (SQ(pkPar->DC_scale[i+1]) - SQ(pkPar->DC_scale[i])));
    }
    else {
      testErrorRet(1, peak_unknown, "Invalid filter type for DC", *err, __LINE__,);
    }
  }
  
  //-- Precomputed part - Maps & histograms
  if (pkPar->field < 2) {
    pkPar->bufferSize        = (int)ceil(bufferSize);                 //-- Buffer size for filter
    pkPar->FFTSize           = MAX(pkPar->resol[0], pkPar->resol[1]); //-- FFT size with buffer area
    pkPar->FFTSize          += 256 - (pkPar->FFTSize % 256);          //-- Round the size to a multiple of 256
    pkPar->FFTNormFactor     = 1.0 / (double)(SQ(pkPar->FFTSize));
    pkPar->peakFieldResol[0] = pkPar->resol[0] - 2 * pkPar->bufferSize;
    pkPar->peakFieldResol[1] = pkPar->resol[1] - 2 * pkPar->bufferSize;
    pkPar->peakListLength    = (pkPar->peakFieldResol[0] + 1) * (pkPar->peakFieldResol[1] + 1) / 4 + 1;
  }
  else {
    pkPar->bufferSize        = (int)ceil(bufferSize);                 //-- Buffer size for filter
    pkPar->FFTSize           = pkPar->HP_resol;
    pkPar->FFTNormFactor     = 1.0 / (double)(SQ(pkPar->FFTSize));
    pkPar->peakFieldResol[0] = pkPar->HP_resol - 2;
    pkPar->peakFieldResol[1] = pkPar->HP_resol - 2;
    pkPar->peakListLength    = (pkPar->peakFieldResol[0] + 1) * (pkPar->peakFieldResol[1] + 1) / 4 + 1;
  }
  
  //-- Precomputed part - Only used if field = 1 (projected HEALPix) or 2 (HEALPix)
  if (pkPar->field == 1 || pkPar->field == 2) {
#ifdef __CAMELUS_USE_HEALPIX__
    testErrorRet(pkPar->nside<8, peak_badValue, "nside too small", *err, __LINE__,);
    ring2nest(pkPar->nside, pkPar->patch, &pkPar->nest);
    pkPar->nsidePix  = pkPar->nside * pkPar->HP_resol;
    pkPar->HP_length = pkPar->HP_resol * pkPar->HP_resol;
    pkPar->HP_first  = pkPar->nest * pkPar->HP_length;
    decompose(pkPar->nside, pkPar->patch, 0, &pkPar->cap, &pkPar->level, &pkPar->length, &pkPar->off, &pkPar->j); //-- nest = 0
    
    pix2ang_ring(pkPar->nside, pkPar->patch, &pkPar->center[1], &pkPar->center[0]); //-- theta, phi in [rad]
    pkPar->z0        = cos(pkPar->center[1]);      //-- Compute z0 = cos(theta)
    pkPar->center[1] = HALF_PI - pkPar->center[1]; //-- Convert theta to DEC in [rad]
    pkPar->center[2] = sin(pkPar->center[1]);
    pkPar->center[3] = cos(pkPar->center[1]);
    if (pkPar->field == 1) pkPar->offset = 0.5 * MAX(pkPar->Omega[0], pkPar->Omega[1]);
    else                   pkPar->offset = 3200.0 * ARCMIN_TO_RADIAN / (double)pkPar->nside; //-- [rad] The field size set to 6400 / nside [arcmin].
#else
    testErrorRet(1, peak_setting, "Camelus not linked to chealpix, field = 1 or 2 unavailable", *err, __LINE__,);
#endif
  }
  
  //-- Precomputed part - Outputs
#ifdef __CAMELUS_USE_FITS__
#else
  if (pkPar->doFITS > 0 && pkPar->verbose < 99) printf("Camelus not linked to cfitsio, doFITS forced to 0\n");
  pkPar->doFITS = 0;
#endif
  
  //-- Running part
  pkPar->doSurvey       = 0;  //-- 0 = no multiplicative correction, 1 = yes (conflit with doKappa < 2)
  pkPar->realizationInd = 0;  //-- In use
  //pkPar->MPISize;  //-- Should not be initialized here
  //pkPar->MPIInd;   //-- Should not be initialized here
  return;
}

//----------------------------------------------------------------------
//-- Functions related to printing

char *printDoKappa(int doKappa)
{
  if (doKappa == 0) return "0 (gamma)";
  if (doKappa == 1) return "1 (kappa)";
  if (doKappa == 2) return "2 (g with linear KS)";
  if (doKappa == 3) return "3 (g with iterative KS)";
  if (doKappa == 4) return "4 (g with SS)";
  return " ";
}

char *printDoSubtraction(int doSubtraction)
{
  if (doSubtraction == 0) return "0 (without)";
  if (doSubtraction == 1) return "1 (subtract mean)";
  if (doSubtraction == 2) return "2 (subtract mass sheet)";
  return " ";
}

char *printField(int field)
{
  if (field == 0) return "0 (rectangle)";
  if (field == 1) return "1 (projected HEALPix)";
  if (field == 2) return "2 (HEALPix)";
  return " ";
}

char *printFilter(int filter)
{
  if (filter == 0) return "0 (Gaussian)";
  if (filter == 1) return "1 (starlet)";
  if (filter == 2) return "2 (M_ap tanh)";
  if (filter == 3) return "3 (M_ap gamma_t)";
  return " ";
}

void printIntArray(int *iArr, int length)
{
  if (length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  printf("%d", iArr[0]);
  for (i=1; i<length; i++) printf(", %d", iArr[i]);
  return;
}

void printDoubleArray(double *lfArr, int length, double factor, int digit)
{
  if (lfArr == NULL || length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  if (digit == 0) {
    printf("%f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %f", lfArr[i]*factor);
  }
  else if (digit == 1) {
    printf("%.1f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.1f", lfArr[i]*factor);
  }
  else if (digit == 2) {
    printf("%.2f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.2f", lfArr[i]*factor);
  }
  else if (digit == 3) {
    printf("%.3f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.3f", lfArr[i]*factor);
  }
  else if (digit == 4) {
    printf("%.4f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.4f", lfArr[i]*factor);
  }
  else if (digit == 5) {
    printf("%.5f", lfArr[0]*factor);
    for (i=1; i<length; i++) printf(", %.5f", lfArr[i]*factor);
  }
  return;
}

void printGalaxyInfo(peak_param *pkPar, redshift_t *sdPar)
{
  if (pkPar->doRandGalPos == 0) printf("Galaxy redshift fixed at %.3f, on a regular grid (%dx%d)\n", pkPar->z_s, pkPar->resol[0], pkPar->resol[1]);
  else if (pkPar->z_s > 0)      printf("Galaxy redshift fixed at %.3f, random angular position\n", pkPar->z_s);
  else {
    printf("z_gal_min  z_gal_max  alpha_gal  beta_gal  z_gal_0\n");
    printf("    %5.3f      %5.3f      %5.3f     %5.3f    %5.3f\n", sdPar->par_nz[0], sdPar->par_nz[1], sdPar->par_nz[2], sdPar->par_nz[3], sdPar->par_nz[4]);
  }
  return;
}

void printFilterArray(int *iArr, int length)
{
  if (iArr == NULL || length == 0) {
    printf("NULL");
    return;
  }
  
  int i;
  printf("%s", printFilter(iArr[0]));
  for (i=1; i<length; i++) printf(", %s", printFilter(iArr[i]));
  return;
}

void printParam(cosmo_hm *chPar, peak_param *pkPar)
{
  printf("Cosmological parameters\n");
  printf("Omega_m  Omega_de  Omega_b    n_s  h_100  sigma_8  w0_de  w1_de\n");
  printf(" %6.4f    %6.4f   %6.4f  %5.3f  %5.3f    %5.3f  %5.2f   %4.2f\n",
	  chPar->cosmo->Omega_m, chPar->cosmo->Omega_de, chPar->cosmo->Omega_b, chPar->cosmo->n_spec, chPar->cosmo->h_100,
	  chPar->cosmo->normalization, chPar->cosmo->w0_de, chPar->cosmo->w1_de);
  printf("\n");
  printf("Galaxy parameters\n");
  printGalaxyInfo(pkPar, chPar->redshift);
  printf("\n");
  
  printf("Peak parameters\n");
  printf("Omega         = (%.1f, %.1f) [arcmin]\n", pkPar->Omega[0]*RADIAN_TO_ARCMIN, pkPar->Omega[1]*RADIAN_TO_ARCMIN);
  printf("theta_pix     = %.3f [arcmin]\n", pkPar->theta_pix*RADIAN_TO_ARCMIN);
  printf("Resolution    = (%d, %d) [pix]\n", pkPar->resol[0], pkPar->resol[1]);
  printf("Buffer size   = %d [pix]\n", pkPar->bufferSize);
  printf("Peak field    = (%d, %d) [pix]\n", pkPar->peakFieldResol[0], pkPar->peakFieldResol[1]);
  printf("field         = %s\n", printField(pkPar->field));
  printf("doKappa       = %d\n", pkPar->doKappa);
  printf("doSubtraction = %d\n", pkPar->doSubtraction);
  
  if (pkPar->FFT_nbFilters) {
  printf("FFT_filter    = "); printFilterArray(pkPar->FFT_filter, pkPar->FFT_nbFilters); printf("\n");
  printf("FFT_scale     = "); printDoubleArray(pkPar->FFT_scale, pkPar->FFT_nbFilters, RADIAN_TO_ARCMIN, 3); printf(" [arcmin]\n");
  }
  if (pkPar->DC_nbFilters) {
  printf("DC_filter     = "); printFilterArray(pkPar->DC_filter, pkPar->DC_nbFilters); printf("\n");
  printf("DC_scale      = "); printDoubleArray(pkPar->DC_scale, pkPar->DC_nbFilters, RADIAN_TO_ARCMIN, 3); printf(" [arcmin]\n");
  }
  return;
}

void printParam_peakCustomized_1(peak_param *pkPar)
{
  printf("Peak parameters, customized part\n");
  printf("------ Parameter files ------\n");
  printf("pkParPath       = \"%s\"\n", pkPar->pkParPath);
  printf("hmParPath       = \"%s\"\n", pkPar->hmParPath);
  printf("cmParPath       = \"%s\"\n", pkPar->cmParPath);
  printf("seed            = %s\n", pkPar->seed);
  printf("verbose         = %d\n", pkPar->verbose);
  printf("------ Field ------\n");
  printf("field           = %s\n", printField(pkPar->field));
  printf("Omega           = (%.1f, %.1f) [arcmin]\n", pkPar->Omega[0]*RADIAN_TO_ARCMIN, pkPar->Omega[1]*RADIAN_TO_ARCMIN);
  printf("theta_pix       = %.3f [arcmin]\n", pkPar->theta_pix*RADIAN_TO_ARCMIN);
  printf("nside           = %ld\n", pkPar->nside);
  printf("patch           = %ld\n", pkPar->patch);
  printf("rotAng          = %f\n", pkPar->rotAng);
  printf("HP_resol        = %d\n", pkPar->HP_resol);
  printf("------ Halos ------\n");
  printf("inHaloCatPath   = \"%s\"\n", pkPar->inHaloCatPath);
  printf("inHaloCatCol    = %d, %d, %d, %d, %d\n", pkPar->inHaloCatCol[0], pkPar->inHaloCatCol[1], pkPar->inHaloCatCol[2], pkPar->inHaloCatCol[3], pkPar->inHaloCatCol[4]);
  printf("z_halo_min      = %.3f [-]\n", pkPar->z_halo_min);
  printf("z_halo_max      = %.3f [-]\n", pkPar->z_halo_max);
  printf("N_z_halo        = %d\n", pkPar->N_z_halo);
  printf("M_min           = %.3e [M_sol/h]\n", pkPar->M_min);
  printf("M_max           = %.3e [M_sol/h]\n", pkPar->M_max);
  printf("dlogM           = %.3f [-]\n", pkPar->dlogM);
  printf("------ Galaxies ------\n");
  printf("inGalCatPath    = \"%s\"\n", pkPar->inGalCatPath);
  printf("inGalCatCol     = %d, %d, %d, %d, %d, %d, %d\n", 
	 pkPar->inGalCatCol[0], pkPar->inGalCatCol[1], pkPar->inGalCatCol[2], pkPar->inGalCatCol[3], pkPar->inGalCatCol[4], pkPar->inGalCatCol[5], pkPar->inGalCatCol[6]);
  printf("z_s             = %.3f [-]\n", pkPar->z_s);
  printf("dz_gal          = %.3f [-]\n", pkPar->dz_gal);
  printf("doRandGalPos    = %d\n", pkPar->doRandGalPos);
  printf("n_gal           = %.3f [arcmin^-2]\n", pkPar->n_gal/RADIAN_SQ_TO_ARCMIN_SQ);
  printf("sigma_eps       = %.3f [-]\n", pkPar->sigma_eps);
  printf("------ Masks ------\n");
  printf("doMask          = %d\n", pkPar->doMask);
  printf("maskPath        = \"%s\"\n", pkPar->maskPath);
  printf("nbHoleTypes     = %d\n", pkPar->nbHoleTypes);
  printf("holeRadius      = "); printDoubleArray(pkPar->holeRadius, pkPar->nbHoleTypes, RADIAN_TO_ARCMIN, 3); printf(" [arcmin]\n");
  printf("holeDensity     = "); printDoubleArray(pkPar->holeDensity, pkPar->nbHoleTypes, 1.0/RADIAN_SQ_TO_DEGREE_SQ, 3); printf(" [deg^-2]\n");
  printf("stripeLRatio    = "); printDoubleArray(pkPar->stripeLRatio, pkPar->nbHoleTypes, 1.0, 3); printf("\n");
  printf("stripeWRatio    = "); printDoubleArray(pkPar->stripeWRatio, pkPar->nbHoleTypes, 1.0, 3); printf("\n");
  return;
}

void printParam_peakCustomized_2(peak_param *pkPar)
{
  printf("------ Lensing ------\n");
  printf("doLensing       = %d\n", pkPar->doLensing);
  printf("doKappa         = %s\n", printDoKappa(pkPar->doKappa));
  printf("doSubtraction   = %s\n", printDoSubtraction(pkPar->doSubtraction));
  printf("------ Filter ------\n");
  printf("doSmoothing     = %d\n", pkPar->doSmoothing);
  printf("FFT_nbFilters   = %d\n", pkPar->FFT_nbFilters);
  printf("FFT_filter      = "); printFilterArray(pkPar->FFT_filter, pkPar->FFT_nbFilters); printf("\n");
  printf("FFT_scale       = "); printDoubleArray(pkPar->FFT_scale, pkPar->FFT_nbFilters+pkPar->FFT_hasGammaT, RADIAN_TO_ARCMIN, 3); printf(" [arcmin]\n");
  printf("DC_nbFilters    = %d\n", pkPar->DC_nbFilters);
  printf("DC_filter       = "); printFilterArray(pkPar->DC_filter, pkPar->DC_nbFilters); printf("\n");
  printf("DC_scale        = "); printDoubleArray(pkPar->DC_scale, pkPar->DC_nbFilters+pkPar->DC_hasGammaT, RADIAN_TO_ARCMIN, 3); printf(" [arcmin]\n");
  printf("------ Peak historgram ------\n");
  printf("doLocalNoise    = %d\n", pkPar->doLocalNoise);
  printf("N_nu            = %d\n", pkPar->N_nu);
  printf("bin_nu          = "); printDoubleArray(pkPar->bin_nu, pkPar->N_nu+1, 1.0, 3); printf(" [-]\n");
  printf("------ Outputs ------\n");
  printf("prefix          = \"%s\"\n", pkPar->prefix);
  printf("doFITS          = %d\n", pkPar->doFITS);
  printf("outHaloCat      = %d\n", pkPar->outHaloCat);
  printf("outGalCat       = %d\n", pkPar->outGalCat);
  printf("outMaps         = %d\n", pkPar->outMaps);
  printf("outTruth        = %d\n", pkPar->outTruth);
  printf("outMask         = %d\n", pkPar->outMask);
  printf("outPeakList     = %d\n", pkPar->outPeakList);
  printf("outHist         = %d\n", pkPar->outHist);
  printf("outMultiscale   = %d\n", pkPar->outMultiscale);
  printf("------ ABC ------\n");
  printf("ABC_f           = %d\n", pkPar->ABC_f);
  printf("ABC_doParam     = "); printIntArray(pkPar->ABC_doParam, pkPar->ABC_f); printf("\n");
  printf("ABC_Q           = %d\n", pkPar->ABC_Q);
  printf("ABC_r_stop      = %f\n", pkPar->ABC_r_stop);
  printf("ABC_obsPath     = \"%s\"\n", pkPar->ABC_obsPath);
  printf("ABC_doCorr      = %d\n", pkPar->ABC_doCorr);
  printf("ABC_invCovPath  = \"%s\"\n", pkPar->ABC_invCovPath);
  return;
}

void printParam_peakPrecomputed(peak_param *pkPar)
{
  printf("Peak parameters, default part\n");
  printf("------ Halos ------\n");
  printf("doInHaloCat     = %d\n", pkPar->doInHaloCat);
  printf("dz_halo         = %.3f [-]\n", pkPar->dz_halo);
  printf("N_M             = %d\n", pkPar->N_M);
  printf("------ Galaxies ------\n");
  printf("doInGalCat      = %d\n", pkPar->doInGalCat);
  printf("w_s             = %.3f [Mpc/h]\n", pkPar->w_s);
  printf("D_s             = %.3f [Mpc/h]\n", pkPar->D_s);
  printf("N_z_gal         = %d\n", pkPar->N_z_gal);
  printf("doNoise         = %d\n", pkPar->doNoise);
  printf("sigma_half      = %.3f [-]\n", pkPar->sigma_half);
  printf("------ Field & maps ------\n");
  printf("area            = %.3f [arcmin^2]\n", pkPar->area * RADIAN_SQ_TO_ARCMIN_SQ);
  printf("theta_pix_inv   = %.3f [arcmin^-1]\n", pkPar->theta_pix_inv / RADIAN_TO_ARCMIN);
  printf("sigma_pix       = %f [-]\n", pkPar->sigma_pix);
  printf("resol           = (%d, %d) [pix]\n", pkPar->resol[0], pkPar->resol[1]);
  printf("------ Filters ------\n");
  printf("nbFilters       = %d\n", pkPar->nbFilters);
  printf("smootherSize    = %d\n", pkPar->smootherSize);
  printf("filter          = "); printFilterArray(pkPar->filter, pkPar->nbFilters); printf("\n");
  printf("FFT_hasGammaT   = %d\n", pkPar->FFT_hasGammaT);
  printf("FFT_scaleInPix  = "); printDoubleArray(pkPar->FFT_scaleInPix, pkPar->FFT_nbFilters+pkPar->FFT_hasGammaT, 1.0, 3); printf(" [pix]\n");
  printf("FFT_sigma_noise = "); printDoubleArray(pkPar->FFT_sigma_noise, pkPar->FFT_nbFilters, 1.0, 0); printf(" [-]\n");
  printf("DC_hasGammaT    = %d\n", pkPar->DC_hasGammaT);
  printf("DC_scaleInPix   = "); printDoubleArray(pkPar->DC_scaleInPix, pkPar->DC_nbFilters+pkPar->DC_hasGammaT, 1.0, 3); printf(" [pix]\n");
  printf("DC_scale_inv    = "); printDoubleArray(pkPar->DC_scale_inv, pkPar->DC_nbFilters+pkPar->DC_hasGammaT, 1.0/RADIAN_SQ_TO_ARCMIN_SQ, 3); printf(" [arcmin^-2]\n");
  printf("DC_cut          = "); printDoubleArray(pkPar->DC_cut, pkPar->DC_nbFilters+pkPar->DC_hasGammaT, RADIAN_SQ_TO_ARCMIN_SQ, 3); printf(" [arcmin^2]\n");
  printf("DC_sigma_noise  = "); printDoubleArray(pkPar->DC_sigma_noise, pkPar->DC_nbFilters, 1.0, 0); printf(" [-]\n");
  printf("bufferSize      = %d [pix]\n", pkPar->bufferSize);
  printf("FFTSize         = %d [pix]\n", pkPar->FFTSize);
  printf("FFTNormFactor   = %.9f [-]\n", pkPar->FFTNormFactor);
  printf("peakListLength  = %d\n", pkPar->peakListLength);
  if (pkPar->field == 1 || pkPar->field == 2) {
  printf("------ Only used if field = 1 (projected HEALPix) or 2 (HEALPix) ------\n");
  printf("nest            = %ld\n", pkPar->nest);
  printf("nsidePix        = %ld\n", pkPar->nsidePix);
  printf("HP_length       = %d\n", pkPar->HP_length);
  printf("HP_first        = %ld\n", pkPar->HP_first);
  printf("cap             = %d\n", pkPar->cap);
  printf("level           = %d\n", pkPar->level);
  printf("length          = %d\n", pkPar->length);
  printf("off             = %d\n", pkPar->off);
  printf("j               = %d\n", pkPar->j);
  printf("z0              = %.3f\n", pkPar->z0);
  printf("center          = %.3f, %.3f, %.3e, %.3e [rad, rad, -, -]\n", pkPar->center[0], pkPar->center[1], pkPar->center[2], pkPar->center[3]);
  printf("offset          = %.1f [arcmin]\n", pkPar->offset*RADIAN_TO_ARCMIN);
  }
  printf("\n");
  printf("Peak parameters, running part\n");
  printf("doSurvey        = %d\n", pkPar->doSurvey);
  printf("realizationInd  = %d\n", pkPar->realizationInd);
  printf("MPISize         = %d\n", pkPar->MPISize);
  printf("MPIInd          = %d\n", pkPar->MPIInd);
  return;
}

//----------------------------------------------------------------------

