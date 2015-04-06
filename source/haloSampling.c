

  /*******************************
   **  haloSampling.c		**
   **  Chieh-An Lin		**
   **  Version 2015.03.02	**
   *******************************/


#include "haloSampling.h"


//----------------------------------------------------------------------
//-- Functions related to halo_t, halo_node, halo_list

void set_halo_t(cosmo_hm *cmhm, halo_t *h, double z, double M, double pos[2], error **err)
{
  testError(h==NULL, peak_null, "halo_t *h is NULL", *err, __LINE__); forwardError(*err, __LINE__,);
  h->z        = z;
  double a    = 1.0/(1.0+z);
  h->a        = a;
  int wOmegar = 0;
  h->w        = w(cmhm->cosmo, a, wOmegar, err); forwardError(*err, __LINE__,);
  double D_l  = a * f_K(cmhm->cosmo, h->w, err); forwardError(*err, __LINE__,);
  h->M        = M;
  double c    = concentration(cmhm, M, a, err);  forwardError(*err, __LINE__,);
  h->c        = c;
  h->c_sq     = SQ(c);
  h->f        = 1.0 / (log(1.0 + c) - c/(1.0 + c));
  if (pos != NULL) {
    h->pos[0] = pos[0];
    h->pos[1] = pos[1];
  }
  
  //-- r_vir is halo's physical size
  //-- M = (4 pi / 3) * rho_vir(z) * r_vir^3
  //--   = (4 pi / 3) * (rho_vir(z) / rho_m(z)) * rho_m(z) * r_vir^3
  //--   = (4 pi / 3) * Delta_vir(z) * rho_m(0) * (1+z)^3 * r_vir^3
  //--   = (4 pi / 3) * Delta_vir(z) * rho_crit(0) * Omega_m * (r_vir / a)^3
  h->r_vir        = a * cbrt(M / (FOUR_PI_OVER_THREE * Delta_vir(cmhm, a) * CRITICAL_DENSITY * cmhm->cosmo->Omega_m));
  h->r_vir_sq     = SQ(h->r_vir);
  //-- theta_vir = r_vir / D_A
  //-- given by angular diameter distance D_A = f_K / (1+z)
  h->theta_vir    = h->r_vir / D_l * RADIAN_TO_ARCMIN;
  h->theta_vir_sq = SQ(h->theta_vir);
  
  //-- factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  //-- See rayTracing.c for details
  h->factor = FOUR_PI_G_OVER_C2 * D_l * M * h->f * h->c_sq / (TWO_PI * h->r_vir_sq);
  return;
}

halo_node *initialize_halo_node(error **err)
{
  halo_node *hNode = (halo_node*)malloc_err(sizeof(halo_node), err); forwardError(*err, __LINE__, NULL);
  hNode->h         = (halo_t*)malloc_err(sizeof(halo_t), err);       forwardError(*err, __LINE__, NULL);;
  hNode->next      = NULL;
  return hNode;
}

halo_list *initialize_halo_list(error **err)
{
  halo_list *hList = (halo_list*)malloc_err(sizeof(halo_list), err); forwardError(*err, __LINE__, NULL);
  hList->length    = 0;
  hList->size      = 0;
  hList->first     = NULL;
  hList->last      = NULL;
  return hList;
}

void free_halo_list(halo_list *hList)
{
  halo_node *hNode;
  while (hList->first != NULL) {
    hNode        = hList->first;
    hList->first = hNode->next;
    if (hNode->h) {free(hNode->h); hNode->h = NULL;}
    free(hNode); hNode = NULL;
  }
  free(hList); hList = NULL;
  return;
}

void reset_halo_list(halo_list *hList)
{
  hList->size = 0;
  return;
}

void append_halo_list(cosmo_hm *cmhm, halo_list *hList, double z, double M, double pos[2], error **err)
{
  halo_t *h;
  if (hList->length == 0) {
    hList->first = initialize_halo_node(err);      forwardError(*err, __LINE__,);
    hList->last  = hList->first;
    hList->length++;
  }
  else if (hList->length == hList->size) {
    hList->last->next = initialize_halo_node(err); forwardError(*err, __LINE__,);
    hList->last       = hList->last->next;
    hList->length++;
  }
  else if (hList->size == 0) hList->last = hList->first;
  else                       hList->last = hList->last->next;
  
  h = hList->last->h;
  set_halo_t(cmhm, h, z, M, pos, err);             forwardError(*err, __LINE__,);
  hList->size++;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to halo_map

halo_map *initialize_halo_map(int N1, int N2, double theta_pix, error **err)
{
  halo_map *hMap  = (halo_map*)malloc_err(sizeof(halo_map), err);                    forwardError(*err, __LINE__,);
  hMap->N1        = N1;
  hMap->N2        = N2;
  hMap->length    = N1 * N2;
  hMap->total     = 0;
  hMap->theta_pix = theta_pix;
  
  hMap->center[0] = 0.0;
  hMap->center[1] = 0.0;
  hMap->limits[0] = -0.5 * N1 * theta_pix;
  hMap->limits[1] = -hMap->limits[0];
  hMap->limits[2] = -0.5 * N2 * theta_pix;
  hMap->limits[3] = -hMap->limits[2];
  
  hMap->map       = (halo_list**)malloc_err(hMap->length * sizeof(halo_list*), err); forwardError(*err, __LINE__,);
  int i;
  for (i=0; i<hMap->length; i++) hMap->map[i] = initialize_halo_list(err);
  return hMap;
}

void free_halo_map(halo_map *hMap)
{
  int i;
  for (i=0; i<hMap->length; i++) {free_halo_list(hMap->map[i]); hMap->map[i] = NULL;}
  free(hMap); hMap = NULL;
  return;
}

void reset_halo_map(halo_map *hMap)
{
  hMap->total = 0;
  int i;
  for (i=0; i<hMap->length; i++) reset_halo_list(hMap->map[i]);
  return;
}

void append_halo_map(cosmo_hm *cmhm, halo_map *hMap, double z, double M, double pos[2], error **err)
{
  double theta_pix = hMap->theta_pix;
  int i = (int)((pos[0] - hMap->limits[0]) / theta_pix);
  int j = (int)((pos[1] - hMap->limits[2]) / theta_pix);
  append_halo_list(cmhm, hMap->map[i+j*hMap->N1], z, M, pos, err); forwardError(*err, __LINE__,);
  hMap->total++;
  return;
}

void read_halo_map(char name[], cosmo_hm *cmhm, halo_map *hMap, error **err)
{
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  printf("Reading...\r");
  fflush(stdout);
  
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double pos[2], w, z, M;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf %lf %*f %lf %lf\n", &pos[0], &pos[1], &z, &M);
      
      append_halo_map(cmhm, hMap, z, M, pos, err); forwardError(*err, __LINE__,);
      count++;
    }
    c = fgetc(file);
  }
  
  fclose(file);
  testError(count!=hMap->total, peak_match, "Halo number match error", *err, __LINE__); forwardError(*err, __LINE__,);
  printf("\"%s\" read       \n", name);
  printf("%d halos generated\n", count);
  return;
}

void output_halo_map(FILE *file, peak_param *peak, halo_map *hMap)
{
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  
  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file, "#       RA        DEC         w        z          M\n");
    fprintf(file, "#     [deg]      [deg]   [Mpc/h]      [-]  [M_sol/h]\n");
  }
  else {
    fprintf(file, "#   theta_x    theta_y        w        z          M\n");
    fprintf(file, "#  [arcmin]   [arcmin]   [Mpc/h]      [-]  [M_sol/h]\n");
  }
  
  halo_list *hList;
  halo_node *hNode;
  halo_t *h;
  int i, j;
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      h = hNode->h;
      fprintf(file, "  %9.3f  %9.3f  %8.3f  %7.5f  %9.3e \n", h->pos[0], h->pos[1], h->w, h->z, h->M);
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mass function

double massFct(cosmo_hm_params *cANDp, double mass, error **err)
{
  double dn = LN_10 * dn_dlnM(mass, (void*)cANDp, err); forwardError(*err, __LINE__,);
  return dn;
}

void fillMassFct(cosmo_hm_params *cANDp, sampler_t *samp, error **err)
{
  double *x    = samp->x;
  double *pdf  = samp->pdf;
  int i;
  
  //-- Fill pdf
  for (i=0; i<samp->length; i++) {
    pdf[i] = massFct(cANDp, x[i], err);
    forwardError(*err, __LINE__,);
  }
  
  //-- Precompute cdf, total pdf, <x>
  int setTotalToOne = 0; //-- This is to keep the total of pdf.
  set_sampler_t(samp, setTotalToOne);
  return;
}

void outputMassFct(char name[], cosmo_hm *cmhm, peak_param *peak, double z, error **err)
{
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = cmhm;
  cANDp->asymptotic = 0;
  cANDp->a = 1.0/(1.0+z);
  
  FILE *file;
  if (name != NULL) file = fopen(name, "w");
  else              file = stdout;
  
  fprintf(file, "# Halo mass function\n");
  fprintf(file, "# Model = %s, redshift = %g\n", smassfct_t(cmhm->massfct), z);
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  fprintf(file, "#          M      dn/dlogM\n");
  fprintf(file, "#   [M_sol/h]  [(Mpc/h)^-3]\n");
  
  double M, dn;
  for (M=peak->M_min; M<=peak->M_max; M*=1.01) {
    dn = massFct(cANDp, M, err); forwardError(*err, __LINE__,);
    fprintf(file, " %12.6e  %12.6e\n", M, dn);
  }
  
  fclose(file);
  free(cANDp);
  if (name != NULL) printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to fast simulations

double dVol(cosmo_hm *cmhm, peak_param *peak, double z, double ww, double dz, error **err)
{
  //-- Comoving volume between z and z+dz, peak->area in [arcmin^2]
  double dV, A;
  double fK = f_K(cmhm->cosmo, ww, err);                         forwardError(*err, __LINE__, 0);
  dV  = HUBBLE_DISTANCE * SQ(fK);
  dV *= peak->area * SQ(ARCMIN_TO_RADIAN);
  int wOmegar = 0;
  dV *= dz / sqrt(Esqr(cmhm->cosmo, 1.0/(1.0+z), wOmegar, err)); forwardError(*err, __LINE__, 0);
  return dV;
}

void randomizePosition(peak_param *peak, double *pos, error **err)
{
  double *Omega      = peak->Omega;
  gsl_rng *generator = peak->generator;
  
  //-- Set random position n in solid angle Omega
  double theta_x, theta_y, theta_sq;
  switch (peak->field) {
    case circle:
      theta_sq = SQ(Omega[0]);
      do {
	theta_x = gsl_ran_flat(generator, -1.0, 1.0) * Omega[0];
	theta_y = gsl_ran_flat(generator, -1.0, 1.0) * Omega[0];
      } while (SUM_SQ_2(theta_x, theta_y) > theta_sq);
      break;
      
    case rectangle:
      theta_x = gsl_ran_flat(generator, -0.5, 0.5) * Omega[0];
      theta_y = gsl_ran_flat(generator, -0.5, 0.5) * Omega[1];
      break;
      
    default :
      peak->randPosFct(generator, pos);
      return;
   }

   pos[0] = theta_x;
   pos[1] = theta_y;
   return;
}

void setMassSamplers(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, error **err)
{
  if (peak->printMode == 0) {
    printf("Computing mass function...\r");
    fflush(stdout);
  }
  
  //-- Redshift
  double z_max;
  if (peak->z_halo_max < 0) z_max = get_zmax(cmhm->redshift, 0);
  else                      z_max = peak->z_halo_max;
  double dz = z_max / (double)(peak->N_z_halo); 
  
  //-- Structure for mass function
  cosmo_hm_params *cANDp = (cosmo_hm_params*)malloc_err(sizeof(cosmo_hm_params), err); forwardError(*err, __LINE__,);
  cANDp->model = cmhm;
  cANDp->asymptotic = 0;
  
  sampler_t *samp;
  double ww, dV;
  double r = pow(10, peak->dlogM); //-- Convert bin width to ratio
  int wOmegar = 0;
  
  //-- Loop over redshifts
  double z;
  int i, j;
  for (i=0, z=0.5*dz; i<sampArr->length; i++, z+=dz) {
    samp       = sampArr->array[i];
    samp->dx   = peak->dlogM;
    samp->x[0] = peak->M_min;
    for (j=1; j<peak->nbMassBins; j++) samp->x[j] = samp->x[j-1] * r; //-- Fill masses
    
    //-- Compute the volume
    cANDp->a = 1.0/(1.0+z);
    ww = w(cmhm->cosmo, cANDp->a, wOmegar, err);  forwardError(*err, __LINE__,);
    dV = dVol(cmhm, peak, z, ww, dz, err);        forwardError(*err, __LINE__,);
    
    //-- Fill mass function and cdf
    fillMassFct(cANDp, samp, err);                forwardError(*err, __LINE__,);
    samp->x_mean *= dV; //-- x_mean becomes total mass in dV
  }
  
  if (peak->printMode != 1) printf("Mass function computation done\n");
  return;
}

void sampleHalos(cosmo_hm *cmhm, peak_param *peak, sampler_t *samp, halo_map *hMap, double z, error **err)
{
  gsl_rng *generator = peak->generator;
  double totalMass   = samp->x_mean;
  double sum = 0.0;
  double p, M, pos[2];
  
  //-- Sample masses
  while (sum < totalMass) {
    p = gsl_ran_flat(generator, 0.0, 1.0);
    M = execute_sampler_t(samp, p);
    
    randomizePosition(peak, pos, err);           forwardError(*err, __LINE__,);
    append_halo_map(cmhm, hMap, z, M, pos, err); forwardError(*err, __LINE__,);
    sum += M;
  }
  return;
}

void makeFastSimul(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, error **err)
{
  //-- Reset
  reset_halo_map(hMap);
  
  //-- Redshift
  double z_max;
  if (peak->z_halo_max < 0) z_max = get_zmax(cmhm->redshift, 0);
  else                      z_max = peak->z_halo_max;
  double dz = z_max / (double)(peak->N_z_halo);
      
  //-- Loop over redshifts
  double z;
  int i;
  for (i=0, z=0.5*dz; i<peak->N_z_halo; i++, z+=dz) {
    sampleHalos(cmhm, peak, sampArr->array[i], hMap, z, err);
    forwardError(*err, __LINE__,);
  }
  
  if (peak->printMode == 1) printf("%6d halos, ", hMap->total);
  else                      printf("%d halos generated         \n", hMap->total);
  return;
}

void outputFastSimul(char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  
  output_halo_map(file, peak, hMap);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doMassFct(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  char name[STRING_LENGTH_MAX];
  double z;
  for (z=0.05; z<1.0; z+=0.1) {
    sprintf(name, "massFct_%s_z%.3f", smassfct_t(cmhm->massfct), z);
    outputMassFct(name, cmhm, peak, z, err);
  }
  printf("------------------------------------------------------------------------\n");
  return;
}

void doFastSimulation(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
  setMassSamplers(cmhm, peak, sampArr, err);                                                        forwardError(*err, __LINE__,);
  halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  
  makeFastSimul(cmhm, peak, sampArr, hMap, err); forwardError(*err, __LINE__,);
  outputFastSimul("haloList", cmhm, peak, hMap);
  
  free_sampler_arr(sampArr);
  free_halo_map(hMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------


