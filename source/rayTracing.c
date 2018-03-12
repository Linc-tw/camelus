

  /***************************************************************
   **  rayTracing.c						**
   **  Chieh-An Lin, Martin Kilbinger, François Lanusse		**
   **  Version 2016.03.22					**
   **								**
   **  References:						**
   **  - Baltz et al. (2009) - JCAP, 1, 15			**
   **  - Bartelmann & Schneider (2001) - Phys. Rep., 340, 291	**
   **  - Oguri & Hamana (2011) - MNRAS, 414, 1851		**
   **  - Takada & Jain (2003a) - MNRAS, 340, 580		**
   **  - Takada & Jain (2003b) - MNRAS, 344, 857		**
   **  - Wright & Brainerd (2000) - ApJ, 534, 34		**
   ***************************************************************/


#include "rayTracing.h"


//----------------------------------------------------------------------
//-- Functions related to gal_t, gal_node, gal_list

void set_gal_t(cosmo_hm *cmhm, gal_t *g, double z, double w_s, double D_s, double pos[2], error **err)
{
  testErrorRet(g==NULL, peak_null, "gal_t *g is NULL", *err, __LINE__,);
  
  g->z        = z;
  double a    = 1.0/(1.0+z);
  g->a        = a;
  g->kappa    = 0.0;
  g->gamma[0] = 0.0;
  g->gamma[1] = 0.0;
  
  if (w_s < 0) {
    //-- Use z to compute w_s and D_s
    g->w   = w(cmhm->cosmo, a, 0, err);       forwardError(*err, __LINE__,); //-- wOmegar = 0
    g->D_s = a * f_K(cmhm->cosmo, g->w, err); forwardError(*err, __LINE__,);
  }
  else {
    g->w   = w_s;
    g->D_s = D_s;
  }     
  
  if (pos != NULL) {
    g->pos[0] = pos[0];
    g->pos[1] = pos[1];
  }
  return;
}

void setWithSignal_gal_t(cosmo_hm *cmhm, gal_t *g, double z, double pos[2], double kappa, double gamma[2], error **err)
{
  //-- WARNING: No w and D_s, only used in read_gal_map and fast pipeline
  testErrorRet(g==NULL, peak_null, "gal_t *g is NULL", *err, __LINE__,);
  
  g->z        = z;
  g->a        = 1.0/(1.0+z);
  g->kappa    = kappa;
  g->gamma[0] = gamma[0];
  g->gamma[1] = gamma[1];
  g->pos[0]   = pos[0];
  g->pos[1]   = pos[1];
  g->w;
  g->D_s;
  return;
}

gal_node *initialize_gal_node(error **err)
{
  gal_node *gNode = (gal_node*)malloc_err(sizeof(gal_node), err); forwardError(*err, __LINE__, NULL);
  gNode->g        = (gal_t*)malloc_err(sizeof(gal_t), err);       forwardError(*err, __LINE__, NULL);;
  gNode->next     = NULL;
  return gNode;
}

gal_list *initialize_gal_list(error **err)
{
  gal_list *gList = (gal_list*)malloc_err(sizeof(gal_list), err); forwardError(*err, __LINE__, NULL);
  gList->length   = 0;
  gList->size     = 0;
  gList->mean[0]  = 0.0;
  gList->mean[1]  = 0.0;
  gList->first    = NULL;
  gList->last     = NULL;
  return gList;
}

void free_gal_list(gal_list *gList)
{
  gal_node *gNode;
  while (gList->first != NULL) {
    gNode        = gList->first;
    gList->first = gNode->next;
    if (gNode->g) {free(gNode->g); gNode->g = NULL;}
    free(gNode); gNode = NULL;
  }
  free(gList); gList = NULL;
  return;
}

void cleanLensing_gal_list(gal_list *gList)
{
  gal_node *gNode = gList->first;
  while (gNode != NULL) {
    gNode->g->kappa    = 0.0;
    gNode->g->gamma[0] = 0.0;
    gNode->g->gamma[1] = 0.0;
    gNode              = gNode->next;
  }
  return;
}

void reset_gal_list(gal_list *gList)
{
  gList->size    = 0;
  gList->mean[0] = 0.0;
  gList->mean[1] = 0.0;
  return;
}

void append_gal_list(cosmo_hm *cmhm, gal_list *gList, double z, double w_s, double D_s, double pos[2], error **err)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node(err);              forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node(err);         forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  set_gal_t(cmhm, gList->last->g, z, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
  gList->size++;
  return;
}

void appendWithSignal_gal_list(cosmo_hm *cmhm, gal_list *gList, double z, double pos[2], double kappa, double gamma[2], error **err)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node(err);                            forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node(err);                       forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  setWithSignal_gal_t(cmhm, gList->last->g, z, pos, kappa, gamma, err); forwardError(*err, __LINE__,);
  gList->size++;
  return;
}

void kappaMean_gal_list(gal_list *gList, double threshold)
{
  int size     = gList->size;
  double *mean = gList->mean;
  if (size == 0 || (double)size < threshold) {
    mean[0] = 0.0;
    mean[1] = 0.0;
    return;
  }
  
  double mean1 = 0.0;
  gal_node *gNode;
  gal_t *g;
  
  int i;
  for (i=0, gNode=gList->first; i<size; i++, gNode=gNode->next) mean1 += gNode->g->kappa;
  mean[0] = mean1 / (double)size;
  mean[1] = 0.0;
  return;
}

void gammaOrGMean_gal_list(gal_list *gList, double threshold)
{
  int size     = gList->size;
  double *mean = gList->mean;
  if (size == 0 || (double)size < threshold) {
    mean[0] = 0.0;
    mean[1] = 0.0;
    return;
  }
  
  double mean1  = 0.0;
  double mean2  = 0.0;
  double factor = 1.0 / (double)size;
  gal_node *gNode;
  gal_t *g;
  
  int i;
  for (i=0, gNode=gList->first; i<size; i++, gNode=gNode->next) {
    mean1 += gNode->g->gamma[0];
    mean2 += gNode->g->gamma[1];
  }
  mean[0] = mean1 * factor;
  mean[1] = mean2 * factor;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to gal_map

gal_map *initialize_gal_map(int N1, int N2, double theta_pix, error **err)
{
  gal_map *gMap          = (gal_map*)malloc_err(sizeof(gal_map), err);                    forwardError(*err, __LINE__, NULL);
  gMap->N1               = N1;
  gMap->N2               = N2;
  gMap->length           = N1 * N2;
  gMap->total            = 0;
  gMap->theta_pix        = theta_pix;
  gMap->theta_pix_inv    = 1.0 / theta_pix;
  
  gMap->limits[0]        = 0;
  gMap->limits[1]        = N1 * theta_pix;
  gMap->limits[2]        = 0;
  gMap->limits[3]        = N2 * theta_pix;
  gMap->center[0]        = 0.5 * gMap->limits[1];
  gMap->center[1]        = 0.5 * gMap->limits[3];
  
  gMap->map              = (gal_list**)malloc_err(gMap->length * sizeof(gal_list*), err); forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<gMap->length; i++) gMap->map[i] = initialize_gal_list(err);
  
  gMap->kappa_mean       = 0.0;
  gMap->fillingThreshold = 0.0;
  gMap->type             = kappa_map;
  return gMap;
}

void free_gal_map(gal_map *gMap)
{
  int i;
  if (gMap->map) {
    for (i=0; i<gMap->length; i++) {free_gal_list(gMap->map[i]); gMap->map[i] = NULL;}
    free(gMap->map); gMap->map = NULL;
  }
  free(gMap); gMap = NULL;
  return;
}

void cleanLensing_gal_map(gal_map *gMap)
{
  gMap->kappa_mean = 0.0;
  int i;
  for (i=0; i<gMap->length; i++) cleanLensing_gal_list(gMap->map[i]);
  return;
}

void reset_gal_map(gal_map *gMap)
{
  gMap->total = 0;
  gMap->kappa_mean = 0.0;
  gMap->fillingThreshold = 0.0;
  int i;
  for (i=0; i<gMap->length; i++) reset_gal_list(gMap->map[i]);
  return;
}

void append_gal_map(cosmo_hm *cmhm, gal_map *gMap, double z, double w_s, double D_s, double pos[2], error **err)
{
  double theta_pix_inv = gMap->theta_pix_inv;
  int i = (int)(pos[0] * theta_pix_inv);
  int j = (int)(pos[1] * theta_pix_inv);
  append_gal_list(cmhm, gMap->map[i+j*gMap->N1], z, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
  gMap->total++;
  return;
}

void appendWithSignal_gal_map(cosmo_hm *cmhm, gal_map *gMap, double z, double pos[2], double kappa, double gamma[2], error **err)
{
  double theta_pix_inv = gMap->theta_pix_inv;
  int i = (int)(pos[0] * theta_pix_inv);
  int j = (int)(pos[1] * theta_pix_inv);
  if (i >= gMap->N1) i -= 1;
  if (j >= gMap->N2) j -= 1;
  appendWithSignal_gal_list(cmhm, gMap->map[i+j*gMap->N1], z, pos, kappa, gamma, err); forwardError(*err, __LINE__,);
  gMap->total++;
  return;
}

void setFillingThreshold(peak_param *peak, gal_map *gMap, error **err)
{
  if (gMap->fillingThreshold != 0.0) return; //-- Already computed
  if (peak->doMask == 0) {
    gMap->fillingThreshold = -1.0;
    if      (peak->printMode < 2)   printf("No mask\n");
    else if (peak->printMode == 2) {printf("no mask, "); fflush(stdout);}
    return;
  }
  
  int bufferSize = peak->bufferSize;
  testErrorRet(bufferSize<=0, peak_badValue, "Buffer size should be at least 1.", *err, __LINE__,);
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  double threshold = 0.0;
  int i, j, jN1;
  
  //-- Set filling threshold
  for (j=bufferSize; j<N2-bufferSize; j++) { 
    jN1 = j * N1;
    for (i=bufferSize; i<N1-bufferSize; i++) threshold += gMap->map[i+jN1]->size;
  }
  threshold /= (double)((N1 - 2*bufferSize) * (N2 - 2*bufferSize)); //-- Average filling factor
  threshold *= FILLING_THRESHOLD_RATIO;                             //-- Set threshold to the half of average
  gMap->fillingThreshold = threshold;
  
  if      (peak->printMode < 2)   printf("Filling threshold = %.2f\n", threshold);
  else if (peak->printMode == 2) {printf("threshold = %.2f, ", threshold); fflush(stdout);}
  return;
}

void signalMean_gal_map(peak_param *peak, gal_map *gMap, error **err)
{
  setFillingThreshold(peak, gMap, err); forwardError(*err, __LINE__,);
  int i;
  if (peak->doKappa == 1) {
    for (i=0; i<gMap->length; i++) kappaMean_gal_list(gMap->map[i], gMap->fillingThreshold);
    gMap->type = kappa_map;
  }
  else {
    for (i=0; i<gMap->length; i++) gammaOrGMean_gal_list(gMap->map[i], gMap->fillingThreshold);
    if (peak->doKappa >= 2) gMap->type = g_map;
    else                    gMap->type = gamma_map;
  }
  return;
}

void output_gal_map(FILE *file, peak_param *peak, gal_map *gMap)
{
  fprintf(file, "# Type = %s, number of galaxies = %d\n", STR_MAPTYPE_T(gMap->type), gMap->total);
  fprintf(file, "#\n");
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int i, j;
  
  if (gMap->type < 12) {
    fprintf(file, "#  theta_x    theta_y       z      kappa    gamma1    gamma2\n");
    fprintf(file, "# [arcmin]   [arcmin]      [-]       [-]       [-]       [-]\n");
    
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	g = gNode->g;
	fprintf(file, " %9.3f  %9.3f  %7.5f  %8.5f  %8.5f  %8.5f\n", g->pos[0], g->pos[1], g->z, g->kappa, g->gamma[0], g->gamma[1]);
      }
    }
  }
  else {
    fprintf(file, "#  theta_x    theta_y       z      kappa       g1        g2\n");
    fprintf(file, "# [arcmin]   [arcmin]      [-]       [-]       [-]       [-]\n");
    
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	g = gNode->g;
	fprintf(file, " %9.3f  %9.3f  %7.5f  %8.5f  %8.5f  %8.5f\n", g->pos[0], g->pos[1], g->z, g->kappa, g->gamma[0], g->gamma[1]);
      }
    }
  }
  return;
}

void read_gal_map(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  //-- WARNING: D_S, w_s are not read.
  
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  printf("Reading...\r");
  fflush(stdout);
  
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double pos[2], gamma[2], z, kappa;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf %lf %lf %lf %lf %lf\n", &pos[0], &pos[1], &z, &kappa, &gamma[0], &gamma[1]);
      
      appendWithSignal_gal_map(cmhm, gMap, z, pos, kappa, gamma, err); forwardError(*err, __LINE__,);
      count++;
    }
    c = fgetc(file);
  }
  
  fclose(file);
  testErrorRet(count!=gMap->total, peak_match, "Galaxy number match error", *err, __LINE__,);
  printf("\"%s\" read       \n", name);
  printf("%d galaxies generated\n", count);
  return;
}


void updateCosmo_gal_map(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  double w_s, D_s;
  if (peak->z_s > 0) {
    w_s  = w(cmhm->cosmo, 1.0/(1.0 + peak->z_s), 0, err); forwardError(*err, __LINE__,); //-- wOmegar = 0
    D_s  = f_K(cmhm->cosmo, w_s, err);                    forwardError(*err, __LINE__,);
    D_s /= 1.0 + peak->z_s;
  }
  else {
    w_s = -1.0; //-- This means that w_s needs to be calculated in set_gal_t
    D_s = -1.0; //-- This means that D_s needs to be calculated in set_gal_t
  }
  
  //-- Loop for setting galaxies
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int i, j;
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      g = gNode->g;
      set_gal_t(cmhm, g, g->z, w_s, D_s, NULL, err); forwardError(*err, __LINE__,);
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to projected mass

#define EPS 1e-10
double G_NFW_kappa(double u_sq, double c, double c_sq, error **err)
{
  //-- See Wright & Brainerd (2000), Eq. (11)
  //-- Sigma = 2 rho_s r_s G(u)
  //-- G(u) = -1/arg1 + arcosh(1/u) / arg1^(3/2)  if u < 1,
  //--         1/3                                if u = 1,
  //--         1/arg1 - arccos(1/u) / arg1^(3/2)  if u > 1,
  //-- where
  //--   arg1 = |1 - u^2|.
  //-- Note that
  //--   arcosh(1 / u) = 2 artanh(p),
  //--   arccos(1 / u) = 2 arctan(p),
  //-- where
  //--   p = sqrt(|1 - u| / (1 + u)) = sqrt(arg1) / (1 + u).
  
  if (u_sq < EPS) u_sq = EPS;
  
  double arg1 = fabs(u_sq - 1.0);
  double value;
  
  if (arg1 < EPS) value = 1.0 / 3.0;
  else {
    if (u_sq < 1.0) {
      value  = -1.0 + acosh(1.0 / sqrt(u_sq)) / sqrt(arg1);
      value /= arg1;
    }
    else {
      value  = 1.0 - acos(1.0 / sqrt(u_sq)) / sqrt(arg1);
      value /= arg1;
    }
  }
  
  return value;
}

double G_NFW_gamma(double u_sq, double c, double c_sq, double f, error **err)
{
  //-- See Wright & Brainerd (2000), Eqs. (15) and (16)
  //-- G(u) = 1/arg1 * (1 - arg3) + 2 u^-2 (arg2 + arg3)    if u < 1,
  //--        5/3 + 2 ln(1/2)                               if u = 1,
  //--        1/arg1 * (-1 + arg4) + 2 u^-2 (arg2 + arg4)]  if u > 1,
  //-- where
  //--   arg1 = |1 - u^2|,
  //--   arg2 = ln(u / 2),
  //--   arg3 = arcosh(1/u) / q,
  //--   arg4 = arccos(1/u) / q,
  //-- and
  //--   q = sqrt(|1 - u^2|) = sqrt(arg1).
  //-- Note that
  //--   arcosh(1 / u) = 2 artanh(p),
  //--   arccos(1 / u) = 2 arctan(p),
  //-- where
  //--   p = sqrt(|1 - u| / (1 + u)) = sqrt(arg1) / (1 + u).
  
  if (u_sq < EPS) u_sq = EPS;
  
  double arg1 = fabs(u_sq - 1.0);
  double value;
  
  if (arg1 < EPS) {
    value = 5.0 / 3.0 + 2.0 * log(0.5);
  }
  else {
    double u    = sqrt(u_sq);
    double q    = sqrt(arg1);
    double arg2 = log(0.5 * u);
    if (u_sq < 1.0) {
      double arg3 = acosh(1.0 / u) / q;
      value = (1.0 - arg3) / arg1 + 2.0 / u_sq * (arg2 + arg3);
    }
    else {
      double arg4 = acos(1.0 / u) / q;
      value = (-1.0 + arg4) / arg1 + 2.0 / u_sq * (arg2 + arg4);
    }
  }
  
  return value;
}

void G_NFW_both(double u_sq, double c, double c_sq, double f, double G_NFW[2], error **err)
{
  if (u_sq < EPS) u_sq = EPS;
  
  double arg1 = fabs(u_sq - 1.0);
  double value;
  
  if (arg1 < EPS) {
    G_NFW[0] = 1.0 / 3.0;
    G_NFW[1] = 5.0 / 3.0 + 2.0 * log(0.5);
  }
  else {
    double u    = sqrt(u_sq);
    double q    = sqrt(arg1);
    double arg2 = log(0.5 * u);
    if (u_sq < 1.0) {
      double arg3 = acosh(1.0 / u) / q;
      G_NFW[0] = (-1.0 + arg3) / arg1;
      G_NFW[1] = -G_NFW[0] + 2.0 / u_sq * (arg2 + arg3);
    }
    else {
      double arg4 = acos(1.0 / u) / q;
      G_NFW[0] = (1.0 - arg4) / arg1;
      G_NFW[1] = -G_NFW[0] + 2.0 / u_sq * (arg2 + arg4);
    }
  }
  
  return;
}

double G_TJ_kappa(double u_sq, double c, double c_sq, error **err)
{
  //-- See Takada & Jain (2003a), Eq. (27)
  //-- G(u) = -arg2/arg1 + arcosh(arg3) / arg1^(3/2)  if u < 1,
  //--         arg2/3 * (1 + 1/p)                     if u = 1,
  //--         arg2/arg1 - arccos(arg3) / arg1^(3/2)  if 1 < u <= c,
  //--        0                                       if u > c,
  //-- where
  //--   arg1 = |q|,
  //--   arg2 = r / p,
  //--   arg3 = (u^2 + c) / u p,
  //-- and
  //--   p = c + 1,
  //--   q = 1 - u^2,
  //--   r = sqrt(c^2 - u^2).
  
  
  if (u_sq > c_sq - EPS) return 0;
  if (u_sq < EPS) u_sq = EPS;
  
  double p_inv = 1.0 / (c + 1.0);
  double arg1  = fabs(u_sq - 1.0);
  double arg2  = sqrt(c_sq - u_sq) * p_inv;
  double value;
  
  if (arg1 < EPS) value = arg2 * (1.0 + p_inv) / 3.0;
  else {
    double arg3 = (u_sq + c) * p_inv / sqrt(u_sq);
    if (u_sq < 1.0) {
      testErrorRet(arg3<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__, -1.0);
      value  = -arg2 + acosh(arg3) / sqrt(arg1);
      value /= arg1;
    }
    else {
      testErrorRet(arg3>1.0, peak_badValue, "Out of range for acos", *err, __LINE__, -1.0);
      value  = arg2 - acos(arg3) / sqrt(arg1);
      value /= arg1;
    }
  }
  
  return value;
}

double G_TJ_gamma(double u_sq, double c, double c_sq, double f, error **err)
{
  //-- See Takada & Jain (2003b), Eq. (17)
  //-- G(x) = u^-2 [arg4 + 2 ln(arg5) + arg6 * arcosh(arg3)]  if u < 1,
  //--        u^-2 [arg7 + 2 ln(arg5)]                        if u = 1,
  //--        u^-2 [arg4 + 2 ln(arg5) + arg6 * arccos(arg3)]  if 1 < u <= c,
  //--        u^-2 * 2/f                                      if u > c,
  //-- where
  //--   arg3 = (u^2 + c) / u p,
  //--   arg4 = [(q + 1) r / q - 2c] / p,
  //--   arg5 = u p / (c + r),
  //--   arg6 = (3q - 1) / q sqrt|q|,
  //--   arg7 = [(c + 10p) r / p - 6c] / 3p,
  //-- and
  //--   p = c + 1,
  //--   q = 1 - u^2,
  //--   r = sqrt(c^2 - u^2).
  
  if (u_sq > c_sq - EPS) return 2.0 / (f * u_sq);
  if (u_sq < EPS) u_sq = EPS;
  
  double p     = c + 1.0;
  double p_inv = 1.0 / p;
  double q     = 1.0 - u_sq;
  double r     = sqrt(c_sq - u_sq);
  double u     = sqrt(u_sq);
  double arg5  = u * p / (c + r);
  double value;
  
  if (fabs(q) < EPS) {
    value  = ((10.0 + c * p_inv) * r - 6.0*c) * p_inv / 3.0 + 2.0 * log(arg5);
    value /= u_sq;
    return value;
  }
  
  double q_inv = 1.0 / q;
  double arg3  = (u_sq + c) * p_inv / u;
  double arg4  = ((1.0 + q_inv) * r - 2.0 * c) * p_inv;
  double arg6  = (3.0 - q_inv) / sqrt(fabs(q));
  
  if (q > 0) {
    testErrorRet(arg3<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__, -1.0);
    value  = arg4 + 2.0 * log(arg5) + arg6 * acosh(arg3);
    value /= u_sq;
  }
  else {
    testErrorRet(arg3>1.0, peak_badValue, "Out of range for acos", *err, __LINE__, -1.0);
    value  = arg4 + 2.0 * log(arg5) + arg6 * acos(arg3);
    value /= u_sq;
  }
  
  return value;
}

void G_TJ_both(double u_sq, double c, double c_sq, double f, double G_TJ[2], error **err)
{
  if (u_sq > c_sq - EPS) {
    G_TJ[0] = 0.0;
    G_TJ[1] = 2.0 / (f * u_sq);
    return;
  }
  if (u_sq < EPS) u_sq = EPS;
  
  double p     = c + 1.0;
  double p_inv = 1.0 / p;
  double q     = 1.0 - u_sq;
  double r     = sqrt(c_sq - u_sq);
  double u     = sqrt(u_sq);
  double arg1  = fabs(q);
  double arg2  = r * p_inv;
  double arg5  = u * p / (c + r);
  double value;
  
  if (arg1 < EPS) {
    value     = p_inv / 3.0;
    G_TJ[0]  = r * (1.0 + p_inv) * value;
    G_TJ[1]  = ((10.0 + c * p_inv) * r - 6.0*c) * value + 2.0 * log(arg5);
    G_TJ[1] /= u_sq;
    return;
  }
  
  double q_inv = 1.0 / q;
  double arg3  = (u_sq + c) * p_inv / u;
  double arg4  = ((1.0 + q_inv) * r - 2.0 * c) * p_inv;
  
  if (q > 0) {
    testErrorRet(arg3<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__,);
    value    = acosh(arg3) / sqrt(arg1);
    G_TJ[0]  = -arg2 + value;
    G_TJ[0] /= arg1;
    G_TJ[1]  = arg4 + 2.0 * log(arg5) + (3.0 - q_inv) * value;
    G_TJ[1] /= u_sq;
  }
  else {
    testErrorRet(arg3>1.0, peak_badValue, "Out of range for acos", *err, __LINE__,);
    value    = acos(arg3) / sqrt(arg1);
    G_TJ[0]  = arg2 - value;
    G_TJ[0] /= arg1;
    G_TJ[1]  = arg4 + 2.0 * log(arg5) + (3.0 - q_inv) * value;
    G_TJ[1] /= u_sq;
  }
  
  return;
}

double G_BMO_kappa(double u_sq, double tau, double tau_sq, error **err)
{
  //-- See Baltz et al. (2009), Eq. (A.28)
  //-- G(x) = factor * [arg5 (-1 + arg6) + 8 arg6 + arg1 - arg2 + (arg3 + arg4) L]  if u < 1,
  //--        factor * [2p / 3           + 8      + arg1 - arg2 + (arg3 + arg4) L]  if u = 1,
  //--        factor * [arg5 (1 - arg7)  + 8 arg7 + arg1 - arg2 + (arg3 + arg4) L]  if u > 1,
  //-- where
  //--   arg1 = q / tau^2 r^2,
  //--   arg2 = pi (4 r^2 + p) / r^3,
  //--   arg3 = q / tau r^3,
  //--   arg4 = (3q - 6p + 8) / tau^3 r,
  //--   arg5 = 2 p / s,
  //--   arg6 = arcosh(1/u) / sqrt(s),
  //--   arg7 = arccos(1/u) / sqrt(s),
  //-- and
  //--   p = tau^2 + 1,
  //--   q = tau^4 - 1,
  //--   r = sqrt(tau^2 + u^2),
  //--   s = |1 - u^2|,
  //--   L = ln(u / (r + tau)),
  //--   factor = (q + 1) / 2 p^3.
  
  if (u_sq < EPS) u_sq = EPS;
  
  double p = tau_sq + 1.0;
  double q = tau_sq * tau_sq - 1.0;
  double r_sq = tau_sq + u_sq;
  double r = sqrt(r_sq);
  double s = r * tau;
  double t = fabs(1.0 - u_sq);
  double u = sqrt(u_sq);
  
  double L    = log(u / (r + tau));
  double arg1 = q / (s * s);
  double arg2 = PI * (4.0*r_sq + p) / (r_sq * r);
  double arg3 = q / (r_sq * s);
  double arg4 = (3.0*q - 6.0*p + 8.0) / (tau_sq * s);
  
  double value = arg1 - arg2 + (arg3 + arg4) * L;
  
  if (t < EPS) value += 2.0*p / 3.0 + 8.0;
  else {
    double arg5 = 2.0 * p / t;
    if (u_sq < 1.0) {
      double arg6 = acosh(1.0 / u) / sqrt(t);
      value += arg5 * (-1.0 + arg6) + 8.0 * arg6;
    }
    else {
      double arg7 = acos(1.0 / u) / sqrt(t);
      value += arg5 * (1.0 - arg7) + 8.0 * arg7;
    }
  }
  
  double factor = 0.5 * (q + 1) / (p * (2.0*p + q));
  value *= factor;
  return value;
}

#undef EPS

double kappa_TJ(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- See Bartelmann & Schneider (2001), Eq. (3.7)
  //-- kappa_TJ = Sigma / Sigma_crit
  //-- where
  //--   Sigma_crit = (c^2 / 4 pi G) * (D_s / D_l D_ls), where (4 pi G / c^2) = 3 / (2 CRITICAL_DENSITY D_H^2 ) = 6.01e-19 [Mpc/M_sol]
  //--   Sigma = projected mass
  //--         = 2 rho_s r_s * G_TJ_kappa(theta / theta_s)
  //--         = 2 * (rho_bar * Delta f_NFW c_NFW^3 / 3) * (r_vir / c_NFW) * G_TJ_kappa(c_NFW theta / theta_vir)
  //--         = 2 * (rho_bar * Delta * 4 pi r_vir^3 / 3) * (f_NFW c_NFW^3 / 4 pi r_vir^3) * (r_vir / c_NFW) * G_TJ_kappa(c_NFW theta / theta_vir)
  //--         = M * (f_NFW c_NFW^2 / 2 pi r_vir^2) * G_TJ_kappa(c_NFW theta / theta_vir)
  //-- 
  //-- Similarly, in Takada & Jain (2003a), Eqs. (26) and (28), one finds
  //-- kappa_TJ = W(w_l, w_s) * (M / rho_bar) * Sigma_m(theta)
  //--          = (3 Omega_m H_0^2 / 2) * (D_l D_ls / D_s a_l) * (M / rho_bar) * Sigma_m
  //--          = (4 pi G) * (D_l D_ls / D_s a_l) * M * (f_NFW c_NFW^2 / 2 pi R_vir^2) * G_TJ_kappa(c_NFW theta / theta_vir)
  //-- However, their convention is different. First, they have omitted c^(-2).
  //-- Second, their D_l, D_ls, D_s are actually f_K: f_l, f_ls, f_s.
  //-- Finally, their R_vir is defined as a comoving distance, not the physical size of the halo, so it should be the real r_vir divided by a_l.
  //-- Therefore, the correct equation for Takada & Jain (2003a) is
  //-- kappa_TJ = (4 pi G / c^2) * (f_l f_ls / f_s a_l) * M * (f_NFW c_NFW^2 a_l^2 / 2 pi r_vir^2) * G_TJ_kappa(c_NFW theta / theta_vir)
  //-- 
  //-- Since (f_l f_ls) / (f_s a_l) = (D_l D_ls) / (D_s a_l^2), one finds
  //-- kappa_TJ = (4 pi G / c^2) * (D_l D_ls / D_s a_l^2) * M * (f_NFW c_NFW^2 a_l^2 / 2 pi r_vir^2) * G_TJ_kappa(c_NFW theta / theta_vir)
  //--          = (4 pi G / c^2) * (D_l D_ls / D_s) * M * (f_NFW c_NFW^2 / 2 pi r_vir^2) * G_TJ_kappa(c_NFW theta / theta_vir)
  //--          = Sigma / Sigma_crit
  //-- Notice also that their Sigma_m is dimensionless, so that Sigma = M * Sigma_m.
  //-- 
  //-- For computation purposes, we write
  //-- kappa_TJ = (4 pi G / c^2) * (D_l D_ls M f_NFW c_NFW^2 / 2 pi D_s r_vir^2) * G_TJ_kappa(c_NFW theta / theta_vir)
  //--          = (4 pi G / c^2) * (D_l M f_NFW c_NFW^2 / 2 pi r_vir^2) * (D_ls / D_s) * G_TJ_kappa(c_NFW theta / theta_vir)
  //--          = factor * (D_ls / D_s) * G_TJ_kappa(c_NFW theta / theta_vir)
  //-- where factor = FOUR_PI_G_OVER_C2 * (D_l M f_NFW c_NFW^2 / 2 pi r_vir^2)
  
  double u_sq  = h->c_sq * theta_sq / h->theta_vir_sq; //-- u^2 = theta^2 / theta_s^2
  double D_ls  = g->a * f_K(cmhm->cosmo, g->w-h->w, err);                          forwardError(*err, __LINE__, -1.0);
  double kappa = h->factor * D_ls / g->D_s * G_TJ_kappa(u_sq, h->c, h->c_sq, err); forwardError(*err, __LINE__, -1.0);
  return kappa;
}

void gamma_TJ(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- See Takada & Jain (2003b), Eqs. (16) and (18)
  //-- |gamma_TJ| = R / Sigma_crit
  //-- where
  //--   Sigma_crit = (c^2 / 4 pi G) * (D_s / D_l D_ls), where (4 pi G / c^2) = 3 / (2 CRITICAL_DENSITY D_H^2 ) = 6.01e-19 [Mpc/M_sol]
  //--   R = M * (f c_NFW^2 / 2 pi r_vir^2) * G_TJ_gamma(c_NFW theta / theta_vir)
  //-- 
  //-- For computation purposes, we write
  //-- |gamma_TJ| = (4 pi G / c^2) * (D_l D_ls M f c_NFW^2 / 2 pi D_s r_vir^2) * G_TJ_gamma(c_NFW theta / theta_vir)
  //--            = factor * (D_ls / D_s) * G_TJ_kappa(c_NFW theta / theta_vir)
  //-- where factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  //--
  //-- gamma1  = -|gamma| cos(2 phi) = -|gamma| (1 - t^2) / (1 + t^2)
  //-- gamma2  = -|gamma| sin(2 phi) = -|gamma|       2t  / (1 + t^2)
  //-- where t = tan(phi)
  
  double u_sq       = h->c_sq * theta_sq / h->theta_vir_sq; //-- u^2 = theta^2 / theta_s^2
  double D_ls       = g->a * f_K(cmhm->cosmo, g->w-h->w, err);                                forwardError(*err, __LINE__,);
  double gamma_norm = h->factor * D_ls / g->D_s * G_TJ_gamma(u_sq, h->c, h->c_sq, h->f, err); forwardError(*err, __LINE__,);
  
  double t     = tan(atan2(g->pos[1]-h->pos[1], g->pos[0]-h->pos[0]));
  double t_sq  = SQ(t);
  double tcos  = (1.0 - t_sq) / (1.0 + t_sq);
  double tsin  =       2 * t  / (1.0 + t_sq);
  g->gamma[0] += -gamma_norm * tcos;
  g->gamma[1] += -gamma_norm * tsin;
  return;
}

void both_TJ(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  double u_sq   = h->c_sq * theta_sq / h->theta_vir_sq; //-- u^2 = theta^2 / theta_s^2
  double G_TJ[2];
  G_TJ_both(u_sq, h->c, h->c_sq, h->f, G_TJ, err);         forwardError(*err, __LINE__,);
  
  double D_ls   = g->a * f_K(cmhm->cosmo, g->w-h->w, err); forwardError(*err, __LINE__,);
  double factor = h->factor * D_ls / g->D_s;
  g->kappa     += factor * G_TJ[0]; //-- kappa_TJ
  G_TJ[1]      *= factor;           //-- gamma_norm
  
  double t     = tan(atan2(g->pos[1]-h->pos[1], g->pos[0]-h->pos[0]));
  double t_sq  = t * t;
  double tcos  = (1.0 - t_sq) / (1.0 + t_sq);
  double tsin  =       2 * t  / (1.0 + t_sq);
  g->gamma[0] += -G_TJ[1] * tcos;
  g->gamma[1] += -G_TJ[1] * tsin;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to drawing a profile

void fillOneHaloTerm(cosmo_hm *cmhm, halo_t *h, gal_t *g, double_mat *profile, error **err)
{
  int N_gal     = profile->N1;
  double D_ls   = g->a * f_K(cmhm->cosmo, g->w-h->w, err); forwardError(*err, __LINE__,);
  double factor = h->factor * D_ls / g->D_s;
  double theta_sq, u_sq, tau, G_halo[2];
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta_sq = SQ(profile->matrix[i]);
    u_sq = h->c_sq * theta_sq / h->theta_vir_sq;    //-- u^2 = theta^2 / theta_s^2
    
    tau = 2.5 * h->c; //-- Value suggested by Oguri & Hamana (2011)
    profile->matrix[i+N_gal]   = factor * G_BMO_kappa(u_sq, tau, SQ(tau), err); forwardError(*err, __LINE__,); //-- kappa_halo
    profile->matrix[i+2*N_gal] = 0.0; //-- gamma_norm
    
    //-- For NFW
    //G_NFW_both(u_sq, h->c, h->c_sq, h->f, G_halo, err); forwardError(*err, __LINE__,);
    //profile->matrix[i+N_gal]   = factor * G_halo[0]; //-- kappa_halo
    //profile->matrix[i+2*N_gal] = factor * G_halo[1]; //-- gamma_norm
    
    if (i % 125 == 0) {
      printf("Computing the one-halo term: %6.2f%% \r", 100.0 * i / (double)N_gal);
      fflush(stdout);
    }
  }
  printf("One-halo term done                       \n");
  return;
}

double factorForTwoHaloTerm(cosmo_hm *cmhm, halo_t *h, gal_t *g, error **err)
{
  //-- Oguri & Hamana (2011), Eq. (13)
  //-- kappa_2h(theta) = \int dl [(l/2pi) * factor * J_0(l*theta) * P_L(k_l, z)]
  //-- where k_l = l / f_K(w(z))
  //-- where factor = (rho_bar(z) b_h(M)) / ((1+z)^3 Sigma_crit D_l^2)
  //--              = (rho_bar(z) / (1+z)^3) * b_h(M) / ((c^2 / 4 pi G) * (D_s / D_l D_ls) * D_l^2)
  //--              = rho_bar(0) * b_h(M) / ((c^2 / 4 pi G) * (D_s D_l / D_ls))
  //--              = (CRITICAL_DENSITY * Omega_m) * b_h(M) * FOUR_PI_G_OVER_C2 * D_ls / (D_s D_l)
  //--              = (FOUR_PI_G_OVER_C2 CRITICAL_DENSITY Omega_m b_h(M) D_ls) / (D_s D_l)
  
  double D_l    = h->a * f_K(cmhm->cosmo, h->w, err);      forwardError(*err, __LINE__, -1.0);
  double D_ls   = g->a * f_K(cmhm->cosmo, g->w-h->w, err); forwardError(*err, __LINE__, -1.0);
  double b_h    = halo_bias(cmhm, h->M, h->a, 2, err);     forwardError(*err, __LINE__, -1.0); //-- order = 2 for bias_sc
  double factor = FOUR_PI_G_OVER_C2 * CRITICAL_DENSITY * cmhm->cosmo->Omega_m * b_h * D_ls / (g->D_s * D_l);
  return factor;
}

double integrandForTwoHaloTerm(double l, void *inteParam, error **err)
{
  //-- Oguri & Hamana (2011), Eq. (13)
  //-- kappa_2h(theta) = \int dl [(l / 2 pi) * factor * J_0(l*theta) * P_L(k_l, z)]
  //--                 = (factor / 2 pi) * \int dl [l * J_0(l*theta) * P_L(k_l, z)]
  //-- where k_l = l / f_K(w(z))
  //-- where factor = (rho_bar(z) b_h(M)) / ((1+z)^3 Sigma_crit D_l^2)
  //-- The integrand is l * J_0(l*theta) * P_L(k_l, z)
  
  profile_inteParam *param = (profile_inteParam*)inteParam;
  double k_l = l / param->f_K;
  double value = l * bessj0(l*param->theta) * P_L(param->cosmo, param->a, k_l, err); forwardError(*err, __LINE__, -1.0);
  return value;
}

double twoHaloTerm(cosmo_hm *cmhm, halo_t *h, double theta, double factor, error **err)
{
  profile_inteParam *inteParam = (profile_inteParam*)malloc_err(sizeof(profile_inteParam), err); forwardError(*err, __LINE__, -1.0);
  inteParam->theta = theta;
  inteParam->a     = h->a;
  inteParam->f_K   = f_K(cmhm->cosmo, h->w, err); forwardError(*err, __LINE__, -1.0);
  inteParam->cosmo = cmhm->cosmo;
  double l_maxx    = 500.0;
  double eps       = 1e-9;
  double totFactor = 0.5 * PI_INV * factor;
  double sum       = sm2_qromberg(integrandForTwoHaloTerm, (void*)inteParam, 0.0, l_maxx, eps, err) * totFactor;
  free(inteParam);
  return sum;
}

void fillTwoHaloTerm(cosmo_hm *cmhm, halo_t *h, double factor, double_mat *profile, error **err)
{
  int N_gal = profile->N1;
  double theta;
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta = profile->matrix[i] * ARCMIN_TO_RADIAN; //-- Conversion
    profile->matrix[i+3*N_gal] = twoHaloTerm(cmhm, h, theta, factor, err); forwardError(*err, __LINE__,); //-- 2-halo term kappa
    if (i % 25 == 0) {
      printf("Computing the two-halo term: %6.2f%% \r", 100.0 * i / (double)N_gal);
      fflush(stdout);
    }
  }
  printf("Two-halo term done                       \n");
  return;
}

void outputProfile(char name[], cosmo_hm *cmhm, peak_param *peak, halo_t *h, double_mat *profile, error **err)
{
  FILE *file = fopen(name, "w");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  fprintf(file, "#     theta           r          kappa        gamma           g      kappa_2h\n");
  fprintf(file, "#   [arcmin]      [Mpc/h]          [-]          [-]          [-]          [-]\n");
  
  int N_gal  = profile->N1;
  double D_l = h->a * f_K(cmhm->cosmo, h->w, err); forwardError(*err, __LINE__,);
  double theta, kappa, gamma, kappa_2h;
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta    = profile->matrix[i];
    kappa    = profile->matrix[i+N_gal];
    gamma    = profile->matrix[i+2*N_gal];
    kappa_2h = profile->matrix[i+3*N_gal];
    fprintf(file, " %11.5e  %11.5e  %11.5e  %11.5e  % 11.5e  %11.5e\n", theta, theta*ARCMIN_TO_RADIAN*D_l, kappa, gamma, gamma/(1.0-kappa), kappa_2h);
  }
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to lensing

void lensingForPair(cosmo_hm *cmhm, halo_t *h, gal_t *g, int doKappa, error **err)
{
  //-- Lensing for a halo-galaxy pair
  testErrorRet(g==NULL, peak_null, "Empty galaxy", *err, __LINE__,);
  
  if (g->z <= h->z) return; //-- Source in front of lens
  double theta_sq = DIST_2D_SQ(h->pos, g->pos);
  if (doKappa == 1) {
    //-- Compute projected kappa
    if (theta_sq > h->theta_vir_sq) return; //-- Too far
    g->kappa += kappa_TJ(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
  }
  else if (doKappa >= 2) {
    //-- Compute projected kappa and gamma
    if (theta_sq > CUTOFF_FACTOR_HALO_SQ * h->theta_vir_sq) return; //-- Too far
    both_TJ(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
  }
  else {
    //-- Compute projected gamma
    if (theta_sq > CUTOFF_FACTOR_HALO_SQ * h->theta_vir_sq) return; //-- Too far
    gamma_TJ(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
  }
  return;
}

void lensingForHalo(cosmo_hm *cmhm, halo_t *h, gal_map *gMap, int doKappa, error **err)
{
  //-- Lensing for a halo-galaxy pair
  testErrorRet(h==NULL, peak_null, "Empty halo", *err, __LINE__,);
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  double theta_pix_inv = 1.0 / gMap->theta_pix;
  double thetaInPix_x  = h->pos[0] * theta_pix_inv;
  double thetaInPix_y  = h->pos[1] * theta_pix_inv;
  double cutInPix_halo = h->theta_vir * theta_pix_inv;
  if (doKappa != 1) cutInPix_halo *= CUTOFF_FACTOR_HALO;
  
  int i_min = (int)(thetaInPix_x - cutInPix_halo);
  int i_max = (int)(thetaInPix_x + cutInPix_halo);
  int j_min = (int)(thetaInPix_y - cutInPix_halo);
  int j_max = (int)(thetaInPix_y + cutInPix_halo);
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int jN1;
  
  int i, j, k;
  for (j=j_min; j<=j_max; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_min; i<=i_max; i++) {
      if (i < 0 || i >= N1) continue;
      gList = gMap->map[i+jN1];
      for (k=0, gNode=gList->first; k<gList->size; k++, gNode=gNode->next) {
	lensingForPair(cmhm, h, gNode->g, doKappa, err);
	forwardError(*err, __LINE__,);
      }
    }
  }
  return;
}

void lensingForMap(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err)
{
  //-- Lensing main function
  
  int count = 0;
  halo_list *hList;
  halo_node *hNode;
  int i, j;
  
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      if ((peak->printMode == 0) && (count % 50 == 0)) {
	printf("Computing lensing signal: %6.2f%% \r", 100.0 * count / (double)hMap->total);
	fflush(stdout);
      }
      lensingForHalo(cmhm, hNode->h, gMap, peak->doKappa, err); forwardError(*err, __LINE__,);
      count++;
    }
  }
  
  if (peak->printMode < 2) printf("Lensing signal computation done           \n");
  if (peak->doKappa == 1)      gMap->type = kappa_map;
  else if (peak->doKappa >= 2) gMap->type = g_map;
  else                         gMap->type = gamma_map;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mask


short_mat *initializeMask(peak_param *peak, error **err)
{
  short_mat *CCDMask;
  if (peak->doMask == 0) {
    CCDMask = initialize_short_mat(2, 2);
    return CCDMask;
  }
  
#ifndef __releaseMenu__
  FITS_t *fits = initializeImageReader_FITS_t(peak->maskPath);
  
  peak->theta_CCD_inv[0] = 1.0 / (60.0 * readKeyword_double(fits, "CD1_1")); //-- [arcmin^-1]
  peak->theta_CCD_inv[1] = 1.0 / (60.0 * readKeyword_double(fits, "CD2_2"));
  int N_x      = readKeyword_int(fits, "NAXIS1");
  int N_y      = readKeyword_int(fits, "NAXIS2");
  int N1_mask  = (int)ceil(peak->Omega[0] * peak->theta_CCD_inv[0] - EPS_MIN);
  int N2_mask  = (int)ceil(peak->Omega[1] * peak->theta_CCD_inv[1] - EPS_MIN);
  testErrorRet(N1_mask > N_x || N2_mask > N_y, peak_geometry, "Mask too small", *err, __LINE__, NULL);
  
  CCDMask = initialize_short_mat(N1_mask, N2_mask);
  
  int limits[4];
  limits[0] = (N_x - N1_mask) / 2;
  limits[1] = limits[0] + N1_mask;
  limits[2] = (N_y - N2_mask) / 2;
  limits[3] = limits[2] + N2_mask;
  
  short *image = (short*)read2DSubImage(fits, limits);
  int i;
  for (i=0; i<CCDMask->length; i++) CCDMask->matrix[i] = image[i];
  free(image);
  free_FITS_t(fits);
  return CCDMask;
  
#else
  CCDMask = initialize_short_mat(2, 2);
  return CCDMask;
#endif
}

#define THETA_PIX 1.666666666666667e-2 //-- [arcmin]
#define THETA_PIX_INV 60.0 //-- [arcmin^-1]

#define W1_OMEGA_X 508 //-- [arcmin]
#define W1_OMEGA_Y 450 //-- [arcmin]
#define W1_I_MIN  1200 //-- [pix]
#define W1_I_MAX 31700 //-- [pix]
#define W1_J_MIN   800 //-- [pix]
#define W1_J_MAX 27800 //-- [pix]
short_mat *initializeMask_CFHTLenS_W1(peak_param *peak, error **err)
{
  short_mat *CCDMask;
  if (peak->doMask == 0) {
    CCDMask = initialize_short_mat(2, 2);
    return CCDMask;
  }
  
#ifndef __releaseMenu__
  FITS_t *fits = initializeImageReader_FITS_t("../../../91_Data/cfhtlens/mask/W1.16bit.small.reg2.fits");
  
  double Omega_x = peak->Omega[0];
  double Omega_y = peak->Omega[1];
  testErrorRet(Omega_x > W1_OMEGA_X || Omega_y > W1_OMEGA_Y, peak_geometry, "Mask too small", *err, __LINE__, NULL);
  
  peak->theta_CCD_inv[0] = THETA_PIX_INV;
  peak->theta_CCD_inv[1] = THETA_PIX_INV;
  int N1_mask = (int)ceil(Omega_x * THETA_PIX_INV - EPS_MIN);
  int N2_mask = (int)ceil(Omega_y * THETA_PIX_INV - EPS_MIN);
  CCDMask = initialize_short_mat(N1_mask, N2_mask);
  
  int limits[4];
  limits[0] = (W1_I_MIN + W1_I_MAX - N1_mask) / 2;
  limits[1] = limits[0] + N1_mask;
  limits[2] = (W1_J_MIN + W1_J_MAX - N2_mask) / 2;
  limits[3] = limits[2] + N2_mask;
  
  short *image = (short*)read2DSubImage(fits, limits);
  int i;
  for (i=0; i<CCDMask->length; i++) CCDMask->matrix[i] = image[i];
  free(image);
  free_FITS_t(fits);
  return CCDMask;
  
#else
  CCDMask = initialize_short_mat(2, 2);
  return CCDMask;
#endif
}
#undef W1_OMEGA_X
#undef W1_OMEGA_Y
#undef W1_I_MIN
#undef W1_I_MAX
#undef W1_J_MIN
#undef W1_J_MAX

#define W3_OMEGA_X 378 //-- [arcmin]
#define W3_OMEGA_Y 387 //-- [arcmin]
#define W3_I_MIN  1900 //-- [pix]
#define W3_I_MAX 24600 //-- [pix]
#define W3_J_MIN  1200 //-- [pix]
#define W3_J_MAX 24400 //-- [pix]
short_mat *initializeMask_CFHTLenS_W3(peak_param *peak, error **err)
{
  short_mat *CCDMask;
  if (peak->doMask == 0) {
    CCDMask = initialize_short_mat(2, 2);
    return CCDMask;
  }
  
#ifndef __releaseMenu__
  FITS_t *fits = initializeImageReader_FITS_t("../../../91_Data/cfhtlens/mask/W3.16bit.small.reg2.fits");
  
  double Omega_x = peak->Omega[0];
  double Omega_y = peak->Omega[1];
  testErrorRet(Omega_x > W3_OMEGA_X || Omega_y > W3_OMEGA_Y, peak_geometry, "Mask too small", *err, __LINE__, NULL);
  
  peak->theta_CCD_inv[0] = THETA_PIX_INV;
  peak->theta_CCD_inv[1] = THETA_PIX_INV;
  int N1_mask = (int)ceil(Omega_x * THETA_PIX_INV - EPS_MIN);
  int N2_mask = (int)ceil(Omega_y * THETA_PIX_INV - EPS_MIN);
  CCDMask = initialize_short_mat(N1_mask, N2_mask);
  
  int limits[4];
  limits[0] = (W3_I_MIN + W3_I_MAX - N1_mask) / 2;
  limits[1] = limits[0] + N1_mask;
  limits[2] = (W3_J_MIN + W3_J_MAX - N2_mask) / 2;
  limits[3] = limits[2] + N2_mask;
  
  short *image = (short*)read2DSubImage(fits, limits);
  int i;
  for (i=0; i<CCDMask->length; i++) CCDMask->matrix[i] = image[i];
  free(image);
  free_FITS_t(fits);
  return CCDMask;
  
#else
  CCDMask = initialize_short_mat(2, 2);
  return CCDMask;
#endif
}
#undef W3_OMEGA_X
#undef W3_OMEGA_Y
#undef W3_I_MIN
#undef W3_I_MAX
#undef W3_J_MIN
#undef W3_J_MAX

#undef THETA_PIX
#undef THETA_PIX_INV

short inCCDMask(short_mat *CCDMask, double pos[2], double theta_CCD_inv[2])
{
  int i = (int)(pos[0] * theta_CCD_inv[0]);
  int j = (int)(pos[1] * theta_CCD_inv[1]);
  return CCDMask->matrix[i+j*CCDMask->N1];
}

//----------------------------------------------------------------------
//-- Functions related to making galaxies

void makeRegularGalaxies(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  //-- Regular implies necessairily fixed redshift.
  
  double *Omega = peak->Omega;
  int M1 = (int)round(Omega[0] * sqrt(peak->n_gal));
  int M2 = (int)round(Omega[1] * sqrt(peak->n_gal));
  double pos[2];
  
  int i, j;
  for (j=0; j<M2; j++) {
    pos[1] = (j + 0.5) / M2 * Omega[1];
    for (i=0; i<M1; i++) {
      pos[0] = (i + 0.5) / M1 * Omega[0];
      append_gal_map(cmhm, gMap, peak->z_s, peak->w_s, peak->D_s, pos, err); forwardError(*err, __LINE__,);
    }
  }
  return;
}

void setGalaxySampler(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, error **err)
{
  galSamp->dx   = peak->dz_gal;
  galSamp->x[0] = get_zmin(cmhm->redshift, 0);
  int i;
  for (i=1; i<galSamp->length; i++) galSamp->x[i] = galSamp->x[i-1] + galSamp->dx;
  
  //-- Structure for redshift law
  redshiftANDint *rANDi = (redshiftANDint*)malloc_err(sizeof(redshiftANDint), err); forwardError(*err, __LINE__,);
  rANDi->self = cmhm->redshift;
  rANDi->i    = 0;
  fillGalaxyLaw(rANDi, galSamp, err); forwardError(*err, __LINE__,);
  free(rANDi);
  return;
}

void fillGalaxyLaw(redshiftANDint *rANDi, sampler_t *galSamp, error **err)
{
  double *z   = galSamp->x;
  double *pdf = galSamp->pdf;
  double *cdf = galSamp->cdf;
  int i;
  
  //-- Fill pdf
  for (i=0; i<galSamp->length; i++) {
    pdf[i] = prob_unnorm(z[i], (void*)rANDi, err);
    forwardError(*err, __LINE__,);
  }
  
  int setTotalToOne = 1;
  set_sampler_t(galSamp, setTotalToOne);
  return;
}

void makeRandomGalaxies(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, gal_map *gMap, error **err)
{
  gsl_rng *generator = peak->generator;
  int N_gal  = (int)round(peak->area * peak->n_gal);
  double z_s = peak->z_s;
  double w_s = peak->w_s;
  double D_s = peak->D_s;
  double p, z, pos[2];
 
  int i;
  if (peak->z_s > 0) {
    for (i=0; i<N_gal; i++) {
      randomizePosition(peak, pos, err);                   forwardError(*err, __LINE__,);
      append_gal_map(cmhm, gMap, z_s, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
    }
  }
  else {
    //-- Loop for generating galaxies
    for (i=0; i<N_gal; i++) { //-- Should change if cmhm->redshift->Nzbin > 1
      p = gsl_ran_flat(generator, 0.0, 1.0);
      z = execute_sampler_t(galSamp, p);
      randomizePosition(peak, pos, err);                   forwardError(*err, __LINE__,);
      append_gal_map(cmhm, gMap, z, -1.0, -1.0, pos, err); forwardError(*err, __LINE__,);
    }
  }
  
  if      (peak->printMode < 2)   printf("%d galaxies generated, no mask\n", N_gal);
  else if (peak->printMode == 2) {printf("%d galaxies, ", N_gal); fflush(stdout);}
  return;
}

void makeMaskedRandomGalaxies(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask, error **err)
{
  gsl_rng *generator = peak->generator;
  int N_gal  = (int)round(peak->area * peak->n_gal);
  double z_s = peak->z_s;
  double w_s = peak->w_s;
  double D_s = peak->D_s;
  double *theta_CCD_inv = peak->theta_CCD_inv;
  double p, z, pos[2];
  
  int i;
  if (peak->z_s > 0) {
    for (i=0; i<N_gal; i++) {
      randomizePosition(peak, pos, err);                   forwardError(*err, __LINE__,);
      if (inCCDMask(CCDMask, pos, theta_CCD_inv) > 0) continue;
      append_gal_map(cmhm, gMap, z_s, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
    }
  }
  else {
    //-- Loop for generating galaxies
    for (i=0; i<N_gal; i++) { //-- Should change if cmhm->redshift->Nzbin > 1
      p = gsl_ran_flat(generator, 0.0, 1.0);
      z = execute_sampler_t(galSamp, p);
      randomizePosition(peak, pos, err);                   forwardError(*err, __LINE__,);
      if (inCCDMask(CCDMask, pos, theta_CCD_inv) > 0) continue;
      append_gal_map(cmhm, gMap, z, -1.0, -1.0, pos, err); forwardError(*err, __LINE__,);
    }
  }
  
  if      (peak->printMode < 2)   printf("%d galaxies not masked, %d in total\n", gMap->total, N_gal);
  else if (peak->printMode == 2) {printf("%d galaxies, ", gMap->total); fflush(stdout);}
  return;
}

void cleanOrMakeOrResample(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask, error **err)
{
  if (peak->doRandGalPos == 1) {
    reset_gal_map(gMap);
    if (peak->doMask == 1) {makeMaskedRandomGalaxies(cmhm, peak, galSamp, gMap, CCDMask, err); forwardError(*err, __LINE__,);}
    else                   {makeRandomGalaxies(cmhm, peak, galSamp, gMap, err);                forwardError(*err, __LINE__,);}
  }
  else if (gMap->total > 0) cleanLensing_gal_map(gMap);
  else                     {makeRegularGalaxies(cmhm, peak, gMap, err);                        forwardError(*err, __LINE__,);}
  return;
}

//----------------------------------------------------------------------
//-- Functions related to galaxy catalogue

void subtractMean(peak_param *peak, gal_map *gMap)
{
  double mean = 0.0;
  gal_list *gList;
  gal_node *gNode;
  int i, j;
  
  //-- Get mean
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      mean += gNode->g->kappa;
    }
  }
  
  mean /= (double)gMap->total;
  gMap->kappa_mean = mean;
  
  //-- Subtract
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      gNode->g->kappa -= mean;
    }
  }
  
  if      (peak->printMode < 2)   printf("kappa mean %g subtracted\n", mean);
  else if (peak->printMode == 2) {printf("kappa mean = %.5f, ", mean); fflush(stdout);}
  return;
}

void makeG(gal_map *gMap)
{
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  double factor;
  int i, j;
  
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      g = gNode->g;
      factor       = 1.0 / (1.0 - g->kappa);
      g->gamma[0] *= factor;
      g->gamma[1] *= factor;
    }
  }
  
  gMap->type = g_map;
  return;
}

void addNoiseToGalaxies(peak_param *peak, gal_map *gMap)
{
  double sigma_half  = peak->sigma_half;
  gsl_rng *generator = peak->generator;
  
  gal_list *gList;
  gal_node *gNode;
  double *gamma;
  double g1, g2, e1, e2, A, B, C, g1_sq, g2_sq, factor;
  int i, j;
  
  if (peak->doKappa == 1) {
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	gNode->g->kappa += gsl_ran_gaussian(generator, sigma_half);
      }
    }
  }
  else if (peak->doKappa >= 2) {
    //-- Let e be source ellipticity, epsilon be observed ellipticity
    //-- epsilon = (e + g) / (1 + g^* e)
    //-- where
    //--   e   = e1 + i e2
    //--   g   = g1 + i g2
    //--   g^* = g1 - i g2
    //--
    //-- If x + iy = (a + ib) / (c + id)
    //-- then x = (ac + bd) / (c^2 + d^2)
    //--      y = (bc - ad) / (c^2 + d^2)
    //-- Here,
    //--   a = e1 + g1
    //--   b = e2 + g2
    //--   c = 1 + g1 e1 + g2 e2
    //--   d = g1 e2 - g2 e1
    //-- such that,
    //--   ac + bd   = e1 + g1 + g1(e1^2 + e2^2) + 2 g1 g2 e2 + e1(g1^2 - g2^2)
    //--   bc - ad   = e2 + g2 + g2(e1^2 + e2^2) + 2 g1 g2 e1 - e2(g1^2 - g2^2)
    //--   c^2 + d^2 = 1 + 2 g1 e1 + 2 g2 e2 + (g1^2 + g2^2)(e1^2 + e2^2)
    //--   
    //-- So,
    //--   epsilon1 = [e1 + g1 + g1 C + g1 B + e1(g1^2 - g2^2)] / [1 + A + B + C(g1^2 + g2^2)]
    //--   epsilon2 = [e2 + g2 + g2 C + g2 A - e2(g1^2 - g2^2)] / [1 + A + B + C(g1^2 + g2^2)]
    //-- where
    //--   A = 2 g1 e1
    //--   B = 2 g2 e2
    //--   C = e1^2 + e2^2
    //-- Finally,
    //--   epsilon1 = [e1(1 + g1^2 - g2^2) + g1(1 + C + B)] / [1 + A + B + C(g1^2 + g2^2)]
    //--   epsilon2 = [e2(1 - g1^2 + g2^2) + g2(1 + C + A)] / [1 + A + B + C(g1^2 + g2^2)]
    
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	gamma = gNode->g->gamma;
	g1 = gamma[0];
	g2 = gamma[1];
	e1 = gsl_ran_gaussian(generator, sigma_half);
	e2 = gsl_ran_gaussian(generator, sigma_half);
	A = 2.0 * g1 * e1;
	B = 2.0 * g2 * e2;
	C = e1 * e1 + e2 * e2;
	g1_sq = g1 * g1;
	g2_sq = g2 * g2;
	factor = 1.0 / (1.0 + A + B + (g1_sq + g2_sq) * C);
	gamma[0] = factor * (e1 * (1.0 + g1_sq - g2_sq) + g1 * (1.0 + C + B));
	gamma[1] = factor * (e2 * (1.0 - g1_sq + g2_sq) + g2 * (1.0 + C + A));
      }
    }
  }
  else {
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	gamma = gNode->g->gamma;
	gamma[0] += gsl_ran_gaussian(generator, sigma_half);
	gamma[1] += gsl_ran_gaussian(generator, sigma_half);
      }
    }
  }
  
  gMap->type += 2;
  return;
}

void lensingCatalogue(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err)
{
  //-- Lensing
  lensingForMap(cmhm, peak, hMap, gMap, err);
  forwardError(*err, __LINE__,);
  
  //-- Subtract mean
  if (peak->doKappa != 0) subtractMean(peak, gMap);
  
  //-- gamma to g
  if (peak->doKappa >= 2) makeG(gMap);
  
  //-- Add noise
  if (peak->doNoise == 1) {
    addNoiseToGalaxies(peak, gMap);
    if (peak->printMode < 2) printf("Added noise to galaxies\n");
  }
  
  return;
}

void lensingCatalogueAndOutputAll(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err)
{
  //-- Lensing
  lensingForMap(cmhm, peak, hMap, gMap, err);
  forwardError(*err, __LINE__,);
  
  //-- Subtract mean
  if (peak->doKappa != 0) subtractMean(peak, gMap);
  
  //-- gamma to g
  if (peak->doKappa >= 2) makeG(gMap);
  outputGalaxies("galCat", cmhm, peak, gMap);
  
  //-- Add noise
  if (peak->doNoise == 1) {
    addNoiseToGalaxies(peak, gMap);
    if (peak->printMode < 2) printf("Added noise to galaxies\n");
    outputGalaxies("galCat_noisy", cmhm, peak, gMap);
  }
  
  return;
}

void outputGalaxies(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Galaxy catalogue\n");
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], randomize positions = %d\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->doRandGalPos);
  fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file, "#\n");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  
  output_gal_map(file, peak, gMap);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Main function

void doRayTracing(char fileName[], cosmo_hm *cmhm, peak_param *peak, int doNoise, error **err)
{
  peak->doNoise = doNoise;
  
  halo_map *hMap     = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                     forwardError(*err, __LINE__,);
  gal_map *gMap      = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask = initializeMask(peak, err);                                                 forwardError(*err, __LINE__,);
  
  if (fileName == NULL) {
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->N_M);
    setMassSamplers(cmhm, peak, sampArr, err);                    forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err);                forwardError(*err, __LINE__,);
    outputFastSimul("haloCat", cmhm, peak, hMap);
    free_sampler_arr(sampArr);
  }
  else read_halo_map(fileName, cmhm, hMap, err);                  forwardError(*err, __LINE__,);
  
  cleanOrMakeOrResample(cmhm, peak, galSamp, gMap, CCDMask, err); forwardError(*err, __LINE__,);
  lensingCatalogueAndOutputAll(cmhm, peak, hMap, gMap, err);      forwardError(*err, __LINE__,);
  //outputGalaxies("galCat", cmhm, peak, gMap);
  
  free_halo_map(hMap);
  free_sampler_t(galSamp);
  free_gal_map(gMap);
  free_short_mat(CCDMask);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doProfile(char fileName[], cosmo_hm *cmhm, peak_param *peak, double z_l, double M, double z_s,  error **err)
{
  printf("z_l           = %f\n", z_l);
  printf("log(Mh/M_sol) = %f\n", log10(M));
  printf("z_s           = %f\n", z_s);
  
  peak->z_halo_max = z_s;
  peak->dz_halo    = peak->z_s / (double)(peak->N_z_halo);
  
  double theta_maxx = 200.0; //-- [arcmin]
  double dtheta     = 0.01;  //-- [arcmin]
  int N_gal         = (int)round(theta_maxx / dtheta);
  double pos[2]     = {0.0, 0.0};
  int i;
  
  halo_t *h       = (halo_t*)malloc_err(sizeof(halo_t), err); forwardError(*err, __LINE__,);
  set_halo_t(cmhm, h, z_l, M, pos, err);                      forwardError(*err, __LINE__,);
  
  gal_t *g        = (gal_t*)malloc_err(sizeof(gal_t), err);   forwardError(*err, __LINE__,);
  set_gal_t(cmhm, g, z_s, -1, -1, NULL, err);                 forwardError(*err, __LINE__,);
  double_mat *profile = initialize_double_mat(N_gal, 4);
  for (i=0; i<N_gal; i++) profile->matrix[i] = (i + 1) * dtheta;
  double factor = factorForTwoHaloTerm(cmhm, h, g, err);      forwardError(*err, __LINE__,);
  
  fillOneHaloTerm(cmhm, h, g, profile, err);                  forwardError(*err, __LINE__,);
  fillTwoHaloTerm(cmhm, h, factor, profile, err);             forwardError(*err, __LINE__,);
  outputProfile(fileName, cmhm, peak, h, profile, err);       forwardError(*err, __LINE__,);
  
  free(h);
  free(g);
  free_double_mat(profile);
  printf("------------------------------------------------------------------------\n");
  return;
}
void lensingCatalogueAndOutputAll2(char fileName[],cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err)
{
  //-- Lensing
  lensingForMap(cmhm, peak, hMap, gMap, err); forwardError(*err, __LINE__,);
  
  //-- Subtract mean
  if (peak->doKappa != 0) subtractMean(peak, gMap);
  
  //-- gamma to g
  if (peak->doKappa >= 2) makeG(gMap);
  outputGalaxies(fileName, cmhm, peak, gMap);
  
  //-- Add noise
  //if (peak->doNoise == 1) {
  //  addNoiseToGalaxies(peak, gMap);
  //  if (peak->printMode < 2) printf("Added noise to galaxies\n");
  //  strcat(fileName,"_noisy");
  //  outputGalaxies(fileName, cmhm, peak, gMap);
  //}
  
  return;
}


//----------------------------------------------------------------------



void read_gal_map2(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  //-- WARNING: D_S, w_s are not read.
  
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  printf("Reading...\r");
  fflush(stdout);
  
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double pos[2], gamma[2], z, kappa;
  
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      buffer2 = fscanf(file, "%lf %lf %lf \n", &pos[0], &pos[1], &z );   
      appendWithSignal_gal_map2(cmhm, gMap, z, pos, err); forwardError(*err, __LINE__,);
      count++;
    }
    c = fgetc(file);
  }
  
  fclose(file);
  testErrorRet(count!=gMap->total, peak_match, "Galaxy number match error", *err, __LINE__,);
  printf("\"%s\" read       \n", name);

  return;
}
void appendWithSignal_gal_map2(cosmo_hm *cmhm, gal_map *gMap, double z, double pos[2], error **err)
{
  double theta_pix_inv = gMap->theta_pix_inv;
  int i = (int)(pos[0] * theta_pix_inv);
  int j = (int)(pos[1] * theta_pix_inv);
  if (i >= gMap->N1) i -= 1;
  if (j >= gMap->N2) j -= 1;
  appendWithSignal_gal_list2(cmhm, gMap->map[i+j*gMap->N1], z, pos, err);
  forwardError(*err, __LINE__,);
  gMap->total++;
  return;
}
void appendWithSignal_gal_list2(cosmo_hm *cmhm, gal_list *gList, double z, double pos[2], error **err)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node(err);                            forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node(err);                       forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  setWithSignal_gal_t2(cmhm, gList->last->g, z, pos, err);
  forwardError(*err, __LINE__,);
  gList->size++;
  return;
}

void setWithSignal_gal_t2(cosmo_hm *cmhm, gal_t *g, double z, double pos[2], error **err)
{
  //-- WARNING: No w and D_s, only used in read_gal_map and fast pipeline
  testErrorRet(g==NULL, peak_null, "gal_t *g is NULL", *err, __LINE__,);
  
  g->z        = z;
  g->a        = 1.0/(1.0+z);
  g->pos[0]   = pos[0];
  g->pos[1]   = pos[1];
  g->w;
  g->D_s;
  return;
}


//----- new func-----

void outputFastSimul_galaxies(char name_cmhm[], char name[], char name2[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{
  error *myerr = NULL, **err = &myerr;
  gal_list *gList = initialize_gal_list(err); forwardError(*err, __LINE__,);
  FILE *file = fopen(name, "w");
  FILE *file2 = fopen(name2, "w");

  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");

  fprintf(file2, "# Galaxy catalogue\n");
  fprintf(file2,  "# Field = %s, Omega = (%g, %g) [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file2, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file2, "#\n");
  outputCosmoParam(file2, cmhm, peak);
  fprintf(file2, "#\n");

  //printf("test3 \n");
  output_halo_map_galaxies(file,file2,cmhm, peak, hMap, gList);
  //printf("test4 \n");
  free_gal_list(gList);
  //printf("Gfreen \n");
  fclose(file);
  fclose(file2);
  //printf("close \n");
  printf("\"%s\" made\n", name);
  printf("\"%s\" made\n", name2);
  return;
}



/// A revoir BEUG
void output2_halo_map_galaxies(FILE *file,FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_list *gList)
{
  halo_list *hList;
  halo_node *hNode;
  error *myerr = NULL, **err = &myerr;
  int i,j,k;

  srand(time(NULL));
  
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");
  
  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file,"#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv   \n");
    fprintf(file, "# [deg]  [deg]    [Mpc/h]     [-]     [M_sol/h]      [-]       [-]     [arcmin]    \n");
  }
  else {
    fprintf(file, "#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv \n");
    fprintf(file,"# [arcmin]  [arcmin]    [Mpc/h]    [-]     [M_sol/h]       [-]       [-]     [arcmin]  \n");
  }
  fprintf(file2, "#\n");
  fprintf(file2, "#\n");
  
 
  fprintf(file2, "#  theta_x   theta_y      z      halo_id    \n  ");
  fprintf(file2,"# [arcmin]  [arcmin]      [-]       [-]      \n");
  
  //printf("test5 \n");

  halo_t *h;

  cmhm->zmin=0.01 ;
		
  double vol_z=vc(cmhm,cmhm->zmin,cmhm->zmax,err);
  double deg2rad=3.14/180.;
  double arc2dec=1./60.;
  double ng= peak->n_gal/(arc2dec*deg2rad)/(arc2dec*deg2rad);
  double n_gal_obs=ng/vol_z ;

  //printf("n_gal_obs volz %9.3E %9.3E   \n",n_gal_obs,vol_z);
  double Mmin = pow(10.,10.)*pow(n_gal_obs,-0.24);
  double M1 = 17.*Mmin;
  cmhm->log10M1=log10(M1) ;
  cmhm->log10M_min = log10(Mmin) ;


  //printf("l = %i  \n",hMap->length);
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    //printf("s = %i  \n",hList->size);
    //printf("i = %i  \n",i);
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      h=hNode->h;
      double Mh=h->M;
      double Ds;

		//printf("Mmin %9.3E    \n",Mmin);
        double ngc = Ngal_c(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
        forwardError(*err, __LINE__,);
        double ngs = Ngal_s(cmhm, Mh, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
        forwardError(*err, __LINE__,);
      
      fprintf(file, "%9.3f  %9.3f    %8.3f  %7.5f  %9.3e   %8.3f  %8.3f   %9.3f  \n", h->pos[0], h->pos[1], h->w, h->z, Mh,ngc,ngs,h->r_vir);

      Ds  = h->a * f_K(cmhm->cosmo, h->w, err);
      
      if( (rand()/(double)RAND_MAX)<ngc) {
	append_gal_list(cmhm, gList,h->z, h->w, Ds, h->pos, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", h->pos[0], h->pos[1], h->z, j);
      }
  
    for (k = 0;k<ngc*ngs+0.5;k++) {

	int bool = 1;
	double r;
	while(bool){
	  double rtest = (rand()/(double)RAND_MAX);
	  if( NFW(5*rtest) > (rand()/(double)RAND_MAX)) {
	    bool = 0;
	    r = rtest*h->r_vir;
	  }
	}
	double theta = 2*M_PI*(rand()/(double)RAND_MAX);
	double phi = acos(2*(rand()/(double)RAND_MAX)-1);
	double pos[2];
        pos[0] = cos(theta) * sin(phi) * r + h->pos[0];
        pos[1] = sin(theta) * sin(phi) * r + h->pos[1];
	append_gal_list(cmhm, gList, h->z, h->w, Ds,h->pos, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", pos[0], pos[1], h->z, j);
      }
      
      //printf("ok \n");
    }
  }
  return;
}



void output_halo_map_galaxies(FILE *file,FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_list *gList)
{
  halo_list *hList;
  halo_node *hNode;
  error *myerr = NULL, **err = &myerr;
  int i,j,ii,k;
  double Ds,Mh;
  ii=0;
  srand(time(NULL));
  fprintf(file, "# Number of halos = %d\n", hMap->total);
  fprintf(file, "#\n");

  if (peak->field == aardvark_hPatch04 || peak->field == aardvark_gPatch086) { //-- For aardvark, positions are RA, DEC in [deg]
    fprintf(file,"#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv   \n");
    fprintf(file, "# [deg]  [deg]    [Mpc/h]     [-]     [M_sol/h]      [-]       [-]     [arcmin]    \n");
  }
  else {
    fprintf(file, "#  theta_x   theta_y      w          z          M         Ngal_c    Ngal_s      Rv \n");
    fprintf(file,"# [arcmin]  [arcmin]    [Mpc/h]    [-]     [M_sol/h]       [-]       [-]     [arcmin]  \n");
  }
  fprintf(file2, "#\n");
  fprintf(file2, "#\n");
  
  fprintf(file2, "#  theta_x   theta_y      z        \n  ");
  fprintf(file2,"# [arcmin]  [arcmin]      [-]       \n");

  halo_t *h;

  //printf("l = %i  \n",hMap->length);
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      h=hNode->h;
 	  Mh=h->M;
        double ngc = Ngal_c(cmhm,h->M, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
        forwardError(*err, __LINE__,);
        double ngs = Ngal_s(cmhm,h->M, cmhm->log10Mstar_min, cmhm->log10Mstar_max, err);
        forwardError(*err, __LINE__,);
        ii=ii+ngc*(1.0+ngs);
      fprintf(file, "%9.3f  %9.3f    %8.3f  %7.5f  %9.3e   %8.3f  %8.3f   %9.3f  \n", h->pos[0], h->pos[1], h->w, h->z, Mh,ngc,ngs,h->r_vir);

       Ds  = h->a * f_K(cmhm->cosmo, h->w, err);
      if( (rand()/(double)RAND_MAX)<ngc) {
	append_gal_list(cmhm, gList,h->z, h->w, Ds, h->pos, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f   \n", h->pos[0], h->pos[1], h->z);
      }
     
      for (k = 0;k<ngc*ngs+0.5;k++) {
	//printf("k = %i  \n",k);
	int bool = 1;
	double r;
	while(bool){
	  double rtest = (rand()/(double)RAND_MAX);
	  if( NFW(5*rtest) > (rand()/(double)RAND_MAX)) {
	    bool = 0;
	    r = rtest*h->r_vir;
	  }
	}
	double theta = 2*M_PI*(rand()/(double)RAND_MAX);
	double phi = acos(2*(rand()/(double)RAND_MAX)-1);
	double pos[2];
    pos[0] = cos(theta) * sin(phi) * r + h->pos[0];
    pos[1] = sin(theta) * sin(phi) * r + h->pos[1];
	append_gal_list(cmhm, gList, h->z, h->w, Ds,h->pos, err); forwardError(*err, __LINE__,);
	fprintf(file2, "%9.3f  %9.3f   %7.5f   \n", pos[0], pos[1], h->z);
      }
    }
  }
  printf("Nb galaxies created : %i \n",ii);
  return;
}


double NFW(double x)
{
  return 1/(x*(1+pow(x,2)));
}


