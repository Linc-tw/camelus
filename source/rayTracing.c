

  /***************************************************************
   **  rayTracing.c						**
   **  Chieh-An Lin, Martin Kilbinger, François Lanusse		**
   **  Version 2015.12.09					**
   **								**
   **  References:						**
   **  - Bartelmann & Schneider (2001) - Phys. Rep., 340, 291	**
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
    int wOmegar = 0;
    g->w   = w(cmhm->cosmo, a, wOmegar, err); forwardError(*err, __LINE__,);
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
  gal_map *gMap          = (gal_map*)malloc_err(sizeof(gal_map), err);                    forwardError(*err, __LINE__,);
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
  
  gMap->map              = (gal_list**)malloc_err(gMap->length * sizeof(gal_list*), err); forwardError(*err, __LINE__,);
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
  threshold *= 0.5;                                                 //-- Set threshold to the half of average
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
    int wOmegar = 0;
    w_s  = w(cmhm->cosmo, 1.0/(1.0 + peak->z_s), wOmegar, err); forwardError(*err, __LINE__,);
    D_s  = f_K(cmhm->cosmo, w_s, err);                          forwardError(*err, __LINE__,);
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
//-- Functions related to projected NFW mass

#define EPS 1e-10
double G_exactNFW_kappa(double x_sq, double c, double c_sq, error **err)
{
  //-- See Wright & Brainerd (2000), Eq. (11)
  //-- Sigma = 2 rho_s r_s G(x)
  //-- G(x) = -1/arg3 + 2 artanh(arg2) / arg3^(3/2), if x < 1;
  //--         1/3,                                  if x = 1;
  //--        +1/arg3 - 2 arctan(arg2) / arg3^(3/2), if x > 1;
  //-- where
  //--   arg2 = sqrt(|1 - x| / (1 + x)) = sqrt(arg3) / (1 + x),
  //--   arg3 = |1 - x^2|.
  
  if (x_sq < EPS) x_sq = EPS;
  
  double arg3 = fabs(x_sq - 1.0);
  double value;
  
  if (arg3 < EPS) value = 1.0 / 3.0;
  else {
    double arg2 = sqrt(arg3) / (1.0 + sqrt(x_sq));
    if (x_sq < 1.0) {
      value  = -1.0 + 2.0 * atanh(arg2) / sqrt(arg3);
      value /= arg3;
    }
    else {
      value  = 1.0 - 2.0 * atan(arg2) / sqrt(arg3);
      value /= arg3;
    }
  }
  
  return value;
}

double G_truncatedNFW_kappa(double x_sq, double c, double c_sq, error **err)
{
  //-- See Takada & Jain (2003a), Eq. (27)
  //-- G(x) = -arg1/arg3 + arcosh(arg2) / arg3^(3/2), if x < 1;
  //--        1/3 * arg1 * (1 + 1 / (c + 1)),         if x = 1;
  //--        +arg1/arg3 - arccos(arg2) / arg3^(3/2), if 1 < x <= c;
  //--        0,                                      if x > c;
  //-- where
  //--   arg1 = sqrt(c^2 - x^2) / (c + 1),
  //--   arg2 = (x^2 + c) / (x(c + 1)),
  //--   arg3 = |1 - x^2|.
  
  if (x_sq > c_sq - EPS) return 0;
  if (x_sq < EPS) x_sq = EPS;
  
  double p_inv = 1.0 / (c + 1.0);
  double arg1  = sqrt(c_sq - x_sq) * p_inv;
  double arg3  = fabs(x_sq - 1.0);
  double value;
  
  if (arg3 < EPS) value = arg1 * (1.0 + p_inv) / 3.0;
  else {
    double arg2 = (x_sq + c) * p_inv / sqrt(x_sq);
    if (x_sq < 1.0) {
      testErrorRet(arg2<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__, 0.0);
      value  = -arg1 + acosh(arg2) / sqrt(arg3);
      value /= arg3;
    }
    else {
      testErrorRet(arg2>1.0, peak_badValue, "Out of range for acos", *err, __LINE__, 0.0);
      value  = arg1 - acos(arg2) / sqrt(arg3);
      value /= arg3;
    }
  }
  
  return value;
}

double G_exactNFW_gamma(double x_sq, double c, double c_sq, double f, error **err)
{
  //-- See Wright & Brainerd (2000), Eqs. (15) and (16)
  //-- G(x) = 1/arg3 * (1 - arg1) + 2 x^-2 (arg1 + arg4),   if x < 1;
  //--        5 / 3 + 2 ln(1/2),                             if x = 1;
  //--        1/arg3 * (-1 + arg5) + 2 x^-2 (arg5 + arg4)], if x > 1;
  //-- where
  //--   arg1 = 2 artanh(arg2) / r,
  //--   arg2 = sqrt(|1 - x| / (1 + x)) = r / (1 + x),
  //--   arg3 = |1 - x^2|,
  //--   arg4 = ln(x / 2),
  //--   arg5 = 2 arctan(arg2) / r,
  //-- and
  //--   r = sqrt(|1 - x^2|) = sqrt(arg3).
  
  if (x_sq < EPS) x_sq = EPS;
  
  double arg3 = fabs(x_sq - 1.0);
  double value;
  
  if (arg3 < EPS) {
    value = 5.0 / 3.0 + 2.0 * log(0.5);
  }
  else {
    double x    = sqrt(x_sq);
    double r    = sqrt(arg3);
    double arg2 = r / (1.0 + x);
    double arg4 = log(0.5 * x);
    if (x_sq < 1.0) {
      double arg1 = 2.0 * atanh(arg2) / r;
      value = (1.0 - arg1) / arg3 + 2.0 / x_sq * (arg1 + arg4);
    }
    else {
      double arg5 = 2.0 * atan(arg2) / r;
      value = (-1.0 + arg5) / arg3 + 2.0 / x_sq * (arg5 + arg4);
    }
  }
  
  return value;
}

double G_truncatedNFW_gamma(double x_sq, double c, double c_sq, double f, error **err)
{
  //-- See Takada & Jain (2003b), Eq. (17)
  //-- G(x) = x^-2 [arg1 + 2 ln(arg2) + arg3 * arcosh(arg4)], if x < 1;
  //--        x^-2 [arg5 + 2 ln(arg2)],                       if x = 1;
  //--        x^-2 [arg1 + 2 ln(arg2) + arg3 * arccos(arg4)], if 1 < x <= c;
  //--        x^-2 * 2/f,                                     if x > c;
  //-- where
  //--   arg1 = [(q + 1) r / q - 2c] / p,
  //--   arg2 = x p / (c + r),
  //--   arg3 = (3q - 1) / q sqrt|q|,
  //--   arg4 = (x^2 + c) / x p,
  //--   arg5 = [(c + 10p) r / p - 6c] / 3p,
  //-- and
  //--   p = c + 1,
  //--   q = 1 - x^2,
  //--   r = sqrt(c^2 - x^2).
  
  if (x_sq > c_sq - EPS) return 2.0 / (f * x_sq);
  if (x_sq < EPS) x_sq = EPS;
  
  double p     = c + 1.0;
  double p_inv = 1.0 / (c + 1.0);
  double q     = 1.0 - x_sq;
  double r     = sqrt(c_sq - x_sq);
  double x     = sqrt(x_sq);
  double arg2  = x * p / (c + r);
  double value;
  
  if (fabs(q) < EPS) {
    value  = ((10.0 + c * p_inv) * r - 6*c) * p_inv / 3.0 + 2 * log(arg2);
    value /= x_sq;
    return value;
  }
  
  double q_inv = 1.0 / q;
  double arg1  = ((1.0 + q_inv) * r - 2 * c) * p_inv;
  double arg3  = (3.0 - q_inv) / sqrt(fabs(q));
  double arg4  = (x_sq + c) * p_inv / x;
  
  if (q > 0) {
    testErrorRet(arg4<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__, 0.0);
    value  = arg1 + 2 * log(arg2) + arg3 * acosh(arg4);
    value /= x_sq;
  }
  else {
    testErrorRet(arg4>1.0, peak_badValue, "Out of range for acos", *err, __LINE__, 0.0);
    value  = arg1 + 2 * log(arg2) + arg3 * acos(arg4);
    value /= x_sq;
  }
  
  return value;
}

void G_exactNFW_both(double x_sq, double c, double c_sq, double f, double G_NFW[2], error **err)
{
  if (x_sq < EPS) x_sq = EPS;
  
  double arg3 = fabs(x_sq - 1.0);
  double value;
  
  if (arg3 < EPS) {
    G_NFW[0] = 1.0 / 3.0;
    G_NFW[1] = 5.0 / 3.0 + 2.0 * log(0.5);
  }
  else {
    double x    = sqrt(x_sq);
    double r    = sqrt(arg3);
    double arg2 = r / (1.0 + x);
    double arg4 = log(0.5 * x);
    if (x_sq < 1.0) {
      double arg1 = 2.0 * atanh(arg2) / r;
      G_NFW[0] = (-1.0 + arg1) / arg3;
      G_NFW[1] = -G_NFW[0] + 2.0 / x_sq * (arg1 + arg4);
    }
    else {
      double arg5 = 2.0 * atan(arg2) / r;
      G_NFW[0] = (1.0 - arg5) / arg3;
      G_NFW[1] = -G_NFW[0] + 2.0 / x_sq * (arg5 + arg4);
    }
  }
  
  return;
}

void G_truncatedNFW_both(double x_sq, double c, double c_sq, double f, double G_NFW[2], error **err)
{
  if (x_sq > c_sq - EPS) {
    G_NFW[0] = 0.0;
    G_NFW[1] = 2.0 / (f * x_sq);
    return;
  }
  if (x_sq < EPS) x_sq = EPS;
  
  double p     = c + 1;
  double p_inv = 1.0 / p;
  double q     = 1 - x_sq;
  double r     = sqrt(c_sq - x_sq);
  double x     = sqrt(x_sq);
  double arg2  = x * p / (c + r);
  double arg5  = r * p_inv;
  double arg6  = fabs(q);
  double value;
  
  if (arg6 < EPS) {
    value     = p_inv / 3.0;
    G_NFW[0]  = r * (1 + p_inv) * value;
    G_NFW[1]  = ((10 + c * p_inv) * r - 6*c) * value + 2 * log(arg2);
    G_NFW[1] /= x_sq;
    return;
  }
  
  double q_inv = 1.0 / q;
  double arg1  = ((1 + q_inv) * r - 2 * c) * p_inv;
  double arg3  = 3 - q_inv;
  double arg4  = (x_sq + c) * p_inv / x;
  
  if (q > 0) {
    testErrorRet(arg4<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__,);
    value     = acosh(arg4) / sqrt(arg6);
    G_NFW[0]  = -arg5 + value;
    G_NFW[0] /= arg6;
    G_NFW[1]  = arg1 + 2 * log(arg2) + arg3 * value;
    G_NFW[1] /= x_sq;
  }
  else {
    testErrorRet(arg4>1.0, peak_badValue, "Out of range for acos", *err, __LINE__,);
    value     = acos(arg4) / sqrt(arg6);
    G_NFW[0]  = arg5 - value;
    G_NFW[0] /= arg6;
    G_NFW[1]  = arg1 + 2 * log(arg2) + arg3 * value;
    G_NFW[1] /= x_sq;
  }
  
  return;
}
#undef EPS

double kappa_truncatedNFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- See Bartelmann & Schneider (2001), Eq. (3.7)
  //-- kappa_truncatedNFW = Sigma / Sigma_crit
  //-- where
  //--   Sigma_crit = (c^2 / 4 pi G) * (D_s / D_l D_ls), where (4 pi G / c^2) = 3 / (2 CRITICAL_DENSITY D_H^2 ) = 6.01e-19 [Mpc/M_sol]
  //--   Sigma = projected mass
  //--         = 2 rho_s r_s * G_truncatedNFW_kappa(theta / theta_s)
  //--         = 2 * (rho_bar * Delta f c_NFW^3 / 3) * (r_vir / c_NFW) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //--         = 2 * (rho_bar * Delta * 4 pi r_vir^3 / 3) * (f c_NFW^3 / 4 pi r_vir^3) * (r_vir / c_NFW) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //--         = M * (f c_NFW^2 / 2 pi r_vir^2) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //-- 
  //-- Similarly, in Takada & Jain (2003a), Eqs. (26) and (28), one finds
  //-- kappa_truncatedNFW = W(w_l, w_s) * (M / rho_bar) * Sigma_m(theta)
  //--                    = (3 Omega_m H_0^2 / 2) * (D_l D_ls / D_s a_l) * (M / rho_bar) * Sigma_m
  //--                    = (4 pi G) * (D_l D_ls / D_s a_l) * M * (f c_NFW^2 / 2 pi R_vir^2) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //-- However, their convention is different. First, they have omitted c^(-2).
  //-- Second, their D_l, D_ls, D_s are actually f_K: f_l, f_ls, f_s.
  //-- Finally, their R_vir is defined as a comoving distance, not the physical size of the halo, so it should be the real r_vir divided by a_l.
  //-- Therefore, the correct equation for Takada & Jain (2003a) is
  //-- kappa_truncatedNFW = (4 pi G / c^2) * (f_l f_ls / f_s a_l) * M * (f c_NFW^2 a_l^2 / 2 pi r_vir^2) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //-- 
  //-- Since (f_l f_ls) / (f_s a_l) = (D_l D_ls) / (D_s a_l^2), one finds
  //-- kappa_truncatedNFW = (4 pi G / c^2) * (D_l D_ls / D_s a_l^2) * M * (f c_NFW^2 a_l^2 / 2 pi r_vir^2) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //--                    = (4 pi G / c^2) * (D_l D_ls / D_s) * M * (f c_NFW^2 / 2 pi r_vir^2) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //--                    = Sigma / Sigma_crit
  //-- Notice also that their Sigma_m is dimensionless, so that Sigma = M * Sigma_m.
  //-- 
  //-- For computation purposes, we write
  //-- kappa_truncatedNFW = (4 pi G / c^2) * (D_l D_ls M f c_NFW^2 / 2 pi D_s r_vir^2) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //--                    = (4 pi G / c^2) * (D_l M f c_NFW^2 / 2 pi r_vir^2) * (D_ls / D_s) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //--                    = factor * (D_ls / D_s) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //-- where factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  
  double x_sq  = h->c_sq * theta_sq / h->theta_vir_sq; //-- x^2 = theta^2 / theta_s^2
  double D_ls  = g->a * f_K(cmhm->cosmo, g->w-h->w, err);                           forwardError(*err, __LINE__, 0.0);
  double kappa = h->factor * D_ls / g->D_s * G_truncatedNFW_kappa(x_sq, h->c, h->c_sq, err); forwardError(*err, __LINE__, 0.0);
  return kappa;
}

void gamma_truncatedNFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- See TJ03b (16, 18)
  //-- |gamma_truncatedNFW| = R / Sigma_crit
  //-- where
  //--   Sigma_crit = (c^2 / 4 pi G) * (D_s / D_l D_ls), where (4 pi G / c^2) = 3 / (2 CRITICAL_DENSITY D_H^2 ) = 6.01e-19 [Mpc/M_sol]
  //--   R = M * (f c_NFW^2 / 2 pi r_vir^2) * G_truncatedNFW_gamma(c_NFW theta / theta_vir)
  //-- 
  //-- For computation purposes, we write
  //-- |gamma_truncatedNFW| = (4 pi G / c^2) * (D_l D_ls M f c_NFW^2 / 2 pi D_s r_vir^2) * G_truncatedNFW_gamma(c_NFW theta / theta_vir)
  //--                      = factor * (D_ls / D_s) * G_truncatedNFW_kappa(c_NFW theta / theta_vir)
  //-- where factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  //--
  //-- gamma1  = -|gamma| cos(2 phi) = -|gamma| (1 - t^2) / (1 + t^2)
  //-- gamma2  = -|gamma| sin(2 phi) = -|gamma|       2t  / (1 + t^2)
  //-- where t = tan(phi)
  
  double x_sq       = h->c_sq * theta_sq / h->theta_vir_sq; //-- x^2 = theta^2 / theta_s^2
  double D_ls       = g->a * f_K(cmhm->cosmo, g->w-h->w, err);                                 forwardError(*err, __LINE__,);
  double gamma_norm = h->factor * D_ls / g->D_s * G_truncatedNFW_gamma(x_sq, h->c, h->c_sq, h->f, err); forwardError(*err, __LINE__,);
  
  double t     = tan(atan2(g->pos[1]-h->pos[1], g->pos[0]-h->pos[0]));
  double t_sq  = SQ(t);
  double tcos  = (1.0 - t_sq) / (1.0 + t_sq);
  double tsin  =       2 * t  / (1.0 + t_sq);
  g->gamma[0] += -gamma_norm * tcos;
  g->gamma[1] += -gamma_norm * tsin;
  return;
}

void both_truncatedNFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- g = gamma / (1 - kappa)
  
  double x_sq   = h->c_sq * theta_sq / h->theta_vir_sq; //-- x^2 = theta^2 / theta_s^2
  double G_NFW[2];
  G_truncatedNFW_both(x_sq, h->c, h->c_sq, h->f, G_NFW, err); forwardError(*err, __LINE__,);
  
  double D_ls   = g->a * f_K(cmhm->cosmo, g->w-h->w, err);    forwardError(*err, __LINE__,);
  double factor = h->factor * D_ls / g->D_s;
  g->kappa     += factor * G_NFW[0]; //-- kappa_proj
  G_NFW[1]     *= factor;            //-- gamma_norm
  
  double t     = tan(atan2(g->pos[1]-h->pos[1], g->pos[0]-h->pos[0]));
  double t_sq  = t * t;
  double tcos  = (1.0 - t_sq) / (1.0 + t_sq);
  double tsin  =       2 * t  / (1.0 + t_sq);
  g->gamma[0] += -G_NFW[1] * tcos;
  g->gamma[1] += -G_NFW[1] * tsin;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to drawing a profile

void computeProfile(cosmo_hm *cmhm, halo_t *h, gal_t *g, double_mat *profile, error **err)
{
  int N_gal     = profile->N1;
  double D_ls   = g->a * f_K(cmhm->cosmo, g->w-h->w, err); forwardError(*err, __LINE__,);
  double factor = h->factor * D_ls / g->D_s;
  double theta_sq, x_sq, G_NFW[2];
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta_sq = SQ(profile->matrix[i]);
    x_sq = h->c_sq * theta_sq / h->theta_vir_sq; //-- x^2 = theta^2 / theta_s^2
    G_exactNFW_both(x_sq, h->c, h->c_sq, h->f, G_NFW, err); forwardError(*err, __LINE__,);
    profile->matrix[i+N_gal]   = factor * G_NFW[0]; //-- kappa_proj
    profile->matrix[i+2*N_gal] = factor * G_NFW[1]; //-- gamma_norm
  }
  return;
}

void outputProfile(char name[], cosmo_hm *cmhm, peak_param *peak, halo_t *h, double_mat *profile, error **err)
{
  FILE *file = fopen(name, "w");
  outputCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  fprintf(file, "#     theta           r          kappa        gamma           g\n");
  fprintf(file, "#   [arcmin]      [Mpc/h]          [-]          [-]          [-]\n");
  
  int N_gal  = profile->N1;
  double D_l = h->a * f_K(cmhm->cosmo, h->w, err); forwardError(*err, __LINE__,);
  double theta, kappa, gamma;
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta = profile->matrix[i];
    kappa = profile->matrix[i+N_gal];
    gamma = profile->matrix[i+2*N_gal];
    fprintf(file, " %11.5e  %11.5e  %11.5e  %11.5e  %11.5e\n", theta, theta*D_l, kappa, gamma, gamma/(1.0-kappa));
  }
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}
/*
double *wtheta_halo(cosmo_hm *model, pofk_t pofk, double *theta, int Ntheta, int i_bin, int j_bin, error **err)
{
  //-- Take the 3D 2PCF (xi), project it using Limber equation (at mean z ?), 
  //-- and return the angular correlation function (w(theta))
  
  //-- WARNING Below to verify
  //-- See Bartelmann & Schneider (2001) eq. (2.79)
  //-- Theta in degree
  //-- N is the number of points that sample xi. The speed of the code is basically inversely proportional to this number.
  
  testErrorRetVA(MIN(i_bin, j_bin)<0 || MAX(i_bin, j_bin)>=model->redshift->Nzbin, redshift_Nzbin,
		 "Requested redshift bins (%d, %d) out of range [0; %d]",
		 *err, __LINE__, NULL, i_bin, j_bin, model->redshift->Nzbin-1);
  
  double *result  = malloc_err(Ntheta*sizeof(double),err); 
  forwardError(*err, __LINE__, NULL);
  
  //-- tabulate xi(r)
  int i, N = 40;
  double umin = 0.001;
  double umax = 800;
  
  double *u    = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double *logu = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double dlogu = log(umax / umin) / (double)N;
  
  for(i=0; i<N; i++) {
    logu[i] = log(umin) + dlogu * (double)i;
    u[i]    = exp(logu[i]);
  }
  
  double *xi = xi_P_NL(model, a, r, N, err); forwardError(*err, __LINE__, NULL);
  //xiofr(model, pofk, u, N, GG, err); forwardError(*err, __LINE__, NULL);
  
  //-- interpolate xi(r)
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline    = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(spline, logu, xi, N);
  
  Nz_hist *nz1 = get_nz(model->redshift, i_bin, err); forwardError(*err, __LINE__, NULL);

  int wOmegar = 0;
  double deg_to_rad_sqr = dsqr(PI / 180.0);
  double ww, r, fK, sum, a;
  int j, k;
  double nsqr_dzdr;
  
  //-- Limber equation - project xi to get w
  for(i=0;i<Ntheta;i++){ //-- loop over theta
    result[i] = 0.0;
    
    //-- loop over z
    for (j=0; j<nz1->nbins; j++) {
      
      a   = 1.0 / (1.0 + nz1->z[j]);
      ww  = w(model->cosmo, a, wOmegar, err);                 forwardError(*err, __LINE__, NULL);
      fK  = f_K(model->cosmo, ww, err);                       forwardError(*err, __LINE__, NULL);
      sum = 0.0;
      
      //-- loop over u
      for (k=0; k<N; k++) { 
	r = sqrt(u[k] * u[k] + fK * fK * theta[i] * theta[i] * deg_to_rad_sqr);
	//-- to unsure log(r) lies within interpolation limits
	if (log(r) < logu[N-1]) sum += u[k] * gsl_spline_eval(spline, log(r), acc);
      }
      nsqr_dzdr = 1.0 / drdz(model->cosmo, a, err); forwardError(*err, __LINE__, NULL);
      result[i] += nsqr_dzdr * sum;
    }

    result[i] *= 2.0 * dlogu;
    testErrorRetVA(!finite(result[i]), ce_infnan, "inf or nan in w(theta_%d=%g)", *err, __LINE__, NULL, i, theta[i]/arcmin);
  }
  
  free(xi);
  free(u);
  free(logu);
  free_nz(nz1);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  
  return result;
}
*/
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
    g->kappa += kappa_truncatedNFW(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
  }
  else if (doKappa >= 2) {
    //-- Compute projected kappa and gamma
    if (theta_sq > CUTOFF_FACTOR_HALO_SQ * h->theta_vir_sq) return; //-- Too far
    both_truncatedNFW(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
  }
  else {
    //-- Compute projected gamma
    if (theta_sq > CUTOFF_FACTOR_HALO_SQ * h->theta_vir_sq) return; //-- Too far
    gamma_truncatedNFW(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
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

#define THETA_PIX 1.666666666666667e-2 //-- [arcmin]
#define THETA_PIX_INV 60.0 //-- [arcmin^-1]

#define W1_OMEGA_X 508 //-- [arcmin]
#define W1_OMEGA_Y 450 //-- [arcmin]
#define W1_I_MIN  1200 //-- [pix]
#define W1_I_MAX 31700 //-- [pix]
#define W1_J_MIN   800 //-- [pix]
#define W1_J_MAX 27800 //-- [pix]
void fillMask_CFHTLenS_W1(peak_param *peak, short_mat *CCDMask)
{
  peak->theta_CCD_inv = THETA_PIX_INV;
  
  double Omega_x = peak->Omega[0];
  double Omega_y = peak->Omega[1];
  if (Omega_x > W1_OMEGA_X || Omega_y > W1_OMEGA_Y) return;
  
  char name[STRING_LENGTH_MAX];
  sprintf(name, "../../../91_Data/CFHTLenS/mask/W1.16bit.small.reg2.fits");
  FITS_t *fits = initializeImageReader_FITS_t(name);
  
  int limits[4];
  limits[0] = (int)(0.5 * (W1_I_MIN + W1_I_MAX - Omega_x * THETA_PIX_INV));
  limits[1] = limits[0] + CCDMask->N1;
  limits[2] = (int)(0.5 * (W1_J_MIN + W1_J_MAX - Omega_y * THETA_PIX_INV));
  limits[3] = limits[2] + CCDMask->N2;
  
  short *image = (short*)read2DSubImage(fits, limits);
  int i;
  for (i=0; i<CCDMask->length; i++) CCDMask->matrix[i] = image[i];
  free(image);
  return;
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
void fillMask_CFHTLenS_W3(peak_param *peak, short_mat *CCDMask)
{
  peak->theta_CCD_inv = THETA_PIX_INV;
  
  double Omega_x = peak->Omega[0];
  double Omega_y = peak->Omega[1];
  if (Omega_x > W3_OMEGA_X || Omega_y > W3_OMEGA_Y) return;
  
  char name[STRING_LENGTH_MAX];
  sprintf(name, "../../../91_Data/CFHTLenS/mask/W3.16bit.small.reg2.fits");
  FITS_t *fits = initializeImageReader_FITS_t(name);
  
  int limits[4];
  limits[0] = (int)(0.5 * (W3_I_MIN + W3_I_MAX - Omega_x * THETA_PIX_INV));
  limits[1] = limits[0] + CCDMask->N1;
  limits[2] = (int)(0.5 * (W3_J_MIN + W3_J_MAX - Omega_y * THETA_PIX_INV));
  limits[3] = limits[2] + CCDMask->N2;
  
  short *image = (short*)read2DSubImage(fits, limits);
  int i;
  for (i=0; i<CCDMask->length; i++) CCDMask->matrix[i] = image[i];
  free(image);
  return;
}
#undef W3_OMEGA_X
#undef W3_OMEGA_Y
#undef W3_I_MIN
#undef W3_I_MAX
#undef W3_J_MIN
#undef W3_J_MAX

#undef THETA_PIX
#undef THETA_PIX_INV

short inCCDMask(short_mat *CCDMask, double pos[2], double theta_CCD_inv)
{
  int i = (int)(pos[0] * theta_CCD_inv);
  int j = (int)(pos[1] * theta_CCD_inv);
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
  double theta_CCD_inv = peak->theta_CCD_inv;
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
  lensingForMap(cmhm, peak, hMap, gMap, err); forwardError(*err, __LINE__,);
  
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
  lensingForMap(cmhm, peak, hMap, gMap, err); forwardError(*err, __LINE__,);
  
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
  int N1_mask = (int)(peak->Omega[0] * peak->theta_CCD_inv);
  int N2_mask = (int)(peak->Omega[1] * peak->theta_CCD_inv);
  peak->doNoise = doNoise;
  
  halo_map *hMap     = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *galSamp = initialize_sampler_t(peak->N_z_gal);
  setGalaxySampler(cmhm, peak, galSamp, err);                                                     forwardError(*err, __LINE__,);
  gal_map *gMap      = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  short_mat *CCDMask = initialize_short_mat(N1_mask, N2_mask);
  
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
  printf("z_l = %f\n", z_l);
  printf("M   = %f\n", log10(M));
  printf("z_s = %f\n", z_s);
  
  peak->z_halo_max = z_s;
  peak->dz_halo    = peak->z_s / (double)(peak->N_z_halo);
  
  /*
  int N_a = 4000;
  double sheet[2];
  sampler_t *samp = initialize_sampler_t(peak->N_M);
  interpolator_t *a_inter = initialize_interpolator_t(N_a+1);
  makeAInterpolator(cmhm, a_inter, peak->z_halo_max, err); forwardError(*err, __LINE__,);
  massSheet(cmhm, peak, samp, a_inter, sheet, err);        forwardError(*err, __LINE__,);
  printf("kappa_0 = %f\n", sheet[0]);
  free_sampler_t(samp);
  free_interpolator_t(a_inter);
  */
  
  double theta_maxx = 200.0; //-- [arcmin]
  double dtheta     = 0.01;  //-- [arcmin]
  int N_gal         = (int)round(theta_maxx / dtheta);
  double pos[2]     = {0.0, 0.0};
  int i;
  
  halo_t *h       = (halo_t*)malloc_err(sizeof(halo_t), err); forwardError(*err, __LINE__,);
  set_halo_t(cmhm, h, z_l, M, pos, err);                      forwardError(*err, __LINE__,);
  
  /*
  double c = 5.0;
  h->c        = c;
  h->c_sq     = SQ(c);
  h->f        = 1.0 / (log(1.0 + c) - c/(1.0 + c));
  */
  
  gal_t *g        = (gal_t*)malloc_err(sizeof(gal_t), err);   forwardError(*err, __LINE__,);
  set_gal_t(cmhm, g, z_s, -1, -1, NULL, err);                 forwardError(*err, __LINE__,);
  double_mat *profile = initialize_double_mat(N_gal, 3);
  for (i=0; i<N_gal; i++) profile->matrix[i] = (i + 1) * dtheta;
  computeProfile(cmhm, h, g, profile, err);                   forwardError(*err, __LINE__,);
  //for (i=N_gal; i<2*N_gal; i++) profile->matrix[i] -= sheet[0];
  outputProfile(fileName, cmhm, peak, h, profile, err);       forwardError(*err, __LINE__,);
  
  free(h);
  free(g);
  free_double_mat(profile);
  printf("------------------------------------------------------------------------\n");
  return;
}
//----------------------------------------------------------------------


