

  /*******************************************************
   **  galaxySampling.c					**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/

   
#include "galaxySampling.h"


//----------------------------------------------------------------------
//-- Functions related to gal_t, gal_node, gal_list

void set_gal_t(cosmo_hm *chPar, peak_param *pkPar, gal_t *g, double pos[2], double z, double weight, error **err)
{
  testErrorRet(g==NULL, peak_null, "gal_t *g is NULL", *err, __LINE__,);
  
  if (pos != NULL) {
    g->pos[0] = pos[0];
    g->pos[1] = pos[1];
    
    if (pkPar->field > 0) {
      g->sinCosDEC[0] = sin(g->pos[1]);
      g->sinCosDEC[1] = cos(g->pos[1]);
    }
  }
  
  g->z     = z;
  double a = 1.0/(1.0+z);
  g->a     = a;
  if (pkPar->w_s < 0) {
    //-- Use z to compute w_s and D_s
    g->w   = w(chPar->cosmo, a, 0, err);       forwardError(*err, __LINE__,); //-- wOmegar = 0
    g->D_s = a * f_K(chPar->cosmo, g->w, err); forwardError(*err, __LINE__,);
  }
  else {
    g->w   = pkPar->w_s;
    g->D_s = pkPar->D_s;
  }
  
  g->weight   = weight;
  g->kappa    = 0.0;
  g->gamma[0] = 0.0;
  g->gamma[1] = 0.0;
  return;
}

void setWithLensing_gal_t(peak_param *pkPar, gal_t *g, double pos[2], double z, double weight, double kappa, double gamma[2], error **err)
{
  //-- Only used in read_gal_map
  testErrorRet(g==NULL, peak_null, "gal_t *g is NULL", *err, __LINE__,);
  
  if (pos != NULL) {
    g->pos[0] = pos[0];
    g->pos[1] = pos[1];
    
    if (pkPar->field > 0) {
      g->sinCosDEC[0] = sin(g->pos[1]);
      g->sinCosDEC[1] = cos(g->pos[1]);
    }
  }
  g->z        = z;
  g->a        = 1.0/(1.0+z);
  g->weight   = weight;
  g->kappa    = kappa;
  g->gamma[0] = gamma[0];
  g->gamma[1] = gamma[1];
  //g->w;   //-- No need to set, since signal is already obtained
  //g->D_s; //-- No need to set, since signal is already obtained
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
  gal_list *gList  = (gal_list*)malloc_err(sizeof(gal_list), err); forwardError(*err, __LINE__, NULL);
  gList->length    = 0;
  gList->size      = 0;
  gList->totWeight = 0.0;
  gList->mean[0]   = 0.0;
  gList->mean[1]   = 0.0;
  gList->first     = NULL;
  gList->last      = NULL;
  return gList;
}

void free_gal_list(gal_list *gList)
{
  gal_node *gNode;
  if (gList) {
    while (gList->first != NULL) {
      gNode        = gList->first;
      gList->first = gNode->next;
      if (gNode->g) {free(gNode->g); gNode->g = NULL;}
      free(gNode); gNode = NULL;
    }
    free(gList); gList = NULL;
  }
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
  gList->size      = 0;
  gList->totWeight = 0.0;
  gList->mean[0]   = 0.0;
  gList->mean[1]   = 0.0;
  return;
}

void append_gal_list(cosmo_hm *chPar, peak_param *pkPar, gal_list *gList, double pos[2], double z, double weight, error **err)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node(err);      forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node(err); forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  set_gal_t(chPar, pkPar, gList->last->g, pos, z, weight, err); forwardError(*err, __LINE__,);
  gList->totWeight += weight;
  gList->size++;
  return;
}

void appendWithLensing_gal_list(peak_param *pkPar, gal_list *gList, double pos[2], double z, double weight, double kappa, double gamma[2], error **err)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node(err);      forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node(err); forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  setWithLensing_gal_t(pkPar, gList->last->g, pos, z, weight, kappa, gamma, err); forwardError(*err, __LINE__,);
  gList->totWeight += weight;
  gList->size++;
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
  
  gMap->kappa_mean       = 0.0;
  gMap->fillingThreshold = 0.0;
  gMap->type             = kappa_map;
  gMap->map              = (gal_list**)malloc_err(gMap->length * sizeof(gal_list*), err); forwardError(*err, __LINE__, NULL);
  int i;
  for (i=0; i<gMap->length; i++) gMap->map[i] = initialize_gal_list(err);
  return gMap;
}

void free_gal_map(gal_map *gMap)
{
  int i;
  if (gMap) {
    if (gMap->map) {
      for (i=0; i<gMap->length; i++) {free_gal_list(gMap->map[i]); gMap->map[i] = NULL;}
      free(gMap->map); gMap->map = NULL;
    }
    free(gMap); gMap = NULL;
  }
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

int append_gal_map(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, double pos[2], double z, double weight, error **err)
{
  double pos2[2];
  if (pkPar->field == 0) {
    pos2[0] = pos[0];
    pos2[1] = pos[1];
  }
  else {
    project(pos, pkPar->center, pos2);
    rotate(pos2, pkPar->rotAng, pos2);
    pos2[0] += pkPar->offset;
    pos2[1] += pkPar->offset;
  }
  
  int i = (int)(pos2[0] * gMap->theta_pix_inv);
  int j = (int)(pos2[1] * gMap->theta_pix_inv);
  if ((i < 0) || (j < 0) || (i >= gMap->N1) || (j >= gMap->N2)) return 0;
  append_gal_list(chPar, pkPar, gMap->map[i+j*gMap->N1], pos, z, weight, err); forwardError(*err, __LINE__, -1);
  gMap->total++;
  return 1;
}

int appendWithLensing_gal_map(peak_param *pkPar, gal_map *gMap, double pos[2], double z, double weight, double kappa, double gamma[2], error **err)
{
  double pos2[2];
  if (pkPar->field == 0) {
    pos2[0] = pos[0];
    pos2[1] = pos[1];
  }
  else {
    project(pos, pkPar->center, pos2);
    rotate(pos2, pkPar->rotAng, pos2);
    pos2[0] += pkPar->offset;
    pos2[1] += pkPar->offset;
  }
  
  int i = (int)(pos2[0] * gMap->theta_pix_inv);
  int j = (int)(pos2[1] * gMap->theta_pix_inv);
  if ((i < 0) || (j < 0) || (i >= gMap->N1) || (j >= gMap->N2)) return 0;
  appendWithLensing_gal_list(pkPar, gMap->map[i+j*gMap->N1], pos, z, weight, kappa, gamma, err); forwardError(*err, __LINE__, -1);
  gMap->total++;
  return 1;
}

void updateCosmo_gal_map(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err)
{
  if (pkPar->z_s > 0) {
    pkPar->w_s  = w(chPar->cosmo, 1.0/(1.0 + pkPar->z_s), 0, err); forwardError(*err, __LINE__,); //-- wOmegar = 0
    pkPar->D_s  = f_K(chPar->cosmo, pkPar->w_s, err);              forwardError(*err, __LINE__,);
    pkPar->D_s /= 1.0 + pkPar->z_s;
  }
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int i, j;
  
  //-- Loop for setting galaxies
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      g = gNode->g;
      set_gal_t(chPar, pkPar, g, NULL, g->z, g->weight, err); forwardError(*err, __LINE__,);
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to signal mean

void kappaMean_gal_list(gal_list *gList, double threshold)
{
  //-- Variance:
  //--   The variance on the first compoenent of the i-th galaxy is the sum of the shape noise and the measurement error.
  //--     sigma_i,1^2 = sigma_half^2 + sigma_m,i,1^2
  //--   Define:
  //--     sigma_i^2   = sigma_i,2^2 + sigma_i,2^2
  //--     sigma_m,i^2 = sigma_m,i,2^2 + sigma_m,i,2^2
  //-- 
  //-- Optimal weight:
  //--   The optimal weight of a galaxy is proportional to the inverse of the variance.
  //--     omega_i \propto sigma_i^-2
  //--   From Lensfit, the weight is set to
  //--     omega_i = 2 / sigma_i^2 or sigma_i^2 = 2 / omega_i
  //--
  //-- Effective number of galaxies:
  //--   The effective number of galaxies is the shape noise times the sum of variances.
  //--   sigma_half^2 / N_eff = sigma_tot^2 = (sum_i omega_i^2 sigma_i^2 / 2) / (sum_i omega_i)^2
  //--                                      = 1 / (sum_i omega_i)
  //--   So, N_eff = sigma_half^2 * (sum_i omega_i).
  //--   If no measurement error, omega_i = 1 / sigma_half^2, and N_eff = N_gal.
  //--
  //-- Total weight and threshold:
  //--   The filling factor is defined as the total weight in that pixel.
  //--   The total weight is proportional to the effective number of galaxies.
  //--   The filling threshold is set as a factor times the average filling factor.
  
  int size         = gList->size;
  double totWeight = gList->totWeight;
  double *mean     = gList->mean;
  if (size == 0 || totWeight < threshold) {
    mean[0] = 0.0;
    mean[1] = 0.0;
    return;
  }
  
  double mean1 = 0.0;
  gal_node *gNode;
  gal_t *g;
  int i;
  
  for (i=0, gNode=gList->first; i<size; i++, gNode=gNode->next) {
    g      = gNode->g;
    mean1 += g->kappa * g->weight;
  }
  mean[0] = mean1 / totWeight;
  mean[1] = 0.0;
  return;
}

void gammaOrGMean_gal_list(gal_list *gList, double threshold, int doSurvey)
{
  //-- See also kappaMean_gal_list for total weight and threshold
  //--
  //-- Mean with multiplicative correction:
  //--   Let m_corr,i be the mulpliticative correction of the i-th galaxy.
  //--   The one-point estimator of the shear is
  //--     <eps_1> = (sum_i omega_i eps_i,1) / (sum_i omega_i (1 + m_corr_i))
  //--   The two-point estimator of the shear is
  //--     <eps_1 eps_1> = (sum_ij omega_i omega_j eps_i,1 eps_j,1) / (sum_ij omega_i omega_j (1 + m_corr_i) (1 + m_corr_j))
  //--   If doSurvey = 0, there is no correction. The mean is a weighted sum.
  //--   If doSurvey = 1, g->kappa is supposed to stock 1 + m_corr, and the mean is computed following the formula above.
  
  int size         = gList->size;
  double totWeight = gList->totWeight;
  double *mean     = gList->mean;
  if (size == 0 || totWeight < threshold) {
    mean[0] = 0.0;
    mean[1] = 0.0;
    return;
  }
  
  double mean1 = 0.0;
  double mean2 = 0.0;
  double mean3 = 0.0;
  gal_node *gNode;
  gal_t *g;
  int i;
  
  if (doSurvey == 1) {
    for (i=0, gNode=gList->first; i<size; i++, gNode=gNode->next) {
      g      = gNode->g;
      mean1 += g->gamma[0] * g->weight;
      mean2 += g->gamma[1] * g->weight;
      mean3 += g->kappa    * g->weight; //-- weight * (1 + m_corr)
    }
    mean[0] = mean1 / mean3;
    mean[1] = mean2 / mean3;
  }
  
  else {
    for (i=0, gNode=gList->first; i<size; i++, gNode=gNode->next) {
      g      = gNode->g;
      mean1 += g->gamma[0] * g->weight;
      mean2 += g->gamma[1] * g->weight;
    }
    mean[0] = mean1 / totWeight;
    mean[1] = mean2 / totWeight;
  }
  return;
}

void setFillingThreshold(peak_param *pkPar, gal_map *gMap, error **err)
{
  //-- The filling factor is defined as the total weight in that pixel.
  //-- The total weight is proportional to the effective number of galaxies.
  //-- The filling threshold is set as a factor times the average filling factor.
  
  if (gMap->fillingThreshold != 0.0) return; //-- Already computed
  if (pkPar->doMask == 0) {
    gMap->fillingThreshold = EPS_NUM; //-- Should not set to zero or negative values, otherwise empty pixels yield +infty in the variance calculation
    if      (pkPar->verbose < 3)   printf("No mask\n");
    else if (pkPar->verbose == 3) {printf("no mask, "); fflush(stdout);}
    return;
  }
  
  int bufferSize = pkPar->bufferSize;
  testErrorRet(bufferSize<0, peak_badValue, "Buffer size should be at least 0.", *err, __LINE__,);
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  double threshold = 0.0;
  int i, j, jN1;
  
  //-- Set filling threshold
  for (j=bufferSize; j<N2-bufferSize; j++) { 
    jN1 = j * N1;
    for (i=bufferSize; i<N1-bufferSize; i++) threshold += gMap->map[i+jN1]->totWeight;
  }
  threshold /= (double)((N1 - 2*bufferSize) * (N2 - 2*bufferSize)); //-- Average filling factor
  threshold *= FILLING_THRESHOLD_RATIO;                             //-- Set threshold to the half of average
  gMap->fillingThreshold = threshold;
  
  //-- Print the filling threshold in terms of effective number of galaxies
  if      (pkPar->verbose < 3)   printf("Filling threshold = %.2f\n", threshold * SQ(pkPar->sigma_half));
  else if (pkPar->verbose == 3) {printf("threshold = %.2f, ", threshold * SQ(pkPar->sigma_half)); fflush(stdout);}
  return;
}

void signalMean_gal_map(peak_param *pkPar, gal_map *gMap, error **err)
{
  setFillingThreshold(pkPar, gMap, err); forwardError(*err, __LINE__,);
  int i;
  if (pkPar->doKappa == 1) {
    for (i=0; i<gMap->length; i++) kappaMean_gal_list(gMap->map[i], gMap->fillingThreshold);
    gMap->type = kappa_map;
  }
  else {
    for (i=0; i<gMap->length; i++) gammaOrGMean_gal_list(gMap->map[i], gMap->fillingThreshold, pkPar->doSurvey);
    if (pkPar->doKappa >= 2) gMap->type = g_map;
    else                     gMap->type = gamma_map;
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to mask

mask_map *initialize_mask_map(int N1, int N2, double theta_mask[2], double offset[2], double center[2], error **err)
{
  mask_map *mask          = (mask_map*)malloc_err(sizeof(mask_map), err);        forwardError(*err, __LINE__, NULL);
  mask->N1                = N1;
  mask->N2                = N2;
  mask->length            = N1 * N2;
  mask->nbChar            = (mask->length - 1) / 8 + 1;
  mask->theta_mask[0]     = theta_mask[0];
  mask->theta_mask[1]     = theta_mask[1];
  mask->theta_mask_inv[0] = 1.0 / theta_mask[0];
  mask->theta_mask_inv[1] = 1.0 / theta_mask[1];
  mask->offset[0]         = (offset == NULL) ? 0.0 : offset[0];
  mask->offset[1]         = (offset == NULL) ? 0.0 : offset[1];
  mask->center[0]         = (center == NULL) ? 0.0 : center[0];
  mask->center[1]         = (center == NULL) ? 0.0 : center[1];
  mask->type              = proj_mask;
  mask->map               = (char*)malloc_err(mask->nbChar * sizeof(char), err); forwardError(*err, __LINE__, NULL);
  return mask;
}

void free_mask_map(mask_map *mask)
{
  if (mask) {
    if (mask->map) {free(mask->map); mask->map = NULL;}
    free(mask); mask = NULL;
  }
  return;
}

mask_map *initializeMask(peak_param *pkPar, error **err)
{
  //-- RADEC_fieldCtr in [deg], can be NULL
  //--
  //-- For CFHTLenS:
  //--   theta_mask      = 0.01667 [arcmin]
  //--   theta_mask_inv  = 60.0 [arcmin^-1]
  //-- W1 field:
  //--   Mask size       = (32845, 28798) [pix]
  //--                   = (547, 480) [arcmin]
  //--   Effective limit = (1200, 31700, 800, 27800) [pix]
  //--   Effective size  = (508, 450) [arcmin]
  //-- W3 field:
  //--   Mask size       = (26454, 25282) [pix]
  //--                   = (440, 421) [arcmin]
  //--   Effective limit = (1900, 24600, 1200, 24400) [pix]
  //--   Effective size  = (378, 387) [arcmin]
  
  mask_map *mask;
  if (pkPar->field == 2) {
    mask = NULL;
    return mask;
  }
  
  if (pkPar->doMask == 0) {
    mask = NULL;
    return mask;
  }
  
  double theta_mask[2];
  int N1, N2;
  
  //-- Random mask
  if (pkPar->doMask == 1) {
    theta_mask[0] = ARCSEC_TO_RADIAN; //-- Pixel size of 1 arcsec in [rad]
    theta_mask[1] = ARCSEC_TO_RADIAN;
    N1 = (int)ceil(pkPar->Omega[0] / theta_mask[0]);
    N2 = (int)ceil(pkPar->Omega[1] / theta_mask[1]);
    mask = initialize_mask_map(N1, N2, theta_mask, NULL, NULL, err); forwardError(*err, __LINE__, NULL); //-- offset = NULL, center = NULL
    //-- Randomization to be done in cleanOrMakeOrResample
    return mask;
  }
  
#ifdef __CAMELUS_USE_FITS__
  //-- CFHTLenS-like settings
  FILE* file = fopen(pkPar->maskPath, "r");
  testErrorRet(file==NULL, peak_unknown, "Cannot open mask file", *err, __LINE__, NULL);
  fclose(file);
  FITS_t *fits = initializeImageReader_FITS_t(pkPar->maskPath);
  
  theta_mask[0] = readKeyword_double(fits, "CD1_1") * DEGREE_TO_RADIAN; //-- [rad]
  theta_mask[1] = readKeyword_double(fits, "CD2_2") * DEGREE_TO_RADIAN;
  //testErrorRet(fabs(theta_mask[0])!=fabs(theta_mask[1]), peak_geometry, "Mask pixel not squared", *err, __LINE__, NULL); //-- Check if pixel is square
  
  double thetaXY_fieldCtr[2] = {0.0, 0.0};
  double ij_maskCtr[2];
  double ij_fieldCtr[2];
  
//   double maskCtr[4];
//   if (RADEC_fieldCtr != NULL) {
//     maskCtr[0]  = readKeyword_double(fits, "CRVAL1"); //-- [deg]
//     maskCtr[1]  = readKeyword_double(fits, "CRVAL2");
//     maskCtr[0] *= DEGREE_TO_RADIAN; //-- [rad]
//     maskCtr[1] *= DEGREE_TO_RADIAN;
//     maskCtr[2]  = sin(maskCtr[1]);
//     maskCtr[3]  = cos(maskCtr[1]);
//     project(RADEC_fieldCtr, maskCtr, thetaXY_fieldCtr); //-- [rad]
//   }
  
  ij_maskCtr[0]  = readKeyword_double(fits, "CRPIX1"); //-- [pix]
  ij_maskCtr[1]  = readKeyword_double(fits, "CRPIX2");
  ij_fieldCtr[0] = ij_maskCtr[0] + thetaXY_fieldCtr[0] / theta_mask[0]; //-- [pix]
  ij_fieldCtr[1] = ij_maskCtr[1] + thetaXY_fieldCtr[1] / theta_mask[1];
  
  theta_mask[0] = fabs(theta_mask[0]);
  theta_mask[1] = fabs(theta_mask[1]);
  
  double i_min = ij_fieldCtr[0] - 0.5 * pkPar->Omega[0] / theta_mask[0];
  double i_max = ij_fieldCtr[0] + 0.5 * pkPar->Omega[0] / theta_mask[0];
  double j_min = ij_fieldCtr[1] - 0.5 * pkPar->Omega[1] / theta_mask[1];
  double j_max = ij_fieldCtr[1] + 0.5 * pkPar->Omega[1] / theta_mask[1];
  
  double offset[2];
  offset[0] = i_min - floor(i_min);
  offset[1] = j_min - floor(j_min);
  
  int limit[4];
  limit[0] = (int)i_min;
  limit[1] = (int)ceil(i_max);
  limit[2] = (int)j_min;
  limit[3] = (int)ceil(j_max);
  
  int N_x = readKeyword_int(fits, "NAXIS1");
  int N_y = readKeyword_int(fits, "NAXIS2");
  testErrorRet(limit[0]<0 || limit[1]>N_x || limit[2]<0 || limit[3]>N_y, peak_geometry, "Mask too small", *err, __LINE__, NULL);
  
  mask = initialize_mask_map(limit[1]-limit[0], limit[3]-limit[2], theta_mask, offset, NULL, err); forwardError(*err, __LINE__, NULL); //-- center = NULL
  readInputMask(fits, pkPar, mask, limit);
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Read mask from input\n");
#endif
  return mask;
}

void resetMask(peak_param *pkPar, mask_map *mask)
{
  int i, q, r;
  for (i=0; i<mask->length; i++) {
    q = i / 8;
    r = i % 8;
    mask->map[q] = CLEAN_BIT(mask->map[q], r); //-- 0 = activated, 1 = masked
  }
  return;
}

void generateHole(peak_param *pkPar, mask_map *mask, double theta, int holeInd, error **err)
{
  double pos[2];
  randomizePosition(pkPar, pos, err); forwardError(*err, __LINE__,);
  
  int N1 = mask->N1;
  int N2 = mask->N2;
  double x_pix = pos[0] * mask->theta_mask_inv[0];
  double y_pix = pos[1] * mask->theta_mask_inv[1];
  int i_pix    = (int)x_pix;
  int j_pix    = (int)y_pix;
  int cutSize  = (int)ceil(theta * mask->theta_mask_inv[0]);
  double thetaSqInPix = SQ(theta * mask->theta_mask_inv[0]);
  char *map = mask->map;
  
  int i, j, q, r, jN1, index;
  
  //-- Create holes
  for (j=j_pix-cutSize; j<=j_pix+cutSize; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_pix-cutSize; i<=i_pix+cutSize; i++) {
      if (i < 0 || i >= N1) continue;
      if (SQ(i-i_pix) + SQ(j-j_pix) <= thetaSqInPix) {
	index = i + jN1;
	q = index / 8;
	r = index % 8;
	map[q] = SET_BIT(map[q], r); //-- 0 = activated, 1 = masked
      }
    }
  }
  
  double factor1 = pkPar->stripeLRatio[holeInd];
  double factor2 = pkPar->stripeWRatio[holeInd];
  int i_min, i_max, j_min, j_max;
  
  //-- Create stripes
  i_min = (int)round(x_pix - theta * mask->theta_mask_inv[0] * factor1);
  i_max = (int)round(x_pix + theta * mask->theta_mask_inv[0] * factor1);
  j_min = (int)round(y_pix - theta * mask->theta_mask_inv[1] * factor2);
  j_max = (int)round(y_pix + theta * mask->theta_mask_inv[1] * factor2);
  
  for (j=j_min; j<=j_max; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_min; i<=i_max; i++) {
      if (i < 0 || i >= N1) continue;
      index = i + jN1;
      q = index / 8;
      r = index % 8;
      map[q] = SET_BIT(map[q], r); //-- 0 = activated, 1 = masked
    }
  }
  
  i_min = (int)round(x_pix - theta * mask->theta_mask_inv[0] * factor2);
  i_max = (int)round(x_pix + theta * mask->theta_mask_inv[0] * factor2);
  j_min = (int)round(y_pix - theta * mask->theta_mask_inv[1] * factor1);
  j_max = (int)round(y_pix + theta * mask->theta_mask_inv[1] * factor1);
  
  for (j=j_min; j<=j_max; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_min; i<=i_max; i++) {
      if (i < 0 || i >= N1) continue;
      index = i + jN1;
      q = index / 8;
      r = index % 8;
      map[q] = SET_BIT(map[q], r); //-- 0 = activated, 1 = masked
    }
  }
  
  return;
}

void makeRandomMask(peak_param *pkPar, mask_map *mask, error **err)
{
  resetMask(pkPar, mask);
  
  double theta, density;
  int i, j, nbHoles;
  
  for (i=0; i<pkPar->nbHoleTypes; i++) {
    theta   = pkPar->holeRadius[i];
    nbHoles = (int)round(pkPar->holeDensity[i] * pkPar->area);
    
    for (j=0; j<nbHoles; j++) {
      generateHole(pkPar, mask, theta, i, err);
      forwardError(*err, __LINE__,);
    }
  }
  
  if (pkPar->verbose < 3) printf("Generated random masks\n");
  return;
}

void readInputMask(FITS_t *fits, peak_param *pkPar, mask_map *mask, int limit[4])
{
  short *image = (short*)read2DSubImage(fits, limit);
  char *map = mask->map;
  int i, q, r;
  
  //-- CFHTLenS accepts flag = 0 or 1.
  //-- WARNING Need to think how to propagate from param file
  for (i=0; i<mask->length; i++) {
    q = i / 8;
    r = i % 8;
    map[q] = (image[i] < 2) ? CLEAN_BIT(map[q], r) : SET_BIT(map[q], r); //-- 0 = activated, 1 = masked
  }
  free(image);
  return;
}

int inMask(mask_map *mask, double pos[2])
{
  if (mask == NULL) return 0;
  int i = (int)(pos[0] * mask->theta_mask_inv[0] + mask->offset[0]);
  int j = (int)(pos[1] * mask->theta_mask_inv[1] + mask->offset[1]);
  int index_mask = i + j * mask->N1;
  int q = index_mask / 8;
  int r = index_mask % 8;
  int bit = (int)((mask->map[q] >> r) & 1);
  return bit;
}

//----------------------------------------------------------------------
//-- Functions related to I/O

void outAscii_gal_map(FILE *file, peak_param *pkPar, gal_map *gMap)
{
  fprintf(file, "# Type = %s, number of galaxies = %d\n", STR_MAP_T(gMap->type), gMap->total);
  fprintf(file, "#\n");
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int i, j;
  
  if (pkPar->field == 0) {
    if (gMap->type < 12) {
      fprintf(file, "#   theta_x    theta_y       z     weight        kappa      gamma_1      gamma_2\n");
      fprintf(file, "#  [arcmin]   [arcmin]      [-]       [-]          [-]          [-]          [-]\n");
    }
    else {
      fprintf(file, "#   theta_x    theta_y       z     weight        kappa          g_1          g_2\n");
      fprintf(file, "#  [arcmin]   [arcmin]      [-]       [-]          [-]          [-]          [-]\n");
    }
  }
  
  else {
    if (gMap->type < 12) {
      fprintf(file, "#       RA        DEC        z     weight        kappa      gamma_1      gamma_2\n");
      fprintf(file, "#     [deg]      [deg]      [-]       [-]          [-]          [-]          [-]\n");
    }
    else {
      fprintf(file, "#       RA        DEC        z     weight        kappa          g_1          g_2\n");
      fprintf(file, "#     [deg]      [deg]      [-]       [-]          [-]          [-]          [-]\n");
    }
  }
  
  if (pkPar->field == 0) {
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	g = gNode->g;
	fprintf(file, "  %9.3f  %9.3f  %7.5f  %8.5f  %11.4e  %11.4e  %11.4e\n", g->pos[0]*RADIAN_TO_ARCMIN, g->pos[1]*RADIAN_TO_ARCMIN, g->z, g->weight, g->kappa, g->gamma[0], g->gamma[1]);
      }
    }
  }
  
  else {
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	g = gNode->g;
	fprintf(file, "  %9.4f  %9.4f  %7.5f  %8.5f  %11.4e  %11.4e  %11.4e\n", g->pos[0]*RADIAN_TO_DEGREE, g->pos[1]*RADIAN_TO_DEGREE, g->z, g->weight, g->kappa, g->gamma[0], g->gamma[1]);
      }
    }
  }
  return;
}

#ifdef __CAMELUS_USE_FITS__
void outFits_gal_t(FITS_t *fits, gal_t *g, double factor)
{
  float fBuff;
  fBuff = (float)(g->pos[0] * factor); writeTableColumn(fits, 0, 1, &fBuff);
  fBuff = (float)(g->pos[1] * factor); writeTableColumn(fits, 1, 1, &fBuff);
  fBuff = (float)g->z;                 writeTableColumn(fits, 2, 1, &fBuff);
  fBuff = (float)g->weight;            writeTableColumn(fits, 3, 1, &fBuff);
  fBuff = (float)g->kappa;             writeTableColumn(fits, 4, 1, &fBuff);
  fBuff = (float)g->gamma[0];          writeTableColumn(fits, 5, 1, &fBuff);
  fBuff = (float)g->gamma[1];          writeTableColumn(fits, 6, 1, &fBuff);
  nextRow(fits);
  return;
}

void outFits_gal_map(FITS_t *fits, peak_param *pkPar, gal_map *gMap)
{
  if (pkPar->field == 0) {
    addColumn(fits, "THX", TFLOAT, "arcmin");
    addColumn(fits, "THY", TFLOAT, "arcmin");
  }
  else {
    addColumn(fits, "RA",  TFLOAT, "deg");
    addColumn(fits, "DEC", TFLOAT, "deg");
  }
  addColumn(fits, "Z",      TFLOAT, "-       ");
  addColumn(fits, "WEIGHT", TFLOAT, "-       ");
  addColumn(fits, "KAPPA",  TFLOAT, "-       ");
  if (gMap->type < 12) {
    addColumn(fits, "GAMMA1", TFLOAT, "-       ");
    addColumn(fits, "GAMMA2", TFLOAT, "-       ");
  }
  else {
    addColumn(fits, "G1",     TFLOAT, "-       ");
    addColumn(fits, "G2",     TFLOAT, "-       ");
  }
  
  addLineSpread(fits);
  addKeyword(fits, TINT, "NBGAL", &gMap->total, "[-] Number of galaxies");
  
  double factor = (pkPar->field == 0) ? RADIAN_TO_ARCMIN : RADIAN_TO_DEGREE;
  gal_list *gList;
  gal_node *gNode;
  int i, j;
  
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) outFits_gal_t(fits, gNode->g, factor);
  }
  return;
}
#endif

void read_gal_map(char name[], cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err)
{
  reset_gal_map(gMap);
  
  //-- Open
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  if (pkPar->verbose < 2) {
    printf("Reading...\r"); fflush(stdout);
  }
  
  //-- Column indices, inGalCatCol successively read as pos1, pos2, z, weight, kappa, e_1, e_2
  int pos1Ind  = pkPar->inGalCatCol[0];
  int pos2Ind  = pkPar->inGalCatCol[1];
  int zInd     = pkPar->inGalCatCol[2];
  int weiInd   = pkPar->inGalCatCol[3];
  int kappaInd = pkPar->inGalCatCol[4];
  int e1Ind    = pkPar->inGalCatCol[5];
  int e2Ind    = pkPar->inGalCatCol[6];
  
  int count  = 0;
  int count2 = 0;
  double factor = (pkPar->field == 0) ? ARCMIN_TO_RADIAN : DEGREE_TO_RADIAN;
  
  char kv[STRING_LENGTH_MAX][256], buffer[STRING_LENGTH_MAX];
  double pos[2], gamma[2], z, weight, kappa;
  
  //-- Read
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      ungetc(c, file);
      fgets(buffer, STRING_LENGTH_MAX, file);
      getKeyAndValues(buffer, kv);
      count++;
      
      pos[0]   = atof(kv[pos1Ind]) * factor; //-- [rad]
      pos[1]   = atof(kv[pos2Ind]) * factor; //-- [rad]
      z        = atof(kv[zInd]);
      weight   = (weiInd < 0)   ? ((pkPar->doNoise == 0) ? 1.0 : 1.0 / (SQ(pkPar->sigma_half))) : atof(kv[weiInd]); //-- Constant weight
      
      if (kappaInd < 0 && e1Ind < 0) {
	count2 += append_gal_map(chPar, pkPar, gMap, pos, z, weight, err);
	forwardError(*err, __LINE__,);
      }
      else {
	kappa    = (kappaInd < 0) ? 0.0 : atof(kv[kappaInd]);
	gamma[0] = (e1Ind < 0)    ? 0.0 : atof(kv[e1Ind]);
	gamma[1] = (e2Ind < 0)    ? 0.0 : atof(kv[e2Ind]);
	count2  += appendWithLensing_gal_map(pkPar, gMap, pos, z, weight, kappa, gamma, err); forwardError(*err, __LINE__,);
      }
    }
    c = fgetc(file);
  }
  fclose(file);
  
  //-- Check
  testErrorRet(count2!=gMap->total, peak_match, "Number not matched after reading files", *err, __LINE__,);
  if (pkPar->verbose < 3) {
    printf("Read \"%s\"   \n", name);
    printf("Found %d galaxies, generated %d\n", count, count2);
  }
  else if (pkPar->verbose == 3) {
    printf("read %d galaxies, ", count2); fflush(stdout);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to making galaxies

void makeRegularGalaxies(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err)
{
  //-- Regular implies necessairily fixed redshift.
  
  reset_gal_map(gMap);
  
  double *Omega = pkPar->Omega;
  int M1 = (int)round(Omega[0] * sqrt(pkPar->n_gal));
  int M2 = (int)round(Omega[1] * sqrt(pkPar->n_gal));
  long nsidePix = pkPar->nsidePix; //-- 16384
  int HP_length = pkPar->HP_length;
  long HP_first = pkPar->HP_first;
  double weight = (pkPar->doNoise == 0) ? 1.0 : 1.0 / (SQ(pkPar->sigma_half));
  
  double pos[2];
  //hpint64 pixNest;
  int i, j;

  if (pkPar->field == 0) {
     for (j=0; j<M2; j++) {
        pos[1] = (j + 0.5) / M2 * Omega[1];
        for (i=0; i<M1; i++) {
           pos[0] = (i + 0.5) / M1 * Omega[0];
           append_gal_map(chPar, pkPar, gMap, pos, pkPar->z_s, weight, err); forwardError(*err, __LINE__,);
        }
     }
  }

  else {
#ifdef __CAMELUS_USE_HEALPIX__
     for (i=0; i<HP_length; i++) {
        pix2ang_nest64((hpint64)nsidePix, HP_first+i, &pos[1], &pos[0]);
        pos[1] = HALF_PI - pos[1];
        append_gal_map(chPar, pkPar, gMap, pos, pkPar->z_s, weight, err); forwardError(*err, __LINE__,);
     }
#endif
  }
  return;
}

void setGalaxySampler(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, error **err)
{
  if (pkPar->doInGalCat == 1) return;
  
  gSamp->dx   = pkPar->dz_gal;
  gSamp->x[0] = chPar->redshift->par_nz[0];
  int i;
  for (i=1; i<gSamp->length; i++) gSamp->x[i] = gSamp->x[i-1] + gSamp->dx;
  
  //-- Structure for redshift law
  redshiftANDint *rANDi = (redshiftANDint*)malloc_err(sizeof(redshiftANDint), err); forwardError(*err, __LINE__,);
  rANDi->self = chPar->redshift;
  rANDi->i    = 0;
  fillGalaxyLaw(rANDi, gSamp, err);                                                 forwardError(*err, __LINE__,);
  free(rANDi);
  return;
}

void fillGalaxyLaw(redshiftANDint *rANDi, sampler_t *gSamp, error **err)
{
  double *z   = gSamp->x;
  double *pdf = gSamp->pdf;
  int i;
  
  //-- Fill pdf
  for (i=0; i<gSamp->length; i++) {
    pdf[i] = prob_unnorm(z[i], (void*)rANDi, err);
    forwardError(*err, __LINE__,);
  }
  
  int setTotalToOne = 1;
  set_sampler_t(gSamp, setTotalToOne);
  return;
}

void makeRandomGalaxies(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, gal_map *gMap, error **err)
{
  reset_gal_map(gMap);
  
  gsl_rng *generator = pkPar->generator;
  int N_gal     = (int)round(pkPar->area * pkPar->n_gal);
  double z_s    = pkPar->z_s;
  double weight = (pkPar->doNoise == 0) ? 1.0 : 1.0 / (SQ(pkPar->sigma_half));
  double p, z, pos[2];
  int i;
  
  if (z_s > 0) {
    for (i=0; i<N_gal; i++) {
      randomizePosition(pkPar, pos, err);                        forwardError(*err, __LINE__,);
      append_gal_map(chPar, pkPar, gMap, pos, z_s, weight, err); forwardError(*err, __LINE__,);
    }
  }
  else {
    //-- Loop for generating galaxies
    for (i=0; i<N_gal; i++) { //-- Should change if chPar->redshift->Nzbin > 1
      p = gsl_ran_flat(generator, 0.0, 1.0);
      z = execute_sampler_t(gSamp, p);
      randomizePosition(pkPar, pos, err);                        forwardError(*err, __LINE__,);
      append_gal_map(chPar, pkPar, gMap, pos, z, weight, err);   forwardError(*err, __LINE__,);
    }
  }
  
  if      (pkPar->verbose < 3)   printf("Generated %d galaxies, no mask\n", N_gal);
  else if (pkPar->verbose == 3) {printf("%d galaxies, ", N_gal); fflush(stdout);}
  return;
}

void makeMaskedRandomGalaxies(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, gal_map *gMap, mask_map *mask, error **err)
{
  reset_gal_map(gMap);
  
  gsl_rng *generator = pkPar->generator;
  int N_gal     = (int)round(pkPar->area * pkPar->n_gal);
  double z_s    = pkPar->z_s;
  double weight = (pkPar->doNoise == 0) ? 1.0 : 1.0 / (SQ(pkPar->sigma_half));
  double p, z, pos[2];
  int i;
  
  if (pkPar->z_s > 0) {
    for (i=0; i<N_gal; i++) {
      randomizePosition(pkPar, pos, err);                        forwardError(*err, __LINE__,);
      if (inMask(mask, pos)) continue;
      append_gal_map(chPar, pkPar, gMap, pos, z_s, weight, err); forwardError(*err, __LINE__,);
    }
  }
  else {
    //-- Loop for generating galaxies
    for (i=0; i<N_gal; i++) { //-- Should change if chPar->redshift->Nzbin > 1
      p = gsl_ran_flat(generator, 0.0, 1.0);
      z = execute_sampler_t(gSamp, p);
      randomizePosition(pkPar, pos, err);                        forwardError(*err, __LINE__,);
      if (inMask(mask, pos)) continue;
      append_gal_map(chPar, pkPar, gMap, pos, z, weight, err);   forwardError(*err, __LINE__,);
    }
  }
  
  if      (pkPar->verbose < 3)   printf("%d galaxies not masked, %d in total\n", gMap->total, N_gal);
  else if (pkPar->verbose == 3) {printf("%d galaxies, ", gMap->total); fflush(stdout);}
  return;
}

void rebinGalaxiesForHEALPix(peak_param *pkPar, gal_map *gMap1, gal_map *gMap2, error **err)
{
#ifdef __CAMELUS_USE_HEALPIX__
  reset_gal_map(gMap2);
  
  long HP_resol = pkPar->HP_resol;
  long nsidePix = pkPar->nsidePix;
  long HP_first = pkPar->HP_first;
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  double theta, phi;
  long pixNest;
  int i, j, carte;
  
  for (i=0; i<gMap1->length; i++) {
    gList = gMap1->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      g = gNode->g;
      theta = HALF_PI - g->pos[1];
      phi   = g->pos[0];
      ang2pix_nest(nsidePix, theta, phi, &pixNest);
      localNestToCartesian(pkPar->HP_resol, (int)(pixNest-HP_first), &carte);
      
      appendWithLensing_gal_list(pkPar, gMap2->map[carte], g->pos, g->z, g->weight, g->kappa, g->gamma, err); forwardError(*err, __LINE__,);
      gMap2->total++;
    }
  }
#endif
  return;
}

void cleanOrMakeOrResample(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, gal_map *gMap, mask_map *mask, error **err)
{
  if (pkPar->doInGalCat == 0) {
    if (pkPar->doRandGalPos == 1) {
      if (pkPar->doMask == 0) {
	makeRandomGalaxies(chPar, pkPar, gSamp, gMap, err);
	forwardError(*err, __LINE__,);
      }
      else if (pkPar->doMask == 1) {
	makeRandomMask(pkPar, mask, err);                               forwardError(*err, __LINE__,);
	makeMaskedRandomGalaxies(chPar, pkPar, gSamp, gMap, mask, err); forwardError(*err, __LINE__,);
      }
      else {
	makeMaskedRandomGalaxies(chPar, pkPar, gSamp, gMap, mask, err);
	forwardError(*err, __LINE__,);
      }
    }
    else if (gMap->total > 0) cleanLensing_gal_map(gMap);
    else {
      makeRegularGalaxies(chPar, pkPar, gMap, err);
      forwardError(*err, __LINE__,);
    }
  }
  else {
    read_gal_map(pkPar->inGalCatPath, chPar, pkPar, gMap, err);
    forwardError(*err, __LINE__,);
  }
  return;
}

//----------------------------------------------------------------------

