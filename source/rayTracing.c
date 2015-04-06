

  /***********************************************************************
   **  rayTracing.c							**
   **  Chieh-An Lin, Martin Kilbinger, François Lanusse			**
   **  Version 2015.04.06						**
   **									**
   **  References:							**
   **  - Bartelmann & Schneider (2001) - Phys. Rep., 340, 291 (BS01)	**
   **  - Takada & Jain (2003a) - MNRAS, 340, 580 (TJ03a)		**
   ***********************************************************************/


#include "rayTracing.h"


//----------------------------------------------------------------------
//-- Functions related to gal_t, gal_node, gal_list

void set_gal_t(cosmo_hm *cmhm, gal_t *g, double z, double w_s, double D_s, double pos[2], error **err)
{
  testError(g==NULL, peak_null, "gal_t *g is NULL", *err, __LINE__); forwardError(*err, __LINE__,);
  g->z     = z;
  double a = 1.0/(1.0+z);
  g->a     = a;
  g->kappa = 0.0;
  
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
  gList->size = 0;
  return;
}

void append_gal_list(cosmo_hm *cmhm, gal_list *gList, double z, double w_s, double D_s, double pos[2], error **err)
{
  if (gList->length == 0) {
    gList->first = initialize_gal_node(err);         forwardError(*err, __LINE__,);
    gList->last  = gList->first;
    gList->length++;
  }
  else if (gList->length == gList->size) {
    gList->last->next = initialize_gal_node(err);    forwardError(*err, __LINE__,);
    gList->last       = gList->last->next;
    gList->length++;
  }
  else if (gList->size == 0) gList->last = gList->first;
  else                       gList->last = gList->last->next;
  
  set_gal_t(cmhm, gList->last->g, z, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
  gList->size++;
  return;
}

double kappaMean_gal_list(gal_list *gList, error **err)
{
  int size = gList->size;
  if (size == 0) return 0.0;
  gal_node *gNode;
  gal_t *g;
  double mean = 0.0;
  int i;
  for (i=0, gNode=gList->first; i<size; i++, gNode=gNode->next) mean += gNode->g->kappa;
  return mean / (double)size;
}

//----------------------------------------------------------------------
//-- Functions related to gal_map

gal_map *initialize_gal_map(int N1, int N2, double theta_pix, error **err)
{
  gal_map *gMap   = (gal_map*)malloc_err(sizeof(gal_map), err);                    forwardError(*err, __LINE__,);
  gMap->N1        = N1;
  gMap->N2        = N2;
  gMap->length    = N1 * N2;
  gMap->total     = 0;
  gMap->theta_pix = theta_pix;
  
  gMap->center[0] = 0.0;
  gMap->center[1] = 0.0;
  gMap->limits[0] = -0.5 * N1 * theta_pix;
  gMap->limits[1] = -gMap->limits[0];
  gMap->limits[2] = -0.5 * N2 * theta_pix;
  gMap->limits[3] = -gMap->limits[2];
  
  gMap->map       = (gal_list**)malloc_err(gMap->length * sizeof(gal_list*), err); forwardError(*err, __LINE__,);
  int i;
  for (i=0; i<gMap->length; i++) gMap->map[i] = initialize_gal_list(err);
  return gMap;
}

void free_gal_map(gal_map *gMap)
{
  int i;
  for (i=0; i<gMap->length; i++) {free_gal_list(gMap->map[i]); gMap->map[i] = NULL;}
  free(gMap); gMap = NULL;
  return;
}

void cleanLensing_gal_map(gal_map *gMap)
{
  int i;
  for (i=0; i<gMap->length; i++) cleanLensing_gal_list(gMap->map[i]);
  return;
}

void reset_gal_map(gal_map *gMap)
{
  gMap->total = 0;
  int i;
  for (i=0; i<gMap->length; i++) reset_gal_list(gMap->map[i]);
  return;
}

void append_gal_map(cosmo_hm *cmhm, gal_map *gMap, double z, double w_s, double D_s, double pos[2], error **err)
{
  double theta_pix = gMap->theta_pix;
  int i = (int)((pos[0] - gMap->limits[0]) / theta_pix);
  int j = (int)((pos[1] - gMap->limits[2]) / theta_pix);
  append_gal_list(cmhm, gMap->map[i+j*gMap->N1], z, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
  gMap->total++;
  return;
}

void output_gal_map(FILE *file, peak_param *peak, gal_map *gMap)
{
  fprintf(file, "# Number of galaxies = %d\n", gMap->total);
  fprintf(file, "#\n");
  fprintf(file, "#  theta_x    theta_y        z     kappa\n");
  fprintf(file, "# [arcmin]   [arcmin]      [-]       [-]\n");
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  int i, j;
  for (i=0; i<gMap->length; i++) {
    gList = gMap->map[i];
    for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
      g = gNode->g;
      fprintf(file, " %9.3f  %9.3f  %7.5f  %8.5f\n", g->pos[0], g->pos[1], g->z, g->kappa);
    }
  }
  return;
}

void updateCosmo_gal_map(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  double w_s, D_s;
  if (peak->z_s > 0) {
    int wOmegar = 0;
    w_s  = w(cmhm->cosmo, 1.0/(1.0 + peak->z_s), wOmegar, err); forwardError(*err, __LINE__,);
    D_s  = f_K(cmhm->cosmo, w_s, err);                           forwardError(*err, __LINE__,);
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
      set_gal_t(cmhm, g, g->z, w_s, D_s, NULL, err);
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to projected NFW mass

#define EPS 1e-10
double G_NFW_kappa(double x_sq, double c, double c_sq, error **err)
{
  //-- See TJ03a (27)
  //-- G(x) = -arg1/arg3 + arg3^(-3/2) * arcosh(arg2), if x < 1;
  //--        1/3 * arg1 * (c + 2) / (c + 1),          if x = 1;
  //--        +arg1/arg3 - arg3^(-3/2) * arccos(arg2), if 1 < x <= c;
  //--        0,                                       if x > c;
  //-- where
  //--   arg1 = sqrt(c^2 - x^2) / (c + 1),
  //--   arg2 = (x^2 + c) / (x(c + 1)),
  //--   arg3 = |1 - x^2|.
  if (x_sq > c_sq - EPS) return 0;
  if (x_sq < EPS) x_sq = EPS;
  
  double p    = c + 1;
  double arg1 = sqrt(c_sq - x_sq) / p;
  double arg3 = fabs(x_sq - 1);
  
  double value;
  if (arg3 < EPS) value = arg1 / 3.0 * (p + 1) / p;
  else {
    double arg2 = (x_sq + c) / (sqrt(x_sq) * p);
    if (x_sq < 1) {
      testErrorRet(arg2<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__, 0.0); forwardError(*err, __LINE__,);
      value  = -arg1 + acosh(arg2) / sqrt(arg3);
      value /= arg3;
    }
    else {
      testErrorRet(arg2>1.0, peak_badValue, "Out of range for acos", *err, __LINE__, 0.0); forwardError(*err, __LINE__,);
      value  = arg1 - acos(arg2) / sqrt(arg3);
      value /= arg3;
    }
  }
  
  return value;
}

double G_NFW_gamma(double x_sq, double c, double c_sq, double f, error **err)
{
  //-- See TJ03b (17)
  //-- G(x) = x^-2 [arg1 + 2 ln(arg2) + arg3 * arccosh(arg4)], if x < 1;
  //--        x^-2 [arg5 + 2 ln(arg2)],                        if x = 1;
  //--        x^-2 [arg1 + 2 ln(arg2) + arg3 * arccos(arg4)],  if 1 < x <= c;
  //--        x^-2 * 2/f,                                      if x > c;
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
  
  double value;
  if (x_sq < EPS)        x_sq = EPS;
  if (c_sq - x_sq < EPS) value = 2.0 / f;
  else {
    double p    = c + 1;
    double q    = 1 - x_sq;
    double r    = sqrt(c_sq - x_sq);
    double x    = sqrt(x_sq);
    double arg2 = x * p / (c + r);
    
    if (fabs(q) < EPS) value = ((c + 10*p) * r - 6 * c * p) / (3 * SQ(p)) + 2 * log(arg2);
    else {
      double arg1 = ((q + 1) * r - 2 * c * q) / (p * q);
      double arg3 = (3 * q - 1) / (q * sqrt(fabs(q)));
      double arg4 = (x_sq + c) / (x * p);
      if (q > 0) {
	testErrorRet(arg4<1.0, peak_badValue, "Out of range for acosh", *err, __LINE__, 0.0); forwardError(*err, __LINE__,);
	value = arg1 + 2 * log(arg2) + arg3 * acosh(arg4);
      }
      else {
	testErrorRet(arg4>1.0, peak_badValue, "Out of range for acos", *err, __LINE__, 0.0);  forwardError(*err, __LINE__,);
	value = arg1 + 2 * log(arg2) + arg3 * acos(arg4);
      }
    }
  }
  value /= x_sq;
  return value;
}
#undef EPS

double kappa_NFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- See BS01 (3.7)
  //-- kappa_NFW = Sigma / Sigma_crit
  //-- where
  //--   Sigma_crit = (c^2 / 4 pi G) * (D_s / D_l D_ls), where (4 pi G / c^2) = 3 / (2 CRITICAL_DENSITY D_H^2 ) = 6.01e-19 [Mpc/M_sol]
  //--   Sigma = projected mass
  //--         = 2 rho_s r_s * G_NFW_kappa(theta / theta_s)
  //--         = 2 * (rho_bar * Delta f c_NFW^3 / 3) * (r_vir / c_NFW) * G_NFW_kappa(c_NFW theta / theta_vir)
  //--         = 2 * (rho_bar * Delta * 4 pi r_vir^3 / 3) * (f c_NFW^3 / 4 pi r_vir^3) * (r_vir / c_NFW) * G_NFW_kappa(c_NFW theta / theta_vir)
  //--         = M * (f c_NFW^2 / 2 pi r_vir^2) * G_NFW_kappa(c_NFW theta / theta_vir)
  //-- 
  //-- Similarly, in TJ03a (26), (28), one finds
  //-- kappa_NFW = W(w_l, w_s) * (M / rho_bar) * Sigma_m(theta)
  //--           = (3 Omega_m H_0^2 / 2) * (D_l D_ls / D_s a_l) * (M / rho_bar) * Sigma_m
  //--           = (4 pi G) * (D_l D_ls / D_s a_l) * M * (f c_NFW^2 / 2 pi R_vir^2) * G_NFW_kappa(c_NFW theta / theta_vir)
  //-- However, their convention is different. First, they have omitted c^(-2).
  //-- Second, their D_l, D_ls, D_s are actually f_K: f_l, f_ls, f_s.
  //-- Finally, their R_vir is defined as a comoving distance, not the physical size of the halo, so it should be the real r_vir divided by a_l.
  //-- Therefore, the correct equation for TJ03a is
  //-- kappa_NFW = (4 pi G / c^2) * (f_l f_ls / f_s a_l) * M * (f c_NFW^2 a_l^2 / 2 pi r_vir^2) * G_NFW_kappa(c_NFW theta / theta_vir)
  //-- 
  //-- Since (f_l f_ls) / (f_s a_l) = (D_l D_ls) / (D_s a_l^2), one finds
  //-- kappa_NFW = (4 pi G / c^2) * (D_l D_ls / D_s a_l^2) * M * (f c_NFW^2 a_l^2 / 2 pi r_vir^2) * G_NFW_kappa(c_NFW theta / theta_vir)
  //--           = (4 pi G / c^2) * (D_l D_ls / D_s) * M * (f c_NFW^2 / 2 pi r_vir^2) * G_NFW_kappa(c_NFW theta / theta_vir)
  //--           = Sigma / Sigma_crit
  //-- Notice also that their Sigma_m is dimensionless, so that Sigma = M * Sigma_m.
  //-- 
  //-- For computation purposes, we write
  //-- kappa_NFW = (4 pi G / c^2) * (D_l D_ls M f c_NFW^2 / 2 pi D_s r_vir^2) * G_NFW_kappa(c_NFW theta / theta_vir)
  //--           = (4 pi G / c^2) * (D_l M f c_NFW^2 / 2 pi r_vir^2) * (D_ls / D_s) * G_NFW_kappa(c_NFW theta / theta_vir)
  //--           = factor * (D_ls / D_s) * G_NFW_kappa(c_NFW theta / theta_vir)
  //-- where factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  
  double x_sq  = h->c_sq * theta_sq / h->theta_vir_sq; //-- x^2 = theta^2 / theta_s^2
  double D_ls  = g->a * f_K(cmhm->cosmo, g->w-h->w, err);                           forwardError(*err, __LINE__, 0.0);
  double kappa = h->factor * D_ls / g->D_s * G_NFW_kappa(x_sq, h->c, h->c_sq, err); forwardError(*err, __LINE__, 0.0);
  return kappa;
}

void gamma_NFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err)
{
  //-- See TJ03b (16, 18)
  //-- |gamma_NFW| = R / Sigma_crit
  //-- where
  //--   Sigma_crit = (c^2 / 4 pi G) * (D_s / D_l D_ls), where (4 pi G / c^2) = 3 / (2 CRITICAL_DENSITY D_H^2 ) = 6.01e-19 [Mpc/M_sol]
  //--   R = M * (f c_NFW^2 / 2 pi r_vir^2) * G_NFW_gamma(c_NFW theta / theta_vir)
  //-- 
  //-- For computation purposes, we write
  //-- |gamma_NFW| = (4 pi G / c^2) * (D_l D_ls M f c_NFW^2 / 2 pi D_s r_vir^2) * G_NFW_gamma(c_NFW theta / theta_vir)
  //--             = factor * (D_ls / D_s) * G_NFW_kappa(c_NFW theta / theta_vir)
  //-- where factor = FOUR_PI_G_OVER_C2 * (D_l M f c_NFW^2 / 2 pi r_vir^2)
  //--
  //-- gamma1  = -|gamma| cos(2 phi) = -|gamma| (1 - t^2) / (1 + t^2)
  //-- gamma2  = -|gamma| sin(2 phi) = -|gamma|       2t  / (1 + t^2)
  //-- where t = tan(phi)
  
  double x_sq       = h->c_sq * theta_sq / h->theta_vir_sq; //-- x^2 = theta^2 / theta_s^2
  double D_ls       = g->a * f_K(cmhm->cosmo, g->w-h->w, err);                                 forwardError(*err, __LINE__,);
  double gamma_norm = h->factor * D_ls / g->D_s * G_NFW_gamma(x_sq, h->c, h->c_sq, h->f, err); forwardError(*err, __LINE__,);
  
  double t    = tan(atan2(g->pos[1]-h->pos[1], g->pos[0]-h->pos[0]));
  double t_sq = SQ(t);
  double tcos = (1.0 - t_sq) / (1.0 + t_sq);
  double tsin =       2 * t  / (1.0 + t_sq);
  g->gamma[0] += -gamma_norm * tcos;
  g->gamma[1] += -gamma_norm * tsin;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to lensing

void lensingForPair(cosmo_hm *cmhm, halo_t *h, gal_t *g, int doKappa, error **err)
{
  //-- Lensing for a halo-galaxy pair
  testError(g==NULL, peak_null, "Empty galaxy", *err, __LINE__); forwardError(*err, __LINE__,);
  
  if (g->z < h->z) return; //-- Source in front of lens
  double theta_sq = DIST_2D_SQ(h->pos, g->pos);
  if (doKappa) {
    //-- Compute projected kappa
    if (theta_sq > h->theta_vir_sq) return;
    g->kappa += kappa_NFW(cmhm, h, g, theta_sq, err); forwardError(*err, __LINE__,);
  }
  else {
    //-- Compute projected gamma
    if (theta_sq > CUTOFF_FACTOR_HALO * h->theta_vir_sq) return;
    gamma_NFW(cmhm, h, g, theta_sq, err);             forwardError(*err, __LINE__,);
  }
  return;
}

void lensingForHalo(cosmo_hm *cmhm, halo_t *h, gal_map *gMap, int doKappa, error **err)
{
  //-- Lensing for a halo-galaxy pair
  testError(h==NULL, peak_null, "Empty halo", *err, __LINE__); forwardError(*err, __LINE__,);
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  double theta_pix = gMap->theta_pix;
  double theta_x   = h->pos[0] - gMap->limits[0];
  double theta_y   = h->pos[1] - gMap->limits[2];
  double cutoff    = h->theta_vir;
  if (!doKappa) cutoff *= CUTOFF_FACTOR_HALO;
  
  int i_min = (int)((theta_x - cutoff) / theta_pix);
  int i_max = (int)((theta_x + cutoff) / theta_pix);
  int j_min = (int)((theta_y - cutoff) / theta_pix);
  int j_max = (int)((theta_y + cutoff) / theta_pix);
  
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  
  int i, j, k;
  for (j=j_min; j<=j_max; j++) {
    if (j < 0)   continue;
    if (j >= N2) continue;
    for (i=i_min; i<=i_max; i++) {
      if (i < 0)   continue;
      if (i >= N1) continue;
      gList = gMap->map[i+j*N1];
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
  cleanLensing_gal_map(gMap);
  
  halo_list *hList;
  halo_node *hNode;
  halo_t *h;
  
  int i, j, count = 0;
  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
    for (j=0, hNode=hList->first; j<hList->size; j++, hNode=hNode->next) {
      if ((peak->printMode == 0) && (count % 50 == 0)) {
	printf("Doing lensing signal computation: %6.2f%% \r", 100.0 * count / (double)hMap->total);
	fflush(stdout);
      }
      lensingForHalo(cmhm, hNode->h, gMap, peak->doKappa, err); forwardError(*err, __LINE__,);
      count++;
    }
  }
  
  if (peak->printMode != 1) printf("Lensing signal computation done                \n");
  return;
}

//----------------------------------------------------------------------
//-- Functions related to making galaxy lists

void makeGalaxiesFixedRedshift(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  double *Omega = peak->Omega;
  int N_gal  = (int)round(peak->area * peak->n_gal);
  int M1     = (int)round(Omega[0] * sqrt(peak->n_gal));
  int M2     = (int)round(Omega[1] * sqrt(peak->n_gal));
  double z_s = peak->z_s;
  double w_s = peak->w_s;
  double D_s = peak->D_s;
  double pos[2];
  
  int i, j;
  if (peak->doRandGalPos) {
    //-- Generate random galaxies
    for (i=0; i<N_gal; i++) {
      randomizePosition(peak, pos, err);                   forwardError(*err, __LINE__,);
      append_gal_map(cmhm, gMap, z_s, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
    }
    printf("%d galaxies generated, redshift fixed at %.3f, randomly distributed\n", gMap->total, peak->z_s);
  }
  else {
    //-- Generate regular galaxies
    for (j=0; j<M2; j++) {
      for (i=0; i<M1; i++) {
	pos[0] = ((i + 0.5) / M1 - 0.5) * Omega[0];
	pos[1] = ((j + 0.5) / M2 - 0.5) * Omega[1];
	append_gal_map(cmhm, gMap, z_s, w_s, D_s, pos, err); forwardError(*err, __LINE__,);
      }
    }
    printf("%d galaxies generated, redshift fixed at %.3f, on a regular grid (%dx%d)\n", gMap->total, peak->z_s, M1, M2);
  }
  return;
}

void fillGalaxyLaw(redshiftANDint *rANDi, sampler_t *samp, error **err)
{
  double *z   = samp->x;
  double *pdf = samp->pdf;
  double *cdf = samp->cdf;
  int i;
  
  //-- Fill pdf
  for (i=0; i<samp->length; i++) {
    pdf[i] = prob_unnorm(z[i], (void*)rANDi, err);
    forwardError(*err, __LINE__,);
  }
  
  int setTotalToOne = 1;
  set_sampler_t(samp, setTotalToOne);
  return;
}

void sampleGalaxies(redshiftANDint *rANDi, cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  //-- Redshift configuration
  double z_min = get_zmin(rANDi->self, 0);
  double z_max = get_zmax(rANDi->self, 0);
  double dz    = 0.001;                           //-- Bin width for redshift
  int N_z_gal  = (int)ceil((z_max - z_min) / dz); //-- Nb of redshift bins
  
  //-- Initialize sampler
  sampler_t *samp = initialize_sampler_t(N_z_gal);
  samp->dx        = dz;
  samp->x[0]      = z_min;
  int i;
  for (i=1; i<N_z_gal; i++) samp->x[i] = samp->x[i-1] + dz; //-- Fill bins
  
  //-- Fill samp
  fillGalaxyLaw(rANDi, samp, err); forwardError(*err, __LINE__,);
  
  gsl_rng *generator = peak->generator;
  int N_gal          = (int)round(peak->area * peak->n_gal);
  double p, z, pos[2];
  
  //-- Loop for generating galaxies
  gal_t *g;
  for (i=0; i<N_gal; i++) { //-- Should change if cmhm->redshift->Nzbin > 1
    p = gsl_ran_flat(generator, 0.0, 1.0);
    z = execute_sampler_t(samp, p);
    randomizePosition(peak, pos, err);                   forwardError(*err, __LINE__,);
    append_gal_map(cmhm, gMap, z, -1.0, -1.0, pos, err); forwardError(*err, __LINE__,);
  }
  free_sampler_t(samp);
  return;
}

void makeRealGalaxies(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  testError(cmhm->redshift->Nzbin!=1, peak_badValue, "cmhm->redshift->Nzbin should be 1", *err, __LINE__); forwardError(*err, __LINE__,);
  
  //-- Structure for redshift law
  redshiftANDint *rANDi = (redshiftANDint*)malloc_err(sizeof(redshiftANDint), err); forwardError(*err, __LINE__,);
  rANDi->self = cmhm->redshift;
  rANDi->i    = 0;
  
  sampleGalaxies(rANDi, cmhm, peak, gMap, err); forwardError(*err, __LINE__,); //-- Only for cmhm->redshift->Nzbin = 1 so far
  printf("%d galaxies generated, %s redshift law, randomly distributed\n", gMap->total, snofz_t(cmhm->redshift->nofz[0]));
  free(rANDi);
  return;
}

void makeGalaxies(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  reset_gal_map(gMap); //-- Reset
  if (peak->z_s > 0) {makeGalaxiesFixedRedshift(cmhm, peak, gMap, err); forwardError(*err, __LINE__,);}
  else               {makeRealGalaxies(cmhm, peak, gMap, err);          forwardError(*err, __LINE__,);}
  return;
}

void outputGalaxies(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap)
{
  FILE *file = fopen(name, "w");
  
  fprintf(file, "# Galaxy list\n");
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

void doRayTracing(char haloFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  halo_map *hMap = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  gal_map *gMap  = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  makeGalaxies(cmhm, peak, gMap, err);                                                        forwardError(*err, __LINE__,);
  
  if (haloFileName == NULL) {
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
    setMassSamplers(cmhm, peak, sampArr, err);       forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err);   forwardError(*err, __LINE__,);
    outputFastSimul("haloList", cmhm, peak, hMap);
    free_sampler_arr(sampArr);
  }
  else read_halo_map(haloFileName, cmhm, hMap, err); forwardError(*err, __LINE__,);
  
  lensingForMap(cmhm, peak, hMap, gMap, err);        forwardError(*err, __LINE__,);
  outputGalaxies("galList", cmhm, peak, gMap);
  
  free_halo_map(hMap);
  free_gal_map(gMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------


