

  /***************************************************************
   **  rayTracing.c						**
   **  Version 2018.03.14					**
   **								**
   **  References:						**
   **  - Baltz et al. (2009) - JCAP, 1, 15			**
   **  - Bartelmann & Schneider (2001) - Phys. Rep., 340, 291	**
   **  - Oguri & Hamana (2011) - MNRAS, 414, 1851		**
   **  - Takada & Jain (2003a) - MNRAS, 340, 580		**
   **  - Takada & Jain (2003b) - MNRAS, 344, 857		**
   **  - Wright & Brainerd (2000) - ApJ, 534, 34		**
   **								**
   **  Copyright (C) 2018 - Chieh-An Lin			**
   **  GNU GPLv3 - https://www.gnu.org/licenses/		**
   ***************************************************************/


#include "rayTracing.h"


//----------------------------------------------------------------------
//-- Functions related to gal_t, gal_node, gal_list

// Lin-tw: set_gal_t, initialize_gal_node, initialize_gal_list, and other functions now in galaxySampling
// Lin-tw: gammaOrGMean_gal_list in galaxySampling, changed (includes weight)

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
  g->w = w(cmhm->cosmo, g->a, 0, err);       forwardError(*err, __LINE__,); //-- wOmegar = 0
  g->D_s = g->a * f_K(cmhm->cosmo, g->w, err); forwardError(*err, __LINE__,);

  return;
}

// MKDBUG: a similar function is galaxySampling:appendWithLensing_gal_list
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

//----------------------------------------------------------------------
//-- Functions related to gal_map

// Lin-tw: similar function galaxySampling:appendWithLensing_gal_map
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

void output_gal_map(FILE *file, peak_param *peak, gal_map *gMap)
{
   fprintf(file, "# Type = %s, number of galaxies = %d\n", STR_MAP_T(gMap->type), gMap->total);
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

//-- Functions related to projected mass

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
  
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
  double arg1 = fabs(u_sq - 1.0);
  double value;
  
  if (arg1 < EPS_NUM) value = 1.0 / 3.0;
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
  
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
  double arg1 = fabs(u_sq - 1.0);
  double value;
  
  if (arg1 < EPS_NUM) {
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
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
  double arg1 = fabs(u_sq - 1.0);
  
  if (arg1 < EPS_NUM) {
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
  
  
  if (u_sq > c_sq - EPS_NUM) return 0;
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
  double p_inv = 1.0 / (c + 1.0);
  double arg1  = fabs(u_sq - 1.0);
  double arg2  = sqrt(c_sq - u_sq) * p_inv;
  double value;
  
  if (arg1 < EPS_NUM) value = arg2 * (1.0 + p_inv) / 3.0;
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
  
  if (u_sq > c_sq - EPS_NUM) return 2.0 / (f * u_sq);
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
  double p     = c + 1.0;
  double p_inv = 1.0 / p;
  double q     = 1.0 - u_sq;
  double r     = sqrt(c_sq - u_sq);
  double u     = sqrt(u_sq);
  double arg5  = u * p / (c + r);
  double value;
  
  if (fabs(q) < EPS_NUM) {
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
  if (u_sq > c_sq - EPS_NUM) {
    G_TJ[0] = 0.0;
    G_TJ[1] = 2.0 / (f * u_sq);
    return;
  }
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
  double p     = c + 1.0;
  double p_inv = 1.0 / p;
  double q     = 1.0 - u_sq;
  double r     = sqrt(c_sq - u_sq);
  double u     = sqrt(u_sq);
  double arg1  = fabs(q);
  double arg2  = r * p_inv;
  double arg5  = u * p / (c + r);
  double value;
  
  if (arg1 < EPS_NUM) {
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
  
  if (u_sq < EPS_NUM) u_sq = EPS_NUM;
  
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
  
  if (t < EPS_NUM) value += 2.0*p / 3.0 + 8.0;
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

double kappa_TJ(cosmo_hm *chPar, halo_t *h, gal_t *g, double theta_sq, error **err)
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
  double D_ls  = g->a * f_K(chPar->cosmo, g->w-h->w, err);                         forwardError(*err, __LINE__, -1.0);
  double kappa = h->factor * D_ls / g->D_s * G_TJ_kappa(u_sq, h->c, h->c_sq, err); forwardError(*err, __LINE__, -1.0);
  return kappa;
}

void gamma_TJ(cosmo_hm *chPar, halo_t *h, gal_t *g, double theta_sq, double phase, error **err)
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
  
  double u_sq    = h->c_sq * theta_sq / h->theta_vir_sq; //-- u^2 = theta^2 / theta_s^2
  double D_ls    = g->a * f_K(chPar->cosmo, g->w-h->w, err);                               forwardError(*err, __LINE__,);
  double gamma_t = h->factor * D_ls / g->D_s * G_TJ_gamma(u_sq, h->c, h->c_sq, h->f, err); forwardError(*err, __LINE__,);
  
  g->gamma[0] += gamma_t * cos(phase);
  g->gamma[1] += gamma_t * sin(phase);
  return;
}

void both_TJ(cosmo_hm *chPar, halo_t *h, gal_t *g, double theta_sq, double phase, error **err)
{
  double u_sq   = h->c_sq * theta_sq / h->theta_vir_sq; //-- u^2 = theta^2 / theta_s^2
  double G_TJ[2];
  G_TJ_both(u_sq, h->c, h->c_sq, h->f, G_TJ, err);          forwardError(*err, __LINE__,);
  
  double D_ls   = g->a * f_K(chPar->cosmo, g->w-h->w, err); forwardError(*err, __LINE__,);
  double factor = h->factor * D_ls / g->D_s;
  g->kappa     += factor * G_TJ[0]; //-- kappa_TJ
  G_TJ[1]      *= factor;           //-- gamma_t
  g->gamma[0]  += G_TJ[1] * cos(phase);
  g->gamma[1]  += G_TJ[1] * sin(phase);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to drawing a profile

void fillOneHaloTerm(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, gal_t *g, double_mat *profile, error **err)
{
  int verbose   = pkPar->verbose;
  int N_gal     = profile->N1;
  double D_ls   = g->a * f_K(chPar->cosmo, g->w-h->w, err); forwardError(*err, __LINE__,);
  double factor = h->factor * D_ls / g->D_s;
  double theta_sq, u_sq, tau; //G_halo[2];
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
    
    if (verbose < 2 && i % 125 == 0) {
      printf("Computing one-halo term: %6.2f%%\r", 100.0 * i / (double)N_gal);
      fflush(stdout);
    }
  }
  if (pkPar->verbose < 3) printf("Computed one-halo term          \n");
  return;
}

double factorForTwoHaloTerm(cosmo_hm *chPar, halo_t *h, gal_t *g, error **err)
{
  //-- Oguri & Hamana (2011), Eq. (13)
  //-- kappa_2h(theta) = \int dl [(l/2pi) * factor * J_0(l*theta) * P_L(k_l, z)]
  //-- where k_l = l / f_K(w(z))
  //-- where factor = (rho_bar(z) b_h(M)) / ((1+z)^3 Sigma_crit D_l^2)
  //--              = (rho_bar(z) / (1+z)^3) * b_h(M) / ((c^2 / 4 pi G) * (D_s / D_l D_ls) * D_l^2)
  //--              = rho_bar(0) * b_h(M) / ((c^2 / 4 pi G) * (D_s D_l / D_ls))
  //--              = (CRITICAL_DENSITY * Omega_m) * b_h(M) * FOUR_PI_G_OVER_C2 * D_ls / (D_s D_l)
  //--              = (FOUR_PI_G_OVER_C2 CRITICAL_DENSITY Omega_m b_h(M) D_ls) / (D_s D_l)
  
  double D_l    = h->a * f_K(chPar->cosmo, h->w, err);      forwardError(*err, __LINE__, -1.0);
  double D_ls   = g->a * f_K(chPar->cosmo, g->w-h->w, err); forwardError(*err, __LINE__, -1.0);
  double b_h    = halo_bias(chPar, h->M, h->a, 2, err);     forwardError(*err, __LINE__, -1.0); //-- order = 2 for bias_sc
  double factor = FOUR_PI_G_OVER_C2 * CRITICAL_DENSITY * chPar->cosmo->Omega_m * b_h * D_ls / (g->D_s * D_l);
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
  double k_l   = l / param->f_K;
  double value = l * bessj0(l*param->theta) * P_L(param->cosmo, param->a, k_l, err); forwardError(*err, __LINE__, -1.0);
  return value;
}

double twoHaloTerm(cosmo_hm *chPar, halo_t *h, double theta, double factor, error **err)
{
  profile_inteParam *inteParam = (profile_inteParam*)malloc_err(sizeof(profile_inteParam), err); forwardError(*err, __LINE__, -1.0);
  inteParam->theta = theta;
  inteParam->a     = h->a;
  inteParam->f_K   = f_K(chPar->cosmo, h->w, err); forwardError(*err, __LINE__, -1.0);
  inteParam->cosmo = chPar->cosmo;
  double l_maxx    = 500.0;
  double eps       = 1e-9;
  double totFactor = 0.5 * PI_INV * factor;
  double sum       = sm2_qromberg(integrandForTwoHaloTerm, (void*)inteParam, 0.0, l_maxx, eps, err) * totFactor;
  free(inteParam);
  return sum;
}

void fillTwoHaloTerm(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, double factor, double_mat *profile, error **err)
{
  int verbose = pkPar->verbose;
  int N_gal   = profile->N1;
  double theta;
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta = profile->matrix[i] * ARCMIN_TO_RADIAN; //-- Conversion
    profile->matrix[i+3*N_gal] = twoHaloTerm(chPar, h, theta, factor, err); forwardError(*err, __LINE__,); //-- 2-halo term kappa
    if (verbose < 2 && i % 25 == 0) {
      printf("Computing the two-halo term: %6.2f%%\r", 100.0 * i / (double)N_gal);
      fflush(stdout);
    }
  }
  if (pkPar->verbose < 3) printf("Computed two-halo term          \n");
  return;
}

void outAsciiProfile(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_t *h, double_mat *profile, error **err)
{
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Halo profile\n");
  fprintf(file, "#\n");
  fprintf(file, "# Key-in parameters\n");
  fprintf(file, "# z_l = %f, M = %8.2e, z_s = %f\n", h->z, h->M, pkPar->z_halo_max);
  fprintf(file, "#\n");
  fprintf(file, "# Halo model parameters = %s\n", pkPar->hmParPath);
  fprintf(file, "# c_0 = %f, beta_NFW = %f, bias = %s\n", chPar->c0, chPar->beta_NFW, shalo_bias_t(chPar->halo_bias));
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "#      theta       r(phys)        kappa        gamma           g      kappa_2h\n");
  fprintf(file, "#    [arcmin]      [Mpc/h]          [-]          [-]          [-]          [-]\n");
  
  int N_gal  = profile->N1;
  double D_l = h->a * f_K(chPar->cosmo, h->w, err); forwardError(*err, __LINE__,);
  double theta, kappa, gamma, kappa_2h;
  int i;
  
  for (i=0; i<N_gal; i++) {
    theta    = profile->matrix[i];
    kappa    = profile->matrix[i+N_gal];
    gamma    = profile->matrix[i+2*N_gal];
    kappa_2h = profile->matrix[i+3*N_gal];
    fprintf(file, "  %11.5e  %11.5e  %11.5e  %11.5e  %11.5e  %11.5e\n", theta, theta*ARCMIN_TO_RADIAN*D_l, kappa, gamma, gamma/(1.0-kappa), kappa_2h);
  }
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to lensing

void lensingForPair(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, gal_t *g, error **err)
{
  //-- Lensing for a halo-galaxy pair
  
  if (g->w <= h->w) return; //-- Source in front of lens
  
  double theta_sq = (pkPar->field == 0)   ? DIST_2D_SQ(h->pos, g->pos) : SQ(SPHE_DIST(h->pos, g->pos));
  double factor   = (pkPar->doKappa == 1) ? 1.0 : CUTOFF_FACTOR_HALO_SQ;
  if (theta_sq > factor * h->theta_vir_sq) return; //-- Too far
  
  double dtheta_x, dtheta_y, dRA, phase;
  
  //-- Compute projected kappa
  if (pkPar->doKappa == 1) {
    g->kappa += kappa_TJ(chPar, h, g, theta_sq, err);
    forwardError(*err, __LINE__,);
  }
  
  else {
    //-- gamma_t = gamma_1 * cos(phase) + gamma_2 * sin(phase)
    //-- gamma_1 = gamma_t * cos(phase)
    //-- gamma_2 = gamma_t * sin(phase)
    //-- phase = pi + 2 angle = pi + 2 arctan(y/x)
  
    if (pkPar->field == 0) {
      phase    = PI + 2.0 * atan2(g->pos[1]-h->pos[1], g->pos[0]-h->pos[0]);
    }
    else {
      dRA      = g->pos[0] - h->pos[0];
      dtheta_x = g->sinCosDEC[1] * sin(dRA);
      dtheta_y = (h->sinCosDEC[1] * g->sinCosDEC[0] - h->sinCosDEC[0] * g->sinCosDEC[1] * cos(dRA));
      phase    = PI + 2.0 * (atan2(dtheta_y, dtheta_x) + pkPar->rotAng);
    }
    
    //-- Compute projected kappa and gamma
    if (pkPar->doKappa == 0) {
      gamma_TJ(chPar, h, g, theta_sq, phase, err);
      forwardError(*err, __LINE__,);
    }
    
    //-- Compute projected gamma
    else {
      both_TJ(chPar, h, g, theta_sq, phase, err);
      forwardError(*err, __LINE__,);
    }
  }
  return;
}

void lensingForHalo(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, gal_map *gMap, int i_h, int j_h, error **err)
{
  //-- i_h, j_h are required by HEALPix.
  
  //-- Lensing for a halo-galaxy pair
  testErrorRet(h==NULL, peak_null, "Empty halo", *err, __LINE__,);
  
  int N1 = gMap->N1;
  int N2 = gMap->N2;
  double theta_pix_inv = 1.0 / gMap->theta_pix;
  double cutInPix = h->theta_vir * theta_pix_inv;
  
  if (pkPar->doKappa != 1) cutInPix *= CUTOFF_FACTOR_HALO;
  if (pkPar->field > 0)    cutInPix *= CUTOFF_FACTOR_HEALPIX; //-- Because of distance dilatation by projection
  
  int cutSize = (int)cutInPix;
  
  gal_list *gList;
  gal_node *gNode;
  int i, j, k, jN1;
  
  for (j=j_h-cutSize; j<=j_h+cutSize; j++) {
    if (j < 0 || j >= N2) continue;
    jN1 = j * N1;
    for (i=i_h-cutSize; i<=i_h+cutSize; i++) {
      if (i < 0 || i >= N1) continue;
      gList = gMap->map[i+jN1];
      for (k=0, gNode=gList->first; k<gList->size; k++, gNode=gNode->next) {
	lensingForPair(chPar, pkPar, h, gNode->g, err);
	forwardError(*err, __LINE__,);
      }
    }
  }
  return;
}

void lensingForMap(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, gal_map *gMap, error **err)
{
  //-- Lensing main function
  
  int N1 = hMap->N1;
  int N2 = hMap->N2;
  int count = 0;
  int verbose = pkPar->verbose;
  
  halo_list *hList;
  halo_node *hNode;
  int i_h, j_h, k;
  
  for (j_h=0; j_h<N2; j_h++) {
    for (i_h=0; i_h<N1; i_h++) {
      hList = hMap->map[i_h+j_h*N1];
      
      for (k=0, hNode=hList->first; k<hList->size; k++, hNode=hNode->next) {
	if (verbose < 2 && count % 50 == 0) {
	  printf("Computing lensing signal: %6.2f%%\r", 100.0 * count / (double)hMap->total);
	  fflush(stdout);
	}
	lensingForHalo(chPar, pkPar, hNode->h, gMap, i_h, j_h, err); forwardError(*err, __LINE__,);
	count++;
      }
    }
  }
  
  if (pkPar->verbose < 3) printf("Computed lensing signal          \n");
  if (pkPar->doKappa == 1)      gMap->type = kappa_map;
  else if (pkPar->doKappa >= 2) gMap->type = g_map;
  else                          gMap->type = gamma_map;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to galaxy catalogue

void subtractMean(peak_param *pkPar, gal_map *gMap, interpolator_t *k1Inter)
{
  double mean      = 0.0;
  double totWeight = 0.0;
  gal_list *gList;
  gal_node *gNode;
  gal_t *g;
  double kappa1;
  int i, k;

  //-- Compute kappa_mean of all galaxies and subtract
  if (pkPar->doSubtraction == 1) {
    //-- Get mean
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (k=0, gNode=gList->first; k<gList->size; k++, gNode=gNode->next) {
        g = gNode->g;
        mean += g->kappa * g->weight;
      }
      totWeight += gList->totWeight;
    }

    mean /= totWeight;
    gMap->kappa_mean = mean;

    //-- Subtract
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (k=0, gNode=gList->first; k<gList->size; k++, gNode=gNode->next) {
        gNode->g->kappa -= mean;
      }
    }
  }

  //-- Compute the mass sheet kappa_1 at z_s of each galaxy and subtract
  else if (pkPar->doSubtraction == 2) {
    //-- Subtract and compute kappa_mean
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (k=0, gNode=gList->first; k<gList->size; k++, gNode=gNode->next) {
        g = gNode->g;
        kappa1    = execute_interpolator_t(k1Inter, g->z, 1); //-- border = 1
        g->kappa -= kappa1;
        mean     += g->kappa * g->weight;
      }
      totWeight += gList->totWeight;
    }

    //-- Stock kappa_mean
    mean /= totWeight;
    gMap->kappa_mean = mean;
  }
  
  if (pkPar->doSubtraction == 1) {
    if      (pkPar->verbose < 3)   printf("Subtracted kappa mean %g\n", mean);
    else if (pkPar->verbose == 3) {printf("kappa mean = %.5f, ", mean); fflush(stdout);}
  }
  else if (pkPar->doSubtraction == 2) {
    if      (pkPar->verbose < 3)   printf("Subtracted kappa_1\n");
    else if (pkPar->verbose == 3) {printf("subtracted kappa_1, "); fflush(stdout);}
  }
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

void addNoiseToGalaxies(peak_param *pkPar, gal_map *gMap)
{
  if (pkPar->doNoise == 0) {
    if (pkPar->verbose < 3) printf("No noise added\n");
    return;
  }
  
  double sigma_half  = pkPar->sigma_half;
  gsl_rng *generator = pkPar->generator;
  
  gal_list *gList;
  gal_node *gNode;
  double *gamma;
  double g1, g2, e1, e2, A, B, C, g1_sq, g2_sq, factor;
  int i, j;
  
  if (pkPar->doKappa == 1) {
    for (i=0; i<gMap->length; i++) {
      gList = gMap->map[i];
      for (j=0, gNode=gList->first; j<gList->size; j++, gNode=gNode->next) {
	gNode->g->kappa += gsl_ran_gaussian(generator, sigma_half);
      }
    }
  }
  else if (pkPar->doKappa >= 2) {
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
	
	do e1 = gsl_ran_gaussian(generator, sigma_half);
	while (fabs(e1) >= 1.0);
	do e2 = gsl_ran_gaussian(generator, sigma_half);
	while (fabs(e2) >= 1.0);
	
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
	do e1 = gsl_ran_gaussian(generator, sigma_half);
	while (fabs(e1) >= 1.0);
	do e2 = gsl_ran_gaussian(generator, sigma_half);
	while (fabs(e2) >= 1.0);
	
	gamma = gNode->g->gamma;
	gamma[0] += e1;
	gamma[1] += e2;
      }
    }
  }
  
  gMap->type += 2;
  if (pkPar->verbose < 3) printf("Added noise to galaxies\n");
  return;
}

void lensingCatalogue(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, gal_map *gMap, interpolator_t *k1Inter, error **err)
{
  if (pkPar->doLensing == 0) return;

  //-- Lensing
  lensingForMap(chPar, pkPar, hMap, gMap, err); forwardError(*err, __LINE__,);

  //-- Subtract mean
  if (pkPar->doKappa > 0 && pkPar->doSubtraction) subtractMean(pkPar, gMap, k1Inter);
  
  //-- gamma to g
  if (pkPar->doKappa >= 2) makeG(gMap);
  
  //-- Add noise
  if (pkPar->doNoise) addNoiseToGalaxies(pkPar, gMap);
  return;
}

void lensingCatalogueAndOutput(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, gal_map *gMap, interpolator_t *k1Inter, error **err)
{
  if (pkPar->doLensing == 0) return;
  
  char name[STRING_LENGTH_MAX];
  
  //-- Lensing
  lensingForMap(chPar, pkPar, hMap, gMap, err);
  forwardError(*err, __LINE__,);
  
  //-- Subtract mean
  if (pkPar->doKappa > 0 && pkPar->doSubtraction) subtractMean(pkPar, gMap, k1Inter);
  
  //-- gamma to g
  if (pkPar->doKappa >= 2) makeG(gMap);
  
  if (pkPar->doFITS == 0) {
    sprintf(name, "%sgalCat_noiseFree", pkPar->prefix);
    outAsciiGalCat(name, chPar, pkPar, gMap, err); forwardError(*err, __LINE__,);
  }
  else {
    sprintf(name, "%sgalCat_noiseFree.fits", pkPar->prefix);
    outFitsGalCat(name, chPar, pkPar, gMap);
  }
  
  //-- Add noise in makeMapsAndOutput
  return;
}

//----------------------------------------------------------------------
//-- Functions related to galaxy output

void outAsciiGalaxyInfo_hm(FILE *file, cosmo_hm *chPar)
{
  fprintf(file, "# z_gal_min = %.3f, z_gal_max = %.3f\n", chPar->redshift->par_nz[0], chPar->redshift->par_nz[1]);
  fprintf(file, "# alpha_gal = %.3f, beta_gal = %.3f, z_gal_0 = %.3f\n", chPar->redshift->par_nz[2], chPar->redshift->par_nz[3], chPar->redshift->par_nz[4]);
  return;
}

void outAsciiGalaxyInfo_pk(FILE *file, peak_param *pkPar)
{
  fprintf(file, "# z_s = %g, dz_gal = %g, doRandGalPos = %d, n_gal = %g [arcmin^-2], sigma_eps = %g\n", 
	  pkPar->z_s, pkPar->dz_gal, pkPar->doRandGalPos, pkPar->n_gal/RADIAN_SQ_TO_ARCMIN_SQ, pkPar->sigma_eps);
  return;
}

void outAsciiLensingInfo(FILE *file, peak_param *pkPar)
{
  if (pkPar->doLensing) fprintf(file, "# doLensing = %d, doKappa = %d, doSubtraction = %d\n", pkPar->doLensing, pkPar->doKappa, pkPar->doSubtraction);
  else                  fprintf(file, "# doLensing = %d\n", pkPar->doLensing);
  return;
}

void outAsciiGalCat(char name[], cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err)
{
  if (pkPar->outGalCat == 0) return;
  
  FILE *file = fopen_err(name, "w", err); forwardError(*err, __LINE__,);
  
  fprintf(file, "# Galaxy catalogue\n");
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, chPar, pkPar);
  fprintf(file, "#\n");
  fprintf(file, "# Halo model parameters = %s\n", pkPar->hmParPath);
  outAsciiGalaxyInfo_hm(file, chPar);
  fprintf(file, "#\n");
  fprintf(file, "# Peak parameters = %s\n", pkPar->pkParPath);
  outAsciiFieldInfo(file, pkPar);
  outAsciiGalaxyInfo_pk(file, pkPar);
  outAsciiLensingInfo(file, pkPar);
  fprintf(file, "#\n");
  
  outAscii_gal_map(file, pkPar, gMap);
  
  fclose(file);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
  return;
}

#ifdef __CAMELUS_USE_FITS__
void outFitsGalaxyInfo(FITS_t *fits, cosmo_hm *chPar, peak_param *pkPar)
{
  addLineSpread(fits);
  
  double buffer_dbl;
  
  addKeyword(fits, TDOUBLE, "ZGALMIN",  &chPar->redshift->par_nz[0], "[-] Minimum galaxy redshift");
  addKeyword(fits, TDOUBLE, "ZGALMAX",  &chPar->redshift->par_nz[1], "[-] Maximum galaxy redshift");
  addKeyword(fits, TDOUBLE, "ALPHAGAL", &chPar->redshift->par_nz[2], "[-]");
  addKeyword(fits, TDOUBLE, "BETAGAL",  &chPar->redshift->par_nz[3], "[-]");
  addKeyword(fits, TDOUBLE, "ZGAL0",    &chPar->redshift->par_nz[4], "[-]");
  
  addKeyword(fits, TDOUBLE, "ZS",       &pkPar->z_s,                 "[-] Source redshift");
  addKeyword(fits, TDOUBLE, "DZGAL",    &pkPar->dz_gal,              "[-] Galaxy redshift binwidth");
  addKeyword(fits, TINT,    "RANDPOS",  &pkPar->doRandGalPos,        "[-] 0 = regular, 1 = random");
  buffer_dbl = pkPar->n_gal / RADIAN_SQ_TO_ARCMIN_SQ;
  addKeyword(fits, TDOUBLE, "NGAL",     &buffer_dbl,                 "[arcmin^-2] Galaxy number density");
  addKeyword(fits, TDOUBLE, "SIGMAEPS", &pkPar->sigma_eps,           "[-] Ellipticity dispersion");
  return;
}

void outFitsLensingInfo(FITS_t *fits, peak_param *pkPar)
{
  addLineSpread(fits);
  addKeyword(fits, TINT, "DOLEN",   &pkPar->doLensing,     "[-] 0 = read from inGalCatPath, 1 = compute lensing");
  if (pkPar->doLensing) {
    addKeyword(fits, TINT, "DOKAPPA", &pkPar->doKappa,       "[-] 0 = gamma, 1 = kappa, 2 = g with linear KS, 3 = g with iterative KS, 4 = g with SS");
    addKeyword(fits, TINT, "DOSUB",   &pkPar->doSubtraction, "[-] 0 = without, 1 = subtract kappa_mean, 2 = subtract the mass sheet kappa_1");
  }
  return;
}
#endif

void outFitsGalCat(char name[], cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap)
{
  if (pkPar->outGalCat == 0) return;
  
#ifdef __CAMELUS_USE_FITS__
  FITS_t *fits = initializeTableWriter_FITS_t(name);
  outFits_gal_map(fits, pkPar, gMap);
  
  outFitsCosmoParam(fits, chPar, pkPar);
  
  addLineSpread(fits);
  addKeyword(fits, TSTRING, "HMPARAM", pkPar->hmParPath,           "[-] Path of the halo model parameter file");
  addKeyword(fits, TSTRING, "PKPARAM", pkPar->pkParPath,           "[-] Path of the peak parameter file");
  
  outFitsFieldInfo(fits, pkPar);
  outFitsGalaxyInfo(fits, chPar, pkPar);
  outFitsLensingInfo(fits, pkPar);
  
  free_FITS_t(fits);
  if (pkPar->verbose < 3) printf("Outputed \"%s\"\n", name);
#endif
  return;
}

//----------------------------------------------------------------------
//-- Main function

void doProfile(cosmo_hm *chPar, peak_param *pkPar, double z_l, double M, double z_s, error **err)
{
  printf("z_l = %f\n", z_l);
  printf("M   = %11.5e [M_sol/h]\n", M);
  printf("z_s = %f\n", z_s);
  
  pkPar->z_halo_max = z_s;
  pkPar->dz_halo    = pkPar->z_s / (double)(pkPar->N_z_halo);
  
  double theta_maxx = 200.0; //-- [arcmin]
  double dtheta     = 0.01;  //-- [arcmin]
  double weight     = (pkPar->sigma_half > 0.0) ? 1.0 / (SQ(pkPar->sigma_half)) : 1.0;
  int N_gal         = (int)round(theta_maxx / dtheta);
  double pos[2]     = {0.0, 0.0};
  int i;
  
  halo_t *h = (halo_t*)malloc_err(sizeof(halo_t), err);   forwardError(*err, __LINE__,);
  set_halo_t(chPar, pkPar, h, pos, z_l, -1.0, M, err);    forwardError(*err, __LINE__,);
  
  gal_t *g  = (gal_t*)malloc_err(sizeof(gal_t), err);     forwardError(*err, __LINE__,);
  set_gal_t(chPar, pkPar, g, NULL, z_s, weight, err);     forwardError(*err, __LINE__,);
  
  double_mat *profile = initialize_double_mat(N_gal, 4);
  for (i=0; i<N_gal; i++) profile->matrix[i] = (i + 1) * dtheta;
  
  double factor = factorForTwoHaloTerm(chPar, h, g, err); forwardError(*err, __LINE__,);
  fillOneHaloTerm(chPar, pkPar, h, g, profile, err);      forwardError(*err, __LINE__,);
  fillTwoHaloTerm(chPar, pkPar, h, factor, profile, err); forwardError(*err, __LINE__,);
  
  char name[STRING_LENGTH_MAX];
  sprintf(name, "%sprofile_zl%.3f_logM%.2f_zs%.3f", pkPar->prefix, z_l, log10(M), z_s);
  outAsciiProfile(name, chPar, pkPar, h, profile, err);   forwardError(*err, __LINE__,);
  
  free(h);
  free(g);
  free_double_mat(profile);
  printf("------------------------------------------------------------------------\n");
  return;
}

// MKDEBUG: removed const from halo_map (to avoid warning); new k1Inter for subtractMean
void lensingCatalogueAndOutputAll2(char fileName[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_map *gMap, interpolator_t *k1Inter, error **err)
{

  //-- Lensing
  lensingForMap(cmhm, peak, hMap, gMap, err); forwardError(*err, __LINE__,);
	
  //-- Subtract mean
  if (peak->doKappa != 0) subtractMean(peak, gMap, k1Inter);
  
  //-- gamma to g
  if (peak->doKappa >= 2) makeG(gMap);
  // MKDEBUG replaced according to Linc-tw
  //outputGalaxies(fileName, cmhm, peak, gMap);
  outAsciiGalCat(fileName, cmhm, peak, gMap, err); forwardError(*err, __LINE__,);
  
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

void doRayTracing(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
  //-- Currently unused
  
  sampler_arr *hSampArr   = initialize_sampler_arr(pkPar->N_z_halo, pkPar->N_M+1);
  interpolator_t *k1Inter = initialize_interpolator_t(pkPar->N_z_gal+1);
  setMassSampAndK1Inter(chPar, pkPar, hSampArr, NULL, k1Inter, err);                                      forwardError(*err, __LINE__,); //-- lambda = NULL
  halo_map *hMap          = initialize_halo_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err); forwardError(*err, __LINE__,);
  sampler_t *gSamp        = initialize_sampler_t(pkPar->N_z_gal+1);
  setGalaxySampler(chPar, pkPar, gSamp, err);                                                             forwardError(*err, __LINE__,);
  gal_map *gMap           = initialize_gal_map(pkPar->resol[0], pkPar->resol[1], pkPar->theta_pix, err);  forwardError(*err, __LINE__,);
  mask_map *mask          = initializeMask(pkPar, err);                                                   forwardError(*err, __LINE__,);
  
  char name[STRING_LENGTH_MAX];
  
  readCatOrMakeSimulAndOutput(chPar, pkPar, hSampArr, hMap, err);    forwardError(*err, __LINE__,);
  cleanOrMakeOrResample(chPar, pkPar, gSamp, gMap, mask, err);       forwardError(*err, __LINE__,);
  lensingCatalogueAndOutput(chPar, pkPar, hMap, gMap, k1Inter, err); forwardError(*err, __LINE__,);
  
  free_sampler_arr(hSampArr);
  free_interpolator_t(k1Inter);
  free_halo_map(hMap);
  free_sampler_t(gSamp);
  free_gal_map(gMap);
  free_mask_map(mask);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------


void read_gal_map2(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err)
{
  //-- WARNING: D_S, w_s are not read.
  
  double factor = (peak->field == 0) ? ARCMIN_TO_RADIAN : DEGREE_TO_RADIAN;

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
      buffer2 = fscanf(file, "%lf %lf %lf %lf %lf %lf\n", &pos[0], &pos[1], &z, &kappa, gamma, gamma+1);   
      pos[0] *= factor;
      pos[1] *= factor;
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
  g->w = w(cmhm->cosmo, g->a, 0, err);       forwardError(*err, __LINE__,); //-- wOmegar = 0
  g->D_s = g->a * f_K(cmhm->cosmo, g->w, err); forwardError(*err, __LINE__,);
  return;
}


//----- new func-----

void outputFastSimul_galaxies(char name_cmhm[], char name[], char name2[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap)
{
  error *myerr = NULL, **err = &myerr;

  gal_map *gMap = initialize_gal_map(hMap->N1,hMap->N2,hMap->theta_pix,err); forwardError(*err, __LINE__,);
  FILE *file = fopen(name, "w");
  FILE *file2 = fopen(name2, "w");

  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");

  fprintf(file2, "# Galaxy catalogue\n");
  fprintf(file2,  "# Field = %s, Omega = (%g, %g) [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file2, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  fprintf(file2, "#\n");
  outAsciiCosmoParam(file2, cmhm, peak);
  fprintf(file2, "#\n");

  //printf("test3 \n");
  output_halo_map_galaxies(file,file2,cmhm, peak, hMap, gMap);
  //printf("test4 \n");
  free_gal_map(gMap);
  //printf("Gfreen \n");
  fclose(file);
  fclose(file2);
  //printf("close \n");
  printf("\"%s\" made\n", name);
  printf("\"%s\" made\n", name2);
  return;
}


void outputFastSimul_galaxies2(char name_cmhm[], char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_map *gMap)
{

  error *myerr = NULL, **err = &myerr;

  FILE *file = fopen(name, "w");

  fprintf(file, "# Halo list, fast simulation\n");
  fprintf(file, "# Model = %s, field = %s, Omega = (%g, %g) [arcmin]\n", smassfct_t(cmhm->massfct), STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1]);
  fprintf(file, "# z_halo_max = %g, N_z_halo = %d, M_min = %8.2e [M_sol/h], M_max = %8.2e\n", peak->z_halo_max, peak->N_z_halo, peak->M_min, peak->M_max);
  fprintf(file, "#\n");
  outAsciiCosmoParam(file, cmhm, peak);
  fprintf(file, "#\n");
  output_halo_map_galaxies2(file,cmhm, peak, hMap, gMap);
  fclose(file);
  printf("\"%s\" made\n", name);

  return;
}


void output_halo_map_galaxies(FILE *file,FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_map *gMap)
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

  // MKDEBUG new, following Linc-tw:doProfile
  double weight     = (peak->sigma_half > 0.0) ? 1.0 / (SQ(peak->sigma_half)) : 1.0;

  //printf("l = %i  \n",hMap->length);
  
  fprintf(file2, "#  theta_x   theta_y      z        \n  ");
  fprintf(file2,"# [arcmin]  [arcmin]      [-]       \n");

  for (i=0; i<hMap->length; i++) {
    hList = hMap->map[i];
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
      peak->w_s = h->w;
      peak->D_s = Ds;

      if( (rand()/(double)RAND_MAX)<ngc) {
        // MKDEBUG diff. arg. in Linc-tw
			append_gal_map(cmhm, peak, gMap, h->pos, h->z, weight, err); forwardError(*err, __LINE__,);
			fprintf(file2, "%9.3f  %9.3f   %7.5f   \n", h->pos[0], h->pos[1], h->z); // TRD version
        //fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", h->pos[0], h->pos[1], h->z, j);  // MKDEBUG new version?
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
        // MKDEBUG diff. arg. in Linc-tw
        //append_gal_list(cmhm, gList, h->z, h->w, Ds,h->pos, err); forwardError(*err, __LINE__,);
        append_gal_map(cmhm, peak, gMap, h->pos, h->z, weight, err); forwardError(*err, __LINE__,);
        //fprintf(file2, "%9.3f  %9.3f   %7.5f  %i  \n", pos[0], pos[1], h->z, j);
        fprintf(file2, "%9.3f  %9.3f   %7.5f\n", pos[0], pos[1], h->z);
      }

      //printf("ok \n");
    }
  }
  printf("Nb galaxies created : %i \n",ii);
  return;
}

double NFW(double x)
{
  return 1/(x*(1+pow(x,2)));
}

void output_halo_map_galaxies2(FILE *file, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_map *gMap)
{
  halo_list *hList;
  halo_node *hNode;
  error *myerr = NULL, **err = &myerr;
  int i,j,ii,ii2,k;
  double Ds,Mh;

  // MKDEBUG new, following Linc-tw:doProfile
  double weight     = (peak->sigma_half > 0.0) ? 1.0 / (SQ(peak->sigma_half)) : 1.0;

  ii=0;
  ii2=0;
  srand(time(NULL));
  printf( "Number of halos = %d\n", hMap->total);

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
  halo_t *h;

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
      peak->w_s = h->w;
      peak->D_s = Ds;
      if( (rand()/(double)RAND_MAX)<ngc) {
        // MKDEBUG diff. arg. in Linc-tw
        append_gal_map(cmhm, peak, gMap, h->pos, h->z, weight, err); forwardError(*err, __LINE__,);
        fprintf(file, "%9.3f  %9.3f   %7.5f   \n", h->pos[0], h->pos[1], h->z);
		  ii2=ii2+1;
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
        double posg[2];
        posg[0] = cos(theta) * sin(phi) * r + h->pos[0];
        posg[1] = sin(theta) * sin(phi) * r + h->pos[1];
        // MKDEBUG diff. arg. in Linc-tw
		  if((posg[0] > 180)||(posg[0]<0)||(posg[1] > 180)||(posg[1]<0)){
		  }else{
			  append_gal_map(cmhm, peak, gMap, h->pos, h->z, weight, err); forwardError(*err, __LINE__,);
			  fprintf(file, "%9.3f  %9.3f   %7.5f   \n", h->pos[0], h->pos[1], h->z);
			  ii2=ii2+1;
		  }
      }
    }
  }
  printf("Nb galaxies created : %i \n",ii);
  printf("Nb galaxies created arrondi : %i \n",ii2);
  return;
}
