

  /*******************************************************
   **  rayTracing.h					**
   **  Version 2018.03.07				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

// MKDEBUG New for STR_MAP...
#include "parameters.h"

#ifndef __CAMELUS_RAY_TRACING__
#define __CAMELUS_RAY_TRACING__

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#include "parameters.h"
#include "haloSampling.h"
#include "galaxySampling.h"


typedef struct {
  double f_K, a, theta;
  cosmo *cosmo;
} profile_inteParam;


//-- Functions related to projected mass
double G_NFW_kappa(double u_sq, double c, double c_sq, error **err);
double G_NFW_gamma(double u_sq, double c, double c_sq, double f, error **err);
void G_NFW_both(double u_sq, double c, double c_sq, double f, double G_NFW[2], error **err);
double G_TJ_kappa(double u_sq, double c, double c_sq, error **err);
double G_TJ_gamma(double u_sq, double c, double c_sq, double f, error **err);
void G_TJ_both(double u_sq, double c, double c_sq, double f, double G_TJ[2], error **err);
double G_BMO_kappa(double u_sq, double tau, double tau_sq, error **err);
double kappa_TJ(cosmo_hm *chPar, halo_t *h, gal_t *g, double theta_sq, error **err);
void gamma_TJ(cosmo_hm *chPar, halo_t *h, gal_t *g, double theta_sq, double phase, error **err);
void both_TJ(cosmo_hm *chPar, halo_t *h, gal_t *g, double theta_sq, double phase, error **err);

//-- Functions related to drawing a profile
void fillOneHaloTerm(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, gal_t *g, double_mat *profile, error **err);
double factorForTwoHaloTerm(cosmo_hm *chPar, halo_t *h, gal_t *g, error **err);
double integrandForTwoHaloTerm(double l, void *inteParam, error **err);
double twoHaloTerm(cosmo_hm *chPar, halo_t *h, double theta, double factor, error **err);
void fillTwoHaloTerm(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, double factor, double_mat *profile, error **err);
void outAsciiProfile(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_t *h, double_mat *profile, error **err);

//-- Functions related to lensing
void lensingForPair(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, gal_t *g, error **err);
void lensingForHalo(cosmo_hm* chPar, peak_param* pkPar, halo_t* h, gal_map* gMap, int i_h, int j_h, error** err);
void lensingForMap(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, gal_map *gMap, error **err);

//-- Functions related to galaxy catalogue
void subtractMean(peak_param *pkPar, gal_map *gMap, interpolator_t *k1Inter);
void makeG(gal_map *gMap);
void addNoiseToGalaxies(peak_param *pkPar, gal_map *gMap);
void lensingCatalogue(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, gal_map *gMap, interpolator_t *k1Inter, error **err);
void lensingCatalogueAndOutput(cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, gal_map *gMap, interpolator_t *k1Inter, error **err);

//-- Functions related to galaxy output
void outAsciiGalaxyInfo_hm(FILE *file, cosmo_hm *chPar);
void outAsciiGalaxyInfo_pk(FILE *file, peak_param *pkPar);
void outAsciiLensingInfo(FILE *file, peak_param *pkPar);
void outAsciiGalCat(char name[], cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err);
#ifdef __CAMELUS_USE_FITS__
void outFitsGalaxyInfo(FITS_t *fits, cosmo_hm *chPar, peak_param *pkPar);
void outFitsLensingInfo(FITS_t *fits, peak_param *pkPar);
#endif
void outFitsGalCat(char name[], cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap);

//-- Main function
void doProfile(cosmo_hm *chPar, peak_param *pkPar, double z_l, double M, double z_s, error **err);
void doRayTracing(cosmo_hm *chPar, peak_param *pkPar, error **err);

//-- New function
void lensingCatalogueAndOutputAll2(char fileName[],cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_map *gMap, interpolator_t *k1Inter, error **err);


void read_gal_map2(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void appendWithSignal_gal_map2(cosmo_hm *cmhm, gal_map *gMap, double z, double pos[2], error **err);
void appendWithSignal_gal_list2(cosmo_hm *cmhm, gal_list *gList, double z, double pos[2], error **err);
void setWithSignal_gal_t2(cosmo_hm *cmhm, gal_t *g, double z, double pos[2], error **err);


void outputFastSimul_galaxies(char name_cmhm[], char name[], char name2[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);
void output_halo_map_galaxies(FILE *file, FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_list *gList);


double NFW(double x);
#endif

