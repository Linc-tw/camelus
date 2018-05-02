

  /*******************************************************
   **  haloSampling.h					**
   **  Version 2018.03.14				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_HALO_SAMPLING__
#define __CAMELUS_HALO_SAMPLING__

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#ifdef __CAMELUS_USE_HEALPIX__
#include "HEALPixFunctions.h"
#endif

#include "parameters.h"


typedef struct {
  double pos[2];       //-- [rad] Angular position
  double z;            //-- [-] Halo redshift
  double a;            //-- [-] 1 / (1+z)
  double w;            //-- [Mpc/h] Comoving radial distance to halo
  double M;            //-- [M_sol/h] Halo mass
  double c;            //-- [-] Halo concentration
  double c_sq;         //-- [-] c^2
  double f;            //-- [-] 1 / (ln(1 + c) - c/(1 + c))
  double r_vir;        //-- [Mpc/h] Physical virial radius
  double r_vir_sq;     //-- [Mpc^2/h^2]
  double theta_vir;    //-- [rad] theta_vir = r_vir / D_A = angular apparent size
  double theta_vir_sq; //-- [rad^2]
  double factor;       //-- factor = FOUR_PI_G_OVER_C2 * (D_A M f c_NFW^2 / 2 pi r_vir^2)
  double sinCosDEC[2]; //-- [-] sinDEC & cosDEC, only used if field > 0
  double n_gal;        //-- N galaxies inside the halo as predicted by HOD model
} halo_t;

typedef struct halo_node {
  halo_t *h;
  struct halo_node *next;
} halo_node;

typedef struct {
  int length;       //-- Number of nodes allocated
  int size;         //-- Number of nodes containing information
  halo_node *first; //-- Begin of the list
  halo_node *last;  //-- End of the list
} halo_list;

typedef struct {
  int N1, N2;           //-- Number of pixels of each side
  int length;           //-- Number of pixels
  int total;            //-- Number of halos
  double theta_pix;     //-- [rad] Pixel size
  double theta_pix_inv; //-- [rad^-1] Inverse of pixel size
  halo_list **map;      //-- Binned halos
} halo_map;



//-- Functions related to halo_t, halo_node, halo_list
void set_halo_t(cosmo_hm *chPar, peak_param *pkPar, halo_t *h, double pos[2], double z, double ww, double M, error **err);
halo_node *initialize_halo_node(error **err);
halo_list *initialize_halo_list(error **err);
void free_halo_list(halo_list *hList);
void reset_halo_list(halo_list *hList);
void append_halo_list(cosmo_hm* chPar, peak_param *pkPar, halo_list* hList, double pos[2], double z, double ww, double M, error** err);

//-- Functions related to halo_map
halo_map *initialize_halo_map(int N1, int N2, double theta_pix, error **err);
void free_halo_map(halo_map *hMap);
void reset_halo_map(halo_map *hMap);
int append_halo_map(cosmo_hm* chPar, peak_param* pkPar, halo_map* hMap, double pos[2], double z, double ww, double M, error** err);
void read_halo_map(char name[], cosmo_hm* chPar, peak_param* pkPar, halo_map* hMap, error** err);



//-- Functions related to mass function
double massFct(cosmo_hm_params *cANDp, double mass, error **err);
void fillMassFct(cosmo_hm_params *cANDp, sampler_t *hSamp, error **err);
void outAsciiMassFct(char name[], cosmo_hm *chPar, peak_param *pkPar, double z, error **err);
void outFitsMassFct(char name[], cosmo_hm *chPar, peak_param *pkPar, double z, error **err);

//-- Functions related to halo output
void outAscii_halo_map(FILE *file, peak_param *pkPar, halo_map *hMap);
void outAsciiFieldInfo(FILE *file, peak_param *pkPar);
void outAsciiHaloInfo(FILE *file, peak_param *pkPar);
void outAsciiHaloCat(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap, error **err);
#ifdef __CAMELUS_USE_FITS__
void outFits_halo_t(FITS_t *fits, halo_t *h, double factor);
void outFits_halo_map(FITS_t *fits, peak_param *pkPar, halo_map *hMap);
void outFitsFieldInfo(FITS_t *fits, peak_param *pkPar);
void outFitsHaloInfo(FITS_t *fits, peak_param *pkPar);
#endif
void outFitsHaloCat(char name[], cosmo_hm *chPar, peak_param *pkPar, halo_map *hMap);

//-- Functions related to fast simulations
double dVol(cosmo_hm* chPar, double z, double ww, double dz, double area, error** err);
void randomizePosition(peak_param *pkPar, double pos[2], error **err);
void setMassSamplers(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, int doVolume, error **err);
void rescaleWithVolume(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, error **err);
void sampleHalos(cosmo_hm *chPar, peak_param *pkPar, sampler_t *hSamp, halo_map *hMap, double z1, double z2, error **err);
void makeFastSimul(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, halo_map *hMap, error **err);
void readCatOrMakeSimulAndOutput(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, halo_map *hMap, error **err);

//-- Functions related to mass sheet
void setLambda(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, double_arr *lambda);
void massSheet(cosmo_hm *chPar, peak_param *pkPar, sampler_t *hSamp, double sheet[2], error **err);
void setKappa1Interpolator(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, interpolator_t *k1Inter, double *k0Arr, error **err);
void setMassSampAndK1Inter(cosmo_hm *chPar, peak_param *pkPar, sampler_arr *hSampArr, double_arr *lambda, interpolator_t *k1Inter, error **err);
void outAsciiLambda(char name[], cosmo_hm *chPar, peak_param *pkPar, double_arr *lambda, error **err);
void outAsciiMassSheet(char name[], cosmo_hm *chPar, peak_param *pkPar, interpolator_t *k1Inter, double *k0Arr, error **err);

//-- Main functions
void doMassFct(cosmo_hm *chPar, peak_param *pkPar, double z, error **err);
void doMassSheet(cosmo_hm *chPar, peak_param *pkPar, double z_s, error **err);
void doFastSimulation(cosmo_hm *chPar, peak_param *pkPar, error **err);

/// NEW HOD
void outputFastSimul_HOD(char name_cmhm[],char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);
void output_halo_map_HOD(char name_cmhm[],FILE *file, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);





#endif

