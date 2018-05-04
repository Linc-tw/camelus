

  /*******************************
   **  haloSampling.h		**
   **  Chieh-An Lin		**
   **  Version 2015.10.30	**
   *******************************/


#ifndef __haloSamp__
#define __haloSamp__

#include "commonHeader.h"
#include "peakParameters.h"
#include <nicaea/hod.h>
#include <nicaea/halomodel.h>


typedef struct {
  double pos[2];       //-- [arcmin/deg] Angular position, (theta_x, theta_y) in arcmin or (RA, DEC) in deg
  double z;            //-- [-] Halo redshift
  double a;            //-- [-] 1 / (1+z)
  double w;            //-- [Mpc/h] Comoving radial distance to halo
  double M;            //-- [M_sol/h] Halo mass
  double c;            //-- [-] Halo concentration
  double c_sq;         //-- [-] c^2
  double f;            //-- [-] 1 / (ln(1 + c) - c/(1 + c))
  double r_vir;        //-- [Mpc/h] Physical virial radius
  double r_vir_sq;     //-- [Mpc^2/h^2]
  double theta_vir;    //-- [arcmin] theta_vir = r_vir / D_A = angular apparent size
  double theta_vir_sq; //-- [arcmin^2]
  double factor;       //-- factor = FOUR_PI_G_OVER_C2 * (D_A M f c_NFW^2 / 2 pi r_vir^2)
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
  int N1, N2;       //-- Number of pixels of each side
  int length;       //-- Number of pixels
  int total;        //-- Number of halos
  double theta_pix; //-- [arcmin/deg] Pixel size
  double theta_pix_inv; //-- [arcmin^-1/deg^-1] Inverse of pixel size
  double limits[4]; //-- [arcmin/deg] (theta_x_min, theta_x_max, theta_y_min, theta_y_max) in arcmin or (RA_min, RA_max, DEC_min, DEC_max) in deg
  double center[2]; //-- [arcmin/deg] Center of the map;
  halo_list **map;  //-- Binned halos
} halo_map;



//-- Functions related to halo_t, halo_node, halo_list
void set_halo_t(cosmo_hm *cmhm, halo_t *h, double z, double M, double pos[2], error **err);
halo_node *initialize_halo_node(error **err);
halo_list *initialize_halo_list(error **err);
void free_halo_list(halo_list *hList);
void reset_halo_list(halo_list *hList);
void append_halo_list(cosmo_hm *cmhm, halo_list *hList, double z, double M, double pos[2], error **err);

//-- Functions related to halo_map
halo_map *initialize_halo_map(int N1, int N2, double theta_pix, error **err);
void free_halo_map(halo_map *hMap);
void reset_halo_map(halo_map *hMap);
void append_halo_map(cosmo_hm *cmhm, halo_map *hMap, double z, double M, double pos[2], error **err);
void read_halo_map(char name[], cosmo_hm *cmhm, halo_map *hMap, error **err);
void read_halo_map2(char name[], cosmo_hm *cmhm, halo_map *hMap, error **err);
void output_halo_map(FILE *file, peak_param *peak, halo_map *hMap);



//-- Functions related to mass function
double massFct(cosmo_hm_params *cANDp, double mass, error **err);
void fillMassFct(cosmo_hm_params *cANDp, sampler_t *samp, error **err);
void outputMassFct(char name[], cosmo_hm *cmhm, peak_param *peak, double z, error **err);

//-- Functions related to fast simulations
double dVol(cosmo_hm *cmhm, peak_param *peak, double z, double ww, double dz, error **err);
void randomizePosition(peak_param *peak, double *pos, error **err);
void setMassSamplers(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, error **err);
void sampleHalos(cosmo_hm *cmhm, peak_param *peak, sampler_t *samp, halo_map *hMap, double z, error **err);
void makeFastSimul(cosmo_hm *cmhm, peak_param *peak, sampler_arr *sampArr, halo_map *hMap, error **err);
void outputFastSimul(char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);

//-- Functions related to mass sheet
void makeAInterpolator(cosmo_hm *cmhm, interpolator_t *a_inter, double z_max, error **err);
void massSheet(cosmo_hm *cmhm, peak_param *peak, sampler_t *samp, interpolator_t *a_inter, double sheet[2], error **err);

//-- Main functions
void doMassFct(cosmo_hm *cmhm, peak_param *peak, error **err);
void doFastSimulation(cosmo_hm *cmhm, peak_param *peak, error **err);
void doMassSheet(cosmo_hm *cmhm, peak_param *peak, double z_halo_max, double M_min, error **err);


/// NEW HOD
void outputFastSimul_HOD(char name_cmhm[],char name[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);
void output_halo_map_HOD(char name_cmhm[],FILE *file, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);





#endif

