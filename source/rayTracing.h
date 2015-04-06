

  /*******************************
   **  rayTracing.h		**
   **  Chieh-An Lin		**
   **  Version 2015.02.25	**
   *******************************/


#ifndef __rayTrac__
#define __rayTrac__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"


typedef struct {
  double pos[2];   //-- [arcmin/deg] Angular position, (theta_x, theta_y) in arcmin or (RA, DEC) in deg
  double z;        //-- [-] Galaxy redshift
  double a;        //-- [-] 1 / z
  double w;        //-- [Mpc/h] Comoving radial distance to galaxy
  double D_s;      //-- [Mpc/h] Angular diameter distance
  double kappa;    //-- [-] Convergence
  double gamma[2]; //-- [-] Shear
  double e[2];     //-- [-] Measured ellipticies
  
  //-- Belows are only used in plugIn.c, to drop in later version
  double weight;   //-- Weight in irrgular smoothing
                   //--   0.0 = unvisited or revisited
                   //--   otherwise, visited
  int indices[2];  //-- Only used for HEALPix grid
} gal_t;

typedef struct gal_node {
  gal_t *g;
  struct gal_node *next;
} gal_node;

typedef struct {
  int length;      //-- Number of nodes allocated
  int size;        //-- Number of nodes containing information
  gal_node *first; //-- Begin of the list
  gal_node *last;  //-- End of the list
} gal_list;

typedef struct {
  int N1, N2;       //-- Number of pixels of each side
  int length;       //-- Number of pixels
  int total;        //-- Number of galaxies
  double theta_pix; //-- [arcmin/deg] Pixel size
  double center[2]; //-- [arcmin/deg] Center of the map;
  double limits[4]; //-- [arcmin/deg] (theta_x_min, theta_x_max, theta_y_min, theta_y_max) in arcmin or (RA_min, RA_max, DEC_min, DEC_max) in deg
  gal_list **map;   //-- Binned galaxies
} gal_map;


//-- Functions related to gal_t, gal_node, gal_list
void set_gal_t(cosmo_hm *cmhm, gal_t *g, double z, double w_s, double D_s, double *pos, error **err);
gal_node *initialize_gal_node(error **err);
gal_list *initialize_gal_list(error **err);
void free_gal_list(gal_list *gList);
void cleanLensing_gal_list(gal_list *gList);
void reset_gal_list(gal_list *gList);
void append_gal_list(cosmo_hm *cmhm, gal_list *gList, double z, double w_s, double D_s, double pos[2], error **err);
double kappaMean_gal_list(gal_list *gList, error **err);

//-- Functions related to gal_map
gal_map *initialize_gal_map(int N1, int N2, double theta_pix, error **err);
void free_gal_map(gal_map *gMap);
void cleanLensing_gal_map(gal_map *gMap);
void reset_gal_map(gal_map *gMap);
void append_gal_map(cosmo_hm *cmhm, gal_map *gMap, double z, double w_s, double D_s, double pos[2], error **err);
void output_gal_map(FILE *file, peak_param *peak, gal_map *gMap);
void updateCosmo_gal_map(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);

//-- Functions related to projected NFW mass
double G_NFW_kappa(double x_sq, double c, double c_sq, error **err);
double G_NFW_gamma(double x_sq, double c, double c_sq, double f, error **err);
double kappa_NFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err);
void gamma_NFW(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err);

//-- Functions related to lensing
void lensingForPair(cosmo_hm *cmhm, halo_t *h, gal_t *g, int doKappa, error **err);
void lensingForHalo(cosmo_hm *cmhm, halo_t *h, gal_map *gMap, int doKappa, error **err);
void lensingForMap(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err);

//-- Functions related to making galaxy lists
void makeGalaxiesFixedRedshift(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void fillGalaxyLaw(redshiftANDint *rANDi, sampler_t *samp, error **err);
void sampleGalaxies(redshiftANDint *rANDi, cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void makeRealGalaxies(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void makeGalaxies(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void outputGalaxies(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap);

//-- Main function
void doRayTracing(char haloList[], cosmo_hm *cmhm, peak_param *peak, error **err);


#endif

