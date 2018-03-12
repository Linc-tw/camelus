

  /*******************************
   **  rayTracing.h		**
   **  Chieh-An Lin		**
   **  Version 2016.03.20	**
   *******************************/


#ifndef __rayTrac__
#define __rayTrac__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"

#ifndef __releaseMenu__
  #include "FITSFunctions.h"
#endif


typedef struct {
  double pos[2];   //-- [arcmin/deg] Angular position, (theta_x, theta_y) in arcmin or (RA, DEC) in deg
  double z;        //-- [-] Galaxy redshift
  double a;        //-- [-] 1 / z
  double w;        //-- [Mpc/h] Comoving radial distance to galaxy
  double D_s;      //-- [Mpc/h] Angular diameter distance
  double kappa;    //-- [-] Convergence
  double gamma[2]; //-- [-] Shear
  
  //-- Belows are only used in paperI.c, to drop in later version
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
  double mean[2];  //-- Mean of the signal
} gal_list;

typedef struct {
  int N1, N2;              //-- Number of pixels of each side
  int length;              //-- Number of pixels
  int total;               //-- Number of galaxies
  double theta_pix;        //-- [arcmin/deg] Pixel size
  double theta_pix_inv;    //-- [arcmin^-1/deg^-1] Inverse of pixel size
  double limits[4];        //-- [arcmin/deg] (theta_x_min, theta_x_max, theta_y_min, theta_y_max) in arcmin or (RA_min, RA_max, DEC_min, DEC_max) in deg
  double center[2];        //-- [arcmin/deg] Center of the map;
  gal_list **map;          //-- Binned galaxies
  double kappa_mean;       //-- Mean of kappa over galaxies
  double fillingThreshold; //-- Threshold from filling factor for masking
  mapType_t type;          //-- Map type (defined in peakParameters.h)
} gal_map;

typedef struct {
  double f_K, a, theta;
  cosmo *cosmo;
} profile_inteParam;


//-- Functions related to gal_t, gal_node, gal_list
void set_gal_t(cosmo_hm *cmhm, gal_t *g, double z, double w_s, double D_s, double *pos, error **err);
void setWithSignal_gal_t(cosmo_hm *cmhm, gal_t *g, double z, double pos[2], double kappa, double gamma[2], error **err);
gal_node *initialize_gal_node(error **err);
gal_list *initialize_gal_list(error **err);
void free_gal_list(gal_list *gList);
void cleanLensing_gal_list(gal_list *gList);
void reset_gal_list(gal_list *gList);
void append_gal_list(cosmo_hm *cmhm, gal_list *gList, double z, double w_s, double D_s, double pos[2], error **err);
void appendWithSignal_gal_list(cosmo_hm *cmhm, gal_list *gList, double z, double pos[2], double kappa, double gamma[2], error **err);
void kappaMean_gal_list(gal_list *gList, double threshold);
void gammaOrGMean_gal_list(gal_list *gList, double threshold);

//-- Functions related to gal_map
gal_map *initialize_gal_map(int N1, int N2, double theta_pix, error **err);
void free_gal_map(gal_map *gMap);
void cleanLensing_gal_map(gal_map *gMap);
void reset_gal_map(gal_map *gMap);
void append_gal_map(cosmo_hm *cmhm, gal_map *gMap, double z, double w_s, double D_s, double pos[2], error **err);
void appendWithSignal_gal_map(cosmo_hm *cmhm, gal_map *gMap, double z, double pos[2], double kappa, double gamma[2], error **err);
void setFillingThreshold(peak_param *peak, gal_map *gMap, error **err);
void signalMean_gal_map(peak_param *peak, gal_map *gMap, error **err);
void output_gal_map(FILE *file, peak_param *peak, gal_map *gMap);
void read_gal_map(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void updateCosmo_gal_map(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);

//-- Functions related to projected mass
double G_NFW_kappa(double u_sq, double c, double c_sq, error **err);
double G_NFW_gamma(double u_sq, double c, double c_sq, double f, error **err);
void G_NFW_both(double u_sq, double c, double c_sq, double f, double G_NFW[2], error **err);
double G_TJ_kappa(double u_sq, double c, double c_sq, error **err);
double G_TJ_gamma(double u_sq, double c, double c_sq, double f, error **err);
void G_TJ_both(double u_sq, double c, double c_sq, double f, double G_TJ[2], error **err);
double G_BMO_kappa(double u_sq, double tau, double tau_sq, error **err);
double kappa_TJ(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err);
void gamma_TJ(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err);
void both_TJ(cosmo_hm *cmhm, halo_t *h, gal_t *g, double theta_sq, error **err);

//-- Functions related to drawing a profile
void fillOneHaloTerm(cosmo_hm *cmhm, halo_t *h, gal_t *g, double_mat *profile, error **err);
double factorForTwoHaloTerm(cosmo_hm *cmhm, halo_t *h, gal_t *g, error **err);
double integrandForTwoHaloTerm(double l, void *inteParam, error **err);
double twoHaloTerm(cosmo_hm *cmhm, halo_t *h, double theta, double factor, error **err);
void fillTwoHaloTerm(cosmo_hm *cmhm, halo_t *h, double factor, double_mat *profile, error **err);
void outputProfile(char name[], cosmo_hm *cmhm, peak_param *peak, halo_t *h, double_mat *profile, error **err);

//-- Functions related to lensing
void lensingForPair(cosmo_hm *cmhm, halo_t *h, gal_t *g, int doKappa, error **err);
void lensingForHalo(cosmo_hm *cmhm, halo_t *h, gal_map *gMap, int doKappa, error **err);
void lensingForMap(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err);

//-- Functions related to mask
short_mat *initializeMask(peak_param *peak, error **err);
short_mat *initializeMask_CFHTLenS_W1(peak_param *peak, error **err);
short_mat *initializeMask_CFHTLenS_W3(peak_param *peak, error **err);
short inCCDMask(short_mat *CCDMask, double pos[2], double theta_CCD_inv[2]);

//-- Functions related to making galaxies
void makeRegularGalaxies(cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void setGalaxySampler(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, error **err);
void fillGalaxyLaw(redshiftANDint *rANDi, sampler_t *samp, error **err);
void makeRandomGalaxies(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, gal_map *gMap, error **err);
void makeMaskedRandomGalaxies(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask, error **err);
void cleanOrMakeOrResample(cosmo_hm *cmhm, peak_param *peak, sampler_t *galSamp, gal_map *gMap, short_mat *CCDMask, error **err);

//-- Functions related to galaxy catalogue
void subtractMean(peak_param *peak, gal_map *gMap);
void makeG(gal_map *gMap);
void addNoiseToGalaxies(peak_param *peak, gal_map *gMap);
void lensingCatalogue(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err);
void lensingCatalogueAndOutputAll(cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err);
void outputGalaxies(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap);

//-- Main function
void doRayTracing(char fileName[], cosmo_hm *cmhm, peak_param *peak, int doNoise, error **err);
void doProfile(char fileName[], cosmo_hm *cmhm, peak_param *peak, double z_l, double M, double z_s,  error **err);


//-- New function
void lensingCatalogueAndOutputAll2(char fileName[],cosmo_hm *cmhm, peak_param *peak, const halo_map *hMap, gal_map *gMap, error **err);


void read_gal_map2(char name[], cosmo_hm *cmhm, peak_param *peak, gal_map *gMap, error **err);
void appendWithSignal_gal_map2(cosmo_hm *cmhm, gal_map *gMap, double z, double pos[2], error **err);
void appendWithSignal_gal_list2(cosmo_hm *cmhm, gal_list *gList, double z, double pos[2], error **err);
void setWithSignal_gal_t2(cosmo_hm *cmhm, gal_t *g, double z, double pos[2], error **err);


void outputFastSimul_galaxies(char name_cmhm[], char name[], char name2[], cosmo_hm *cmhm, peak_param *peak, halo_map *hMap);
void output_halo_map_galaxies(FILE *file, FILE *file2, cosmo_hm *cmhm, peak_param *peak, halo_map *hMap, gal_list *gList);


double NFW(double x);
#endif

