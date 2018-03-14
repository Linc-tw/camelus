

  /*******************************************************
   **  galaxySampling.h					**
   **  Version 2018.03.12				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_GALAXY_SAMPLING__
#define __CAMELUS_GALAXY_SAMPLING__

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#include "parameters.h"
#include "haloSampling.h"


typedef struct {
  double pos[2];       //-- [rad] Angular position
  double z;            //-- [-] Galaxy redshift
  double a;            //-- [-] 1 / (1+z)
  double w;            //-- [Mpc/h] Comoving radial distance to galaxy
  double D_s;          //-- [Mpc/h] Angular diameter distance
  double weight;       //-- [-] Weight from shape measurement
  double kappa;        //-- [-] Convergence or multiplicative correction
  double gamma[2];     //-- [-] Shear or reduced shear
  double sinCosDEC[2]; //-- [-] sinDEC & cosDEC, only used if field > 0
} gal_t;

typedef struct gal_node {
  gal_t *g;
  struct gal_node *next;
} gal_node;

typedef struct {
  int length;       //-- Number of nodes allocated
  int size;         //-- Number of nodes containing information
  double totWeight; //-- Sum of weights
  double mean[2];   //-- Mean of the signal
  gal_node *first;  //-- Begin of the list
  gal_node *last;   //-- End of the list
} gal_list;

typedef struct {
  int N1, N2;              //-- Number of pixels of each side
  int length;              //-- Number of pixels
  int total;               //-- Number of galaxies
  double theta_pix;        //-- [rad] Pixel size
  double theta_pix_inv;    //-- [rad^-1] Inverse of pixel size
  double kappa_mean;       //-- Mean of kappa over galaxies
  double fillingThreshold; //-- Threshold from filling factor for masking
  map_t type;              //-- Map type (defined in peakParameters.h)
  gal_list **map;          //-- Binned galaxies
} gal_map;

typedef struct {
  int N1, N2;               //-- Number of pixels of each side
  int length;               //-- Number of pixels
  int nbChar;               //-- Number of char, smallest integer larger than number of pixels divided by 8
  double theta_mask[2];     //-- [rad] Pixel size
  double theta_mask_inv[2]; //-- [rad^-1] Inverse of pixel size
  double offset[2];         //-- [pix] Offset of the left bottom corner of the field compared to the mask, always positive or 0, only used for proj_mask
  double center[2];         //-- [deg] (RA, DEC) for the center of the map, only used for HP_mask
  map_t type;               //-- Map type (defined in peakParameters.h)
  char *map;                //-- Array of char each contains 8 pixels. The pixels are ordered as a flattened Cartesian map: (0, 0), (1, 0), ..., (0, 1), (1, 1), ...
} mask_map;


//-- Functions related to gal_t, gal_node, gal_list
void set_gal_t(cosmo_hm* chPar, peak_param* pkPar, gal_t* g, double pos[2], double z, double weight, error** err);
void setWithLensing_gal_t(peak_param *pkPar, gal_t* g, double pos[2], double z, double weight, double kappa, double gamma[2], error** err);
gal_node *initialize_gal_node(error **err);
gal_list *initialize_gal_list(error **err);
void free_gal_list(gal_list *gList);
void cleanLensing_gal_list(gal_list *gList);
void reset_gal_list(gal_list *gList);
void append_gal_list(cosmo_hm* chPar, peak_param* pkPar, gal_list* gList, double pos[2], double z, double weight, error** err);
void appendWithLensing_gal_list(peak_param *pkPar, gal_list* gList, double pos[2], double z, double weight, double kappa, double gamma[2], error** err);

//-- Functions related to gal_map
gal_map *initialize_gal_map(int N1, int N2, double theta_pix, error **err);
void free_gal_map(gal_map *gMap);
void cleanLensing_gal_map(gal_map *gMap);
void reset_gal_map(gal_map *gMap);
int append_gal_map(cosmo_hm* chPar, peak_param* pkPar, gal_map* gMap, double pos[2], double z, double weight, error** err);
int appendWithLensing_gal_map(peak_param* pkPar, gal_map* gMap, double pos[2], double z, double weight, double kappa, double gamma[2], error** err);
void updateCosmo_gal_map(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err);

//-- Functions related to signal mean
void kappaMean_gal_list(gal_list *gList, double threshold);
void gammaOrGMean_gal_list(gal_list *gList, double threshold, int doSurvey);
void setFillingThreshold(peak_param *pkPar, gal_map *gMap, error **err);
void signalMean_gal_map(peak_param *pkPar, gal_map *gMap, error **err);

//-- Functions related to mask
mask_map *initialize_mask_map(int N1, int N2, double theta_mask[2], double offset[2], double center[2], error **err);
void free_mask_map(mask_map *mask);
mask_map *initializeMask(peak_param *pkPar, error **err);
void resetMask(peak_param *pkPar, mask_map *mask);
void generateHole(peak_param *pkPar, mask_map *mask, double theta, int holeInd, error **err);
void makeRandomMask(peak_param *pkPar, mask_map *mask, error **err);
void readInputMask(FITS_t *fits, peak_param *pkPar, mask_map *mask, int limit[4]);
int inMask(mask_map *mask, double pos[2]);

//-- Functions related to I/O
void outAscii_gal_map(FILE *file, peak_param *pkPar, gal_map *gMap);
#ifdef __CAMELUS_USE_FITS__
void outFits_gal_t(FITS_t *fits, gal_t *g, double factor);
void outFits_gal_map(FITS_t *fits, peak_param *pkPar, gal_map *gMap);
#endif
void read_gal_map(char name[], cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err);

//-- Functions related to making galaxies
void makeRegularGalaxies(cosmo_hm *chPar, peak_param *pkPar, gal_map *gMap, error **err);
void setGalaxySampler(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, error **err);
void fillGalaxyLaw(redshiftANDint *rANDi, sampler_t *gSamp, error **err);
void makeRandomGalaxies(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, gal_map *gMap, error **err);
void makeMaskedRandomGalaxies(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, gal_map *gMap, mask_map *mask, error **err);
void rebinGalaxiesForHEALPix(peak_param *pkPar, gal_map *gMap1, gal_map *gMap2, error **err);
void cleanOrMakeOrResample(cosmo_hm *chPar, peak_param *pkPar, sampler_t *gSamp, gal_map *gMap, mask_map *mask, error **err);

#endif

