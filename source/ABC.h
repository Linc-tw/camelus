

  /*******************************
   **  ABC.h			**
   **  Chieh-An Lin		**
   **  Version 2015.04.05	**
   *******************************/


#ifndef __ABC__
#define __ABC__

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <fitsio.h>

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"


typedef enum {abd6=0,  pct6=1,  cut6=2,
              abd5=3,  pct5=4,  cut5=5,
                       pct4=7,
              abd2=12, pct2=13, cut2=14,
	      
  pct998=20, pct996=21,
  cut900=30
} summary_t;
#define NB_SUMMARY_T 32
#define STR_SUMMARY_T(i) ( \
  i==0  ? "abd6" : \
  i==1  ? "pct6" : \
  i==2  ? "cut6" : \
  i==3  ? "abd5" : \
  i==4  ? "pct5" : \
  i==5  ? "cut5" : \
  i==7  ? "pct4" : \
  i==12 ? "abd2" : \
  i==13 ? "pct2" : \
  i==14 ? "cut2" : \
  i==20 ? "pct998" : \
  i==21 ? "pct996" : \
  i==30 ? "cut900" : \
  "")

typedef struct {
  int d;
  double *param;
  gsl_vector *param_gsl;
  double diff;
  double weight;
} particle_t;

typedef struct {
  int d;
  int p;
  int nbAttempts;
  double epsilon;
  double_arr *mean;
  gsl_matrix *cov;
  gsl_matrix *cov2;
  gsl_matrix *invCov;
  gsl_permutation *perm;
  particle_t **array;
} particle_arr;

typedef int prior_fct(particle_t*);
typedef void summary_fct(double_arr*, hist_t*, double*);
typedef double distance_fct(double*, double*);

typedef struct {
  //-- Set to customized values
  int d;          //-- Parameter dimension
  int p;          //-- Number of particles
  double r_stop;  //-- Shutoff success rate
  summary_t summ; //-- Summary type
  
  //-- Set to default values
  double epsilon_0; //-- Initial accept tolerance
  prior_fct *priorFct;
  
  //-- Iteration-dependent elements
  int t; //-- Index of iteration
  particle_arr *oldPart;
  particle_arr *newPart;
  double_arr *diffList;
  
  //-- Summary-dependent
  double_arr *obsSummary;
  double_arr *simulSummary;
  hist_t *peakHist;
  summary_fct *summaryFct;
  distance_fct *distFct;
  
  //-- Camelus pipeline
  sampler_arr *sampArr;
  halo_map *hMap;
  gal_map *gMap;
  map_t *kMap, *nMap;
  FFT_t *transformer;
  double_arr *peakList;
} SMC_ABC_t;


//-- Functions related to particle_t
particle_t *initialize_particle_t(int d, error **err);
void free_particle_t(particle_t *pa);
void print_particle_t(particle_t *pa);

//-- Functions related to particle_arr
particle_arr *initialize_particle_arr(int d, int p, error **err);
void free_particle_arr(particle_arr *part);
void print_particle_arr(particle_arr *part);
void output1_particle_arr(FILE *file, particle_arr *part);
void output2_particle_arr(FILE *file, particle_arr *part);
void updateEpsilon_particle_arr(particle_arr *part, double *diffArr);
void updateMean_particle_arr(particle_arr *part);
void updateCovariance_particle_arr(particle_arr *part);
void read_particle_arr(char name[], particle_arr *part, double *diffArr, error **err);

//-- Functions related to SMC_ABC_t
SMC_ABC_t *initialize_SMC_ABC_t(peak_param *peak, error **err);
void free_SMC_ABC_t(SMC_ABC_t *ABC);
void fillObservation_SMC_ABC_t(char name[], SMC_ABC_t *ABC, error **err);

//-- Functions related to accept-reject algorithm
void generateParam(peak_param *peak, particle_arr *oldPart, particle_t *newPa, prior_fct *prior);
cosmo_hm *initialize_cosmo_hm_ABC(double Omega_M, double sigma_8, error **err);
void generateObs(peak_param *peak, SMC_ABC_t *ABC, particle_t *newPa, error **err);
void acceptParticleFromPrior(peak_param *peak, SMC_ABC_t *ABC, particle_t *newPa, error **err);
void acceptParticle(peak_param *peak, SMC_ABC_t *ABC, particle_t *newPa, error **err);

//-- Functions related to loop
void swapOldAndNew(SMC_ABC_t *ABC);
double argOfExp(particle_t *oldPa, particle_t *newPa, gsl_matrix *invCov);
void setWeights(SMC_ABC_t *ABC);
void loopZero(peak_param *peak, SMC_ABC_t *ABC, error **err);
void loop(peak_param *peak, SMC_ABC_t *ABC, error **err);
void readAndLoop(char name[], peak_param *peak, SMC_ABC_t *ABC, error **err);

//-- Functions related to prior
void priorGenerator(gsl_rng *generator, particle_t *pa);
int prior_rectangle(particle_t *pa);
int prior_pentagon(particle_t *pa);

//-- Functions related to summary statistics
void summary_abd_all(double_arr *peakList, hist_t *peakHist, double *summ);
void summary_pct6(double_arr *peakList, hist_t *peakHist, double *summ);
void summary_cut6(double_arr *peakList, hist_t *hist, double *summ);
void summary_pct5(double_arr *peakList, hist_t *hist, double *summ);
void summary_cut5(double_arr *peakList, hist_t *hist, double *summ);
void summary_pct4(double_arr *peakList, hist_t *hist, double *summ);
void summary_pct998(double_arr *peakList, hist_t *hist, double *summ);
void summary_pct996(double_arr *peakList, hist_t *peakHist, double *summ);
void summary_cut900(double_arr *peakList, hist_t *hist, double *summ);

//-- Functions related to distance
double dist_abd6(double *a, double *b);
double dist_abd5(double *a, double *b);
double dist_pct5(double *a, double *b);
double dist_cut5(double *a, double *b);
double dist_6D(double *a, double *b);
double dist_5D(double *a, double *b);
double dist_4D(double *a, double *b);
double dist_1D(double *a, double *b);

//-- Main functions
void doABC(cosmo_hm *cmhm, peak_param *peak, error **err);


#endif


