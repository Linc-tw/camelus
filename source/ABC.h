

  /*******************************
   **  ABC.h			**
   **  Chieh-An Lin		**
   **  Version 2015.12.09	**
   *******************************/


#ifndef __ABC__
#define __ABC__

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <mpi.h>

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"


typedef enum {
  summ_gauss=0, summ_star=1, summ_mrlens=2
} summary_t;
#define NB_SUMMARY_T 3
#define STR_SUMMARY_T(i) ( \
  i==0  ? "summ_gauss" : \
  i==1  ? "summ_star" : \
  i==2  ? "summ_mrlens" : \
  "")

typedef struct {
  int f;                 //-- Parameter dimension
  int Q;                 //-- Number of particles
  int Q_MPI;             //-- Number of particles per processor to compute when using MPI
  int nbAttempts;        //-- Number of attempts
  double epsilon;        //-- Tolerance level for next iteration, median of differences for this iteration
  gsl_vector *mean;      //-- Mean of parameters
  gsl_vector *buffer;    //-- Buffer vector
  gsl_matrix *cov;       //-- Covariance matrix of parameters
  gsl_matrix *cov2;      //-- Copy of *cov, destroyed in inversion
  gsl_matrix *invCov;    //-- Inverted covariance matrix
  gsl_matrix *cholesky;  //-- Cholesky decomposition used in generating multivariate normal random number
  gsl_permutation *perm; //-- Permutation used in inversion
  double_mat *matrix;    //-- A (d+2)*(p_MPI*MPISize) matrix for particles, delta, and weight: 
                         //--   pi^0_0,   pi^0_1,   ... , pi^0_N-1,   delta_0,   weight_0,
                         //--   pi^1_0,   pi^1_1,   ... , pi^1_N-1,   delta_1,   weight_1,
                         //--     ...      ...      ...     ...         ...        ...
                         //--   pi^p-1_0, pi^p-1_1, ... , pi^p-1_N-1, delta_p-1, wieght_p-1
                         //--     ...      ...      ...     ...         ...        ...
} iteration_t;

typedef int prior_fct(double*);
typedef void summary_fct(double_arr*, hist_t*, double_mat*, double*);
typedef double distance_fct(double*, double*);

typedef struct {
  //----------------------------------------------------------------------
  //-- Customized part
  
  int Q;                   //-- Number of particles
  double r_stop;           //-- Shutoff success rate
  summary_t summ;          //-- Summary type
  
  //----------------------------------------------------------------------
  //-- Default part
  
  int f;                   //-- Parameter dimension
  prior_fct *priorFct;     //-- Prior distribution
  
  //----------------------------------------------------------------------
  //-- Precomputed part
  
  int Q_MPI;               //-- Number of particles per processor to compute when using MPI
  
  //----------------------------------------------------------------------
  //-- Iteration-dependent elements
  
  int t;                   //-- Index of iteration
  iteration_t *oldIter;    //-- Old iteration
  iteration_t *newIter;    //-- New iteration
  double_arr *deltaList;   //-- Buffer structure for the median of delta
  
  //----------------------------------------------------------------------
  //-- Summary-dependent
  
  double_arr *x_obs;     //-- Summary statistics from observation
  double_arr *x_mod;     //-- Summary statistics from model
  summary_fct *summFct;  //-- Summary statistics function
  distance_fct *distFct; //-- Distance function
  
  //----------------------------------------------------------------------
  //-- Camelus pipeline
  
  sampler_arr *sampArr;
  halo_map *hMap;
  sampler_t *galSamp;
  gal_map *gMap;
  short_mat *CCDMask;
  FFT_arr *smoother;
  map_t *kMap;
  FFT_arr *variance;
  double_arr *peakList;
  hist_t *hist;
  hist_t *hist2;
  double_mat *multiscale;
  
  //----------------------------------------------------------------------
} PMC_ABC_t;


//-- Functions related to iteration_t
iteration_t *initialize_iteration_t(int f, int Q, int MPISize, error **err);
void free_iteration_t(iteration_t *iter);
void print_iteration_t(iteration_t *iter);
void output_gsl_matrix(FILE *file, gsl_matrix *mat);
void output1_iteration_t(FILE *file, iteration_t *iter);
void output2_iteration_t(FILE *file, iteration_t *iter);
void updateMean_iteration_t(iteration_t *iter);
void updateCovariance_iteration_t(iteration_t *iter);
void updateCholesky_iteration_t(iteration_t *iter);
void read_iteration_t(char name[], iteration_t *iter, double *deltaArr, error **err);

//-- Functions related to PMC_ABC_t;
PMC_ABC_t *initialize_PMC_ABC_t(peak_param *peak, error **err);
void free_PMC_ABC_t(PMC_ABC_t *ABC);
void fillObservation(char name[], peak_param *peak, PMC_ABC_t *ABC, error **err);

//-- Functions related to accept-reject algorithm
void generateParam(peak_param *peak, iteration_t *oldIter, double *newPart, prior_fct *prior);
cosmo_hm *initialize_cosmo_hm_ABC(double Omega_M, double sigma_8, double w0_de, error **err);
void generateModel(peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err);
void acceptParticleFromPrior(peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err);
void acceptParticle(peak_param *peak, PMC_ABC_t *ABC, double *newPart, error **err);

//-- Functions related to loop
void swapOldAndNew(PMC_ABC_t *ABC);
double argOfExp(double *oldPa, double *newPart, gsl_matrix *invCov, gsl_vector *Delta_param, gsl_vector *intermediate);
void setWeights(PMC_ABC_t *ABC);
void setEpsilon(PMC_ABC_t *ABC);
void loopZero(peak_param *peak, PMC_ABC_t *ABC, error **err);
void loop(peak_param *peak, PMC_ABC_t *ABC, error **err);
void readAndLoop(char name[], peak_param *peak, PMC_ABC_t *ABC, error **err);

//-- Functions related to prior
void priorGenerator(gsl_rng *generator, double *part);
int prior_rectangle(double *part);
int prior_pentagon(double *part);

//-- Functions related to summary statistics
void summary_multiscale(double_arr *peakList, hist_t *hist, double_mat *multiscale, double *summ);

//-- Functions related to distance
double dist_gauss(double *a, double *b);
double dist_star(double *a, double *b);
double dist_6D(double *a, double *b);
double dist_5D(double *a, double *b);
double dist_4D(double *a, double *b);
double dist_1D(double *a, double *b);

//-- Main functions
void doABC(cosmo_hm *cmhm, peak_param *peak, error **err);


#endif

