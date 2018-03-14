

  /*******************************************************
   **  ABC.h						**
   **  Version 2018.03.06				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_ABC__
#define __CAMELUS_ABC__

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#ifdef __CAMELUS_USE_MPI__
#include <mpi/mpi.h>
#endif

#include "parameters.h"
#include "haloSampling.h"
#include "galaxySampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"
#include "multiscale.h"


#define ABC_TOLERANCE_MAX 1e+15

//-- Constants for the prior, to change if the parameter space is enlarged
#define OMEGA_M_MIN   0.10
#define OMEGA_M_MAX   0.90
#define OMEGA_DE_MIN  0.10
#define OMEGA_DE_MAX  0.90
#define OMEGA_B_MIN   0.0001
#define OMEGA_B_MAX   0.90
#define N_S_MIN       0.01
#define N_S_MAX       3.00
#define H_100_MIN     0.01
#define H_100_MAX     3.00 
#define SIGMA_8_MIN   0.30
#define SIGMA_8_MAX   1.60
#define W0_DE_MIN    -1.80
#define W0_DE_MAX     0.00
#define W1_DE_MIN    -2.00
#define W1_DE_MAX     2.00
#define C_0_MIN       1.00
#define C_0_MAX      20.00
#define BETA_NFW_MIN -1.00
#define BETA_NFW_MAX  1.00

//-- To change if the parameter space is enlarged
typedef enum {
  param_Omega_m=0, param_Omega_de=1, param_Omega_b=2, param_n_s=3, param_h_100=4,
  param_sigma_8=5, param_w0_de=6,    param_w1_de=7,   param_c_0=8, param_beta_NFW=9
} param_t;
#define NB_PARAM_T 10
#define STR_PARAM_T(i) ( \
  i==0 ? "Omega_m" : \
  i==1 ? "Omega_de" : \
  i==2 ? "Omega_b" : \
  i==3 ? "n_s" : \
  i==4 ? "h_100" : \
  i==5 ? "sigma_8" : \
  i==6 ? "w0_de" : \
  i==7 ? "w1_de" : \
  i==8 ? "c_0" : \
  i==9 ? "beta_NFW" : \
  "")

typedef enum {
  summ_gauss=0, summ_star=1, summ_tanh=2, summ_mrlens=3, 
  summ_FDR02=4, summ_FDR05=5
} summ_t;
#define NB_SUMM_T 6
#define STR_SUMM_T(i) ( \
  i==0  ? "summ_gauss" : \
  i==1  ? "summ_star" : \
  i==2  ? "summ_tanh" : \
  i==3  ? "summ_mrlens" : \
  i==4  ? "summ_FDR_0.2" : \
  i==5  ? "summ_FDR_0.5" : \
  "")
  
typedef struct {
  int f;                 //-- Dimension of parameter set
  param_t *doParam;      //-- Parameters to include into constraints
                         //-- 0 = Omega_m, 1 = Omega_de, 2 = Omega_b, 3 = n_s, 4 = h_100,
                         //-- 5 = sigma_8, 6 = w0_de,    7 = w1_de,   8 = c_0, 9 = beta_NFW
  int nbParticles;       //-- Number of particles contained in this structure
  int nbAttempts;        //-- Number of attempts
  double epsilon;        //-- Tolerance level for this iteration, median of differences for the last iteration
  gsl_vector *mean;      //-- Mean of parameters
  gsl_vector *buffer;    //-- Buffer vector
  gsl_matrix *cov;       //-- Covariance matrix of parameters
  gsl_matrix *cov2;      //-- Copy of *cov, destroyed in inversion
  gsl_matrix *invCov;    //-- Inverted covariance matrix
  gsl_matrix *cholesky;  //-- Cholesky decomposition used in generating multivariate normal random number
  gsl_permutation *perm; //-- Permutation used in inversion
  double_mat *matrix;    //-- A (f+2)*(Q_MPI*MPISize) matrix for particles, delta, and weight: 
                         //--   pi^0_0,   pi^0_1,   ... , pi^0_f-1,   delta^0,   weight^0,
                         //--   pi^1_0,   pi^1_1,   ... , pi^1_f-1,   delta^1,   weight^1,
                         //--     ...      ...      ...     ...         ...        ...
                         //--   pi^Q-1_0, pi^Q-1_1, ... , pi^Q-1_f-1, delta^Q-1, weight^Q-1
                         //--     ...      ...      ...     ...         ...        ...
} iteration_t;

typedef int prior_fct(double_mat*, double*);
typedef void summary_fct(void*, void*, int, double*);
typedef double distance_fct(double*, double*, gsl_matrix*, gsl_vector*, gsl_vector*);

typedef struct {
  //----------------------------------------------------------------------
  //-- Customized part
  
  int f;                    //-- Dimension of parameter set
  int Q;                    //-- Number of particles
  double r_stop;            //-- Shutoff success rate
  int nbAttempts_max;       //-- Maximum number of attempts
  
  //-- Precomputed
  int Q_MPI;                //-- Number of particles per processor to compute when using MPI
  
  //----------------------------------------------------------------------
  //-- Iteration part
  
  int t;                    //-- Index of iteration
  iteration_t *oldIter;     //-- Old iteration
  iteration_t *newIter;     //-- New iteration
  double_arr *deltaList;    //-- Buffer structure for the median of delta
  
  //----------------------------------------------------------------------
  //-- Prior part
  
  prior_fct *priorFct;      //-- Prior distribution
  double_mat *limit;        //-- Lower and upper limits of the prior
  int flatness;             //-- +1 = Omega_de set to 1 - Omega_b
                            //-- -1 = Omega_b  set to 1 - Omega_de
                            //--  0 = no flatness garanteed
  double value[NB_PARAM_T]; //-- Values for parameters suceptible to be constrained
  
  //----------------------------------------------------------------------
  //-- Summary part
  
  double_arr *x_obs;        //-- Summary statistics from observation
  double_arr *x_mod;        //-- Summary statistics from model
  summary_fct *summFct;     //-- Summary statistics function
  distance_fct *distFct;    //-- Distance function
  gsl_matrix *invCov;
  gsl_vector *Delta_x;
  gsl_vector *intermediate;
  
  //----------------------------------------------------------------------
} PMC_ABC_t;


//-- Functions related to print
void valueStringForPrint(char buffer2[], param_t p, double value);
void paramStringForOutput(char buffer2[], param_t p);
void valueStringForOutput(char buffer2[], param_t p, double value);
void particleString(char buffer1[], char buffer2[], int f, param_t *doParam, double *part);
void completeParticleString(char buffer1[], char buffer2[], double *value);

//-- Functions related to iteration_t
iteration_t *initialize_iteration_t(int f, int *doParam, int Q, error **err);
void free_iteration_t(iteration_t *iter);
void print_gsl_matrix(gsl_matrix *mat);
void print_iteration_t(iteration_t *iter);
void output_gsl_matrix(FILE *file, gsl_matrix *mat);
void output_iteration_t(FILE *file, iteration_t *iter);
void updateMean_iteration_t(iteration_t *iter);
void updateCovariance_iteration_t(iteration_t *iter);
void updateCholesky_iteration_t(iteration_t *iter);

//-- Functions related to PMC_ABC_t;
PMC_ABC_t *initialize_PMC_ABC_t(cosmo_hm *chPar, peak_param *pkPar, error **err);
void free_PMC_ABC_t(PMC_ABC_t *ABC);
void print_PMC_ABC_t(peak_param *pkPar, PMC_ABC_t *ABC);
int garanteeFlatness(int f, param_t *doParam);
void fillLimitAndValue(cosmo_hm *chPar, PMC_ABC_t *ABC);
void readObservation(char name[], peak_param *pkPar, PMC_ABC_t *ABC, error **err);
void readIteration(char name[], peak_param *pkPar, PMC_ABC_t *ABC, int t, error **err);

//-- Functions related to accept-reject algorithm
void generateParam(PMC_ABC_t *ABC, gsl_rng *generator, double *newPart, prior_fct *prior);
cosmo_hm *initialize_cosmo_hm_ABC(cosmo_hm *oldChPar, PMC_ABC_t *ABC, double *newPart, error **err);
void generateModel(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, double *newPart, error **err);
int oneSampleTest(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, double *newPart, error **err);

//-- Functions related to iteration
void swapOldAndNew(PMC_ABC_t *ABC);
double chiSquared(double *x1, double *x2, gsl_matrix *invCov, gsl_vector *Delta_x, gsl_vector *intermediate);
void setWeights(PMC_ABC_t *ABC);
void setEpsilon(PMC_ABC_t *ABC);
void runIteration(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, error **err);
void outputIteration(char name[], PMC_ABC_t *ABC);

//-- Functions related to subset
void runIterationWithoutUpdate(cosmo_hm *chPar, peak_param *pkPar, pipeline_t *pipe, PMC_ABC_t *ABC, error **err);
void readSubset(char name[], iteration_t *iter, error **err);
void print_PMC_ABC_t_subset(peak_param *pkPar, PMC_ABC_t *ABC);

//-- Functions related to prior
void priorGenerator(gsl_rng *generator, double_mat *limit, double *part);
int prior_rectangle(double_mat *limit, double *part);
int prior_pentagon(double_mat *limit, double *part);

//-- Functions related to summary statistics
void summary_int(void *arg0, void *arg1, int length, double *summ);
void summary_double(void *arg0, void *arg1, int length, double *summ);

//-- Functions related to distance
double dist_chi(double *x1, double *x2, gsl_matrix *invCov, gsl_vector *Delta_x, gsl_vector *intermediate);
double dist_uncorrChi(double *x1, double *x2, gsl_matrix *invCov, gsl_vector *Delta_x, gsl_vector *intermediate);
void readInvCov(char name[], peak_param *pkPar, gsl_matrix *invCov, error **err);

//-- Main functions
void doABC(cosmo_hm *chPar, peak_param *pkPar, error **err);
void doABC_subset(cosmo_hm *chPar, peak_param *pkPar, int t, int nbAttempts, int subsetInd, error **err);
void doABC_gather(cosmo_hm *chPar, peak_param *pkPar, int t, int nbSubsets, error **err);

#endif

