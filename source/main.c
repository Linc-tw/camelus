

  /*******************************************************
   **  main.c						**
   **  Version 2018.03.15				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "main.h"


int main(int argc, char *argv[])
{
  int MPISize = 1;
  int MPIInd  = 0;
  
#ifdef __CAMELUS_USE_MPI__
  //-- MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIInd);
#endif
  
  //-- Start stopwatch
  clock_t start = clock();
  if (MPIInd == 0) {
    printf("\n");
    printf("        *************************************************\n");
    printf("        **                Camelus v2.0                 **\n");
    printf("        **                                             **\n");
    printf("        **  Copyright (C) 2018 - Chieh-An Lin          **\n");
    printf("        **  GNU GPLv3 - https://www.gnu.org/licenses/  **\n");
    printf("        *************************************************\n");
    printf("\n");
  }
  
  //-- Read inputs
  char *pkParPath = (argc >= 2) ? argv[1] : "../param/peakParam.par";
  if (!strcmp(pkParPath, "default"))     pkParPath = "../param/peakParam.par";
  else if (!strcmp(pkParPath, "LHF"))    pkParPath = "../param/peakParam_LHF.par";
  else if (!strcmp(pkParPath, "peakVI")) pkParPath = "../param/peakParam_peakVI.par";
  int task = (argc >= 3) ? atoi(argv[2]) : -1;
  
  //-- Initialization
  error *theError = NULL;
  error **err = &theError;
  peak_param *pkPar = initialize_peak_param(err);                  quitOnError(*err, __LINE__, stderr);
  cosmo_hm *chPar = initialize_cosmo_hm_default(err);              quitOnError(*err, __LINE__, stderr);
  read_peak_param(pkParPath, pkPar, err);                          quitOnError(*err, __LINE__, stderr); //-- Read parameters from .par
  read_cosmo_hm(pkPar->hmParPath, chPar, pkPar, err);              quitOnError(*err, __LINE__, stderr); //-- Read parameters from .par
  int help = updateFromCommandLine(argc, argv, chPar, pkPar, err); quitOnError(*err, __LINE__, stderr); //-- Update parameters from the command line
  checkParam(pkPar, err);                                          quitOnError(*err, __LINE__, stderr); //-- Check if some are missing
  chPar = reinitialize_cosmo_hm(chPar, err);                       quitOnError(*err, __LINE__, stderr); //-- Precalculate some parameters
  set_peak_param(chPar, pkPar, err);                               quitOnError(*err, __LINE__, stderr); //-- Precalculate some parameters
  pkPar->MPISize = MPISize;
  pkPar->MPIInd  = MPIInd;
  
  //-- Print, only for MPI root processor
  if (MPIInd == 0) {
    printf("Initialized Camelus\n");
    printf("------------------------------------------------------------------------\n");
      
    if (!strcmp(pkParPath, "../param/peakParam.par")) {
      printParam(chPar, pkPar);
//       printParam_peakCustomized_1(pkPar);
//       printf("\n");
//       printParam_peakCustomized_2(pkPar);
//       printf("\n");
//       printParam_peakPrecomputed(pkPar);
      printf("------------------------------------------------------------------------\n");
    }
#ifdef __CAMELUS_USE_LHF__
    else if (!strcmp(pkParPath, "../param/peakParam_LHF.par")) {
      _LHF__printParam(pkPar);
      printf("------------------------------------------------------------------------\n");
    }
#endif
#ifdef __CAMELUS_BETA_MODE__
    else if (!strcmp(pkParPath, "../param/peakParam_peakVI.par")) {
      _peakVI__printParam(chPar, pkPar);
      printf("------------------------------------------------------------------------\n");
    }
#endif
  }
  
#ifdef __CAMELUS_USE_MPI__
  //-- Stopping lock for no MPI task
  if (MPISize > 1 && !(task == 0 || task == 8)) {
    printf("Cannot use multiple processors for task %d\n", task);
    MPI_Finalize();
    return 1;
  }
  
  //-- Synchronize all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  //-- Sandbox
  if (task == 0) {
    sandbox(chPar, pkPar, err);
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Mass function
  else if (task == 1) {
    if (argc < 4 || help) printInstructions(task, 1);
    else {
      double z = atof(argv[3]);
      pkPar->verbose = (pkPar->verbose == 0) ? 1 : pkPar->verbose;
      doMassFct(chPar, pkPar, z, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- Mass sheet
  else if (task == 2) {
    if (argc < 4 || help) printInstructions(task, 1);
    else {
      double z_s = atof(argv[3]);
      doMassSheet(chPar, pkPar, z_s, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- Profile
  else if (task == 3) {
    if (argc < 6 || help) printInstructions(task, 1);
    else {
      double z_l = atof(argv[3]);
      double M   = atof(argv[4]);
      double z_s = atof(argv[5]);
      doProfile(chPar, pkPar, z_l, M, z_s, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- Fast simulation
  else if (task == 4) {
    if (help) printInstructions(task, 1);
    else {
      pkPar->verbose = (pkPar->verbose == 0) ? 1 : pkPar->verbose;
      if (pkPar->outHaloCat == 0) printf("Found outHaloCat = 0, forced to 1\n");
      pkPar->outHaloCat = 1;
      doFastSimulation(chPar, pkPar, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- Compute lensing maps
  else if (task == 5) {
    if (help) printInstructions(task, 1);
    else {
      pkPar->verbose = (pkPar->verbose == 0) ? 1 : pkPar->verbose;
      doKMap(chPar, pkPar, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- Multiscale data
  else if (task == 6) {
    if (help) printInstructions(task, 1);
    else {
      pkPar->verbose = (pkPar->verbose == 0) ? 1 : pkPar->verbose;
      doMultiscale(chPar, pkPar, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- Data matrix
  else if (task == 7) {
    if (argc < 4 || help) printInstructions(task, 1);
    else {
      pkPar->verbose = (pkPar->verbose == 0) ? 3 : pkPar->verbose;
      int N = atoi(argv[3]);
      doDataMatrix(chPar, pkPar, N, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  //-- ABC
  else if (task == 8) {
    if (!(argc == 3 || argc == 5 || argc == 6) || help) printInstructions(task, 1);
    else {
      int t          = (argc >= 4) ? atoi(argv[3]) : 0;
      int nbAttempts = (argc >= 5) ? atoi(argv[4]) : 0;
      int subsetInd  = (argc >= 6) ? atoi(argv[5]) : 0;
      if      (argc == 3) doABC(chPar, pkPar, err); 
      else if (argc == 6) doABC_subset(chPar, pkPar, t, nbAttempts, subsetInd, err);
      else                doABC_gather(chPar, pkPar, t, nbAttempts, err);
      quitOnError(*err, __LINE__, stderr);
    }
  }

#ifdef __CAMELUS_USE_LHF__
  else if (task == 51) {
    if (!(argc == 5) || help) printInstructions(task, 1);
    else {
      pkPar->verbose = (pkPar->verbose == 0) ? 1 : pkPar->verbose;
      int N1 = atoi(argv[3]);
      int N2 = atoi(argv[4]);
      _LHF__doTrainingSet(N1, N2);
    }
  }
#endif
  
#ifdef __CAMELUS_BETA_MODE__
  else if (task == 61) {
    if (!(argc == 4) || help) printInstructions(task, 1);
    else {
      pkPar->verbose = (pkPar->verbose == 0) ? 3 : pkPar->verbose;
      _peakVI__doLSSMap_fullSky(chPar, pkPar, argv[3], err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  else if (task == 62) {
    if (!(argc == 5) || help) printInstructions(task, 1);
    else {
      int N1 = atoi(argv[3]);
      int N2 = atoi(argv[4]);
      _peakVI__doHaloAssignment(chPar, pkPar, N1, N2, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  else if (task == 63) {
    if (!(argc == 5) || help) printInstructions(task, 1);
    else {
      int N1 = atoi(argv[3]);
      int N2 = atoi(argv[4]);
      _peakVI__doKMap_regular(chPar, pkPar, N1, N2, err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  else if (task == 64) {
    if (!(argc == 4) || help) printInstructions(task, 1);
    else {
      _peakVI__doRayTracing(chPar, pkPar, argv[3], err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  else if (task == 65) {
    if (!(argc == 4) || help) printInstructions(task, 1);
    else {
      _peakVI__doKMap(chPar, pkPar, argv[3], err); quitOnError(*err, __LINE__, stderr);
    }
  }
  
  else if (task == 66) {
    if (!(argc == 6) || help) printInstructions(task, 1);
    else {
      int N1 = atoi(argv[3]);
      int N2 = atoi(argv[4]);
      int doRandom = atoi(argv[5]);
      _peakVI__doMultiscale_clustered(chPar, pkPar, N1, N2, doRandom, err); quitOnError(*err, __LINE__, stderr);
    }
  }
#endif
  
  else if (task == -1) printInstructions(-1, 1);
  else if (task == -2) printInstructions(-2, 0);
  else                 printInstructions(-1, 0);
  
  free_parameters_hm(&chPar);
  free_peak_param(pkPar);
  
  //-- Stop stopwatch
  clock_t finish = clock();
  if (MPIInd == 0) {
    printTime(start, finish);
    printf("\n");
  }
#ifdef __CAMELUS_USE_MPI__
  MPI_Finalize();
#endif
  return 0;
}

void printInstructions(int task, int doHelp)
{
  if (task >= 1 && task <= 30) {
    printDetails(1, 1, 0);
    printDetails(task, 0, doHelp);
  }
  
#ifdef __CAMELUS_USE_LHF__
  else if (task >= 51 && task <= 60) {
    printDetails(51, 1, 0);
    printDetails(51, 0, 0);
  }
#endif
  
#ifdef __CAMELUS_BETA_MODE__
  else if (task >= 61 && task <= 70) {
    printDetails(61, 1, 0);
    printDetails(61, 0, 0);
  }
#endif

  else if (task == -2) { //-- Simple menu
    printDetails(1, 1, doHelp);
    printDetails(1, 0, 0); printDetails(4, 0, 0); printDetails(5, 0, 0);
    printDetails(6, 0, 0); printDetails(7, 0, 0); printDetails(8, 0, 0);
#ifdef __CAMELUS_USE_LHF__
    printf("\n");
    printDetails(51, 1, 0);
    printDetails(51, 0, 0);
#endif
#ifdef __CAMELUS_BETA_MODE__
    printf("\n");
    printDetails(61, 1, 0);
    printDetails(61, 0, 0);
#endif
  }
  
  else { //-- Complete menu
    printDetails(1, 1, doHelp);
    printDetails(1, 0, 0); printDetails(2, 0, 0); printDetails(3, 0, 0); printDetails(4, 0, 0); printDetails(5, 0, 0); 
    printDetails(6, 0, 0); printDetails(7, 0, 0); printDetails(8, 0, 0);
#ifdef __CAMELUS_USE_LHF__
    printf("\n");
    printDetails(51, 1, 0);
    printDetails(51, 0, 0);
#endif
#ifdef __CAMELUS_BETA_MODE__
    printf("\n");
    printDetails(61, 1, 0);
    printDetails(61, 0, 0);
#endif
  }
  
  printf("------------------------------------------------------------------------\n");
  return;
}

void printDetails(int task, int doHeader, int doHelp)
{
  if (doHeader) {
    if (task >= 1 && task <= 30) {
      if (doHelp) {
      printf("Commands:\n");
      printf("  ./camelus PATH TASK\n");
      printf("  ./camelus PATH TASK -h\n");
      printf("  ./camelus PATH TASK KEY=VALUE\n");
      printf("\n");
      printf("Notes:\n");
      printf("  - PATH = path to the peakParam.par file\n");
      printf("  - Can replace PATH with the string \"default\", equivalent to \"../param/peakParam.par\"\n");
      printf("  - TASK = number of the task to do\n");
      printf("  - Add -h after TASK for more detailed instructions\n");
      printf("  - For some tasks, parameters can be updated by KEY=VALUE.\n");
      printf("\n");
      printf("Examples:\n");
      printf("  ./camelus default 1\n");
      printf("  ./camelus ../param/peakParam.par 4\n");
      printf("  ./camelus default 7 sigma_8=0.82 outMaps=1\n");
      printf("\n");
      }
      printf("Tasks:\n");
    }
    
    else if (task >= 51 && task <= 60) printf("For LHF:\n");
    else if (task >= 61 && task <= 70) printf("For peakVI:\n");
    return; 
  }
  
  if (task == 1) {
      printf("  1 = Compute mass function\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 1  Z    # Mass function at redshift Z\n");
      printf("  ./camelus PATH 1 -1    # Output more than one redshift at a time (interactive mode)\n");
    }
  }
  else if (task == 2) {
      printf("  2 = Mass sheet values for convergence\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 2 ZS    # Compute kappa_1 at redshift ZS\n");
      printf("\n");
      printf("Notes:\n");
      printf("  - TODO\n");
    }
  }
  else if (task == 3) {
      printf("  3 = Halo lensing profile with 1-h and 2-h terms\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 3 ZL M ZS    # ZL = lens redshift, ZS = source redshift, M = lens mass in [M_sol/h]\n");
    }
  }
  else if (task == 4) {
      printf("  4 = Fast simulation\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 4                # Create a halo catalogue by sampling a mass function\n");
      printf("  ./camelus PATH 4 [KEY=VALUE]    # Dynamically update parmeters\n");
    }
  }
  else if (task == 5) {
      printf("  5 = Lensing maps\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 5                # Generate differnt map products\n");
      printf("  ./camelus PATH 5 [KEY=VALUE]    # Dynamically update parmeters\n");
    }
  }
  else if (task == 6) {
      printf("  6 = One realization of multiscale peak counts\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 6                # Peak histograms from all filters\n");
      printf("  ./camelus PATH 6 [KEY=VALUE]    # Dynamically update parmeters\n");
    }
  }
  else if (task == 7) {
      printf("  7 = Data matrix of several realizations\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH 7 N                # Make N realizations of a field and count multiscale peaks\n");
      printf("  ./camelus PATH 7 N [KEY=VALUE]    # Dynamically update parmeters\n");
      printf("\n");
      printf("Notes:\n");
      printf("  - Cannot use external halo catalogues\n");
      printf("  - Cannot output intermediate files\n");
    }
  }
  else if (task == 8) {
      printf("  8 = ABC analysis\n");
    if (doHelp) {
      printf("\n");
      printf("Commands:\n");
      printf("  ./camelus PATH  8          # Do the entire ABC run\n");
      printf("  ./camelus PATH  8 T N I    # For the I-th subset at iteration T, do N attempts of ABC\n");
      printf("  ./camelus PATH  8 T S      # Gather all subsets from I = 0 to S-1 at iteration T\n");
      printf("\n");
      printf("Notes:\n");
      printf("  - Not safe to do the entire run\n");
    }
  }
  else if (task >= 51 && task <= 60) {
      printf("  ./camelus default  51 N1 N2    # Do training set\n");
  }
  else if (task >= 61 && task <= 70) {
      printf("  ./camelus peakVI  61 str             # LSS maps\n");
      printf("  ./camelus peakVI  62 N1 N2           # Assign halos\n");
      printf("  ./camelus peakVI  63 N1 N2           # Noise-free lensing map at z_s = 1\n");
      printf("  ./camelus peakVI  64 str             # Full curved-sky RT\n");
      printf("  ./camelus peakVI  65 N1 N2           # Full curved-sky filtering\n");
      printf("  ./camelus peakVI  66 N1 N2 doRand    # Fast peak counts\n");
  }
  return;
}

void sandbox(cosmo_hm *chPar, peak_param *pkPar, error **err)
{
#ifdef __CAMELUS_BETA_MODE__
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  printf("This is a sandbox of Camelus.\n\n");
  
  
  
  
//   for (i=0; i<4; i++) {
//     printf("pix[%d] = %ld\n", i, pix[i]);
//   }
  
  
  /*
  double nbHalos[9] = {         839,         1346,          851,           1052,         1681,         1010,            812,         1158,          665};
  double zLArr[9]   = {    0.195620,     0.305570,     0.374820,       0.192260,     0.303390,     0.374685,       0.187840,     0.302895,     0.376040};
  double logMArr[9] = {   14.046000,    14.047000,    14.043000,      14.185000,    14.182000,    14.183000,      14.438500,    14.430000,    14.438000};
  double MArr[9]    = {1.111731e+14, 1.114294e+14, 1.104079e+14,   1.531089e+14, 1.520548e+14, 1.524051e+14,   2.744732e+14, 2.691537e+14, 2.741572e+14};
  double z_s = 1.1689;
  
  int ind = 6;
  for (ind=0; ind<9; ind++) {
    sprintf(name, "profile_z%.4f_hBin%d", z_s, ind);
    doProfile(name, chPar, pkPar, zLArr[ind], MArr[ind], z_s, err); forwardError(*err, __LINE__,);
  }
  //sprintf(name, "profile_BMO_z%.4f", z_s);
  //doProfile(name, chPar, pkPar, zLArr[ind], MArr[ind], z_s, err); forwardError(*err, __LINE__,);
  */
  
  //mpirun -n 3 ./camelus default 35
  //-- Hold root processor if others have not finished
#endif
  printf("------------------------------------------------------------------------\n");
  return;
}

