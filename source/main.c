

  /*******************************
   **  main.c			**
   **  Chieh-An Lin		**
   **  Version 2015.12.09	**
   *******************************/


#include "main.h"


#ifdef __releaseMenu__
int main(int argc, char *argv[])
{
  int MPISize, MPIInd;
  
  //-- MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIInd);
  
  //-- Start stopwatch
  clock_t start = clock();
  
  //-- Read inputs
  char *arg2, *arg3, *arg4, *arg5;
  int task;
  if (argc>=2) task = atoi(argv[1]);
  if (argc>=3) arg2 = argv[2];
  if (argc>=4) arg3 = argv[3];
  if (argc>=5) arg4 = argv[4];
  if (argc>=6) arg5 = argv[5];
  
  //-- Stopping lock for no MPI task
  if (MPISize > 1 && !(task == -1 || task == 7)) {
    printf("Cannot use multiple processors for task %d\n", task);
    MPI_Finalize();
    return 1;
  }
  
  error *myerr = NULL, **err = &myerr;
  cosmo_hm *cmhm;
  peak_param *peak;
  
  //-- Initialize cosmology and halo model
  read_cosmo_hm("../param/hmParam.par", &cmhm, err);       quitOnError(*err, __LINE__, stderr);
  peak = malloc_err(sizeof(peak_param), err);              quitOnError(*err, __LINE__, stderr);
  read_peak_param("../param/peakParam.par", peak, err);    quitOnError(*err, __LINE__, stderr);
  
  //-- Reinitialization, if necessary
  if (argc == 6 && task == 6) {
    double Omega_m = atof(arg3); 
    double sigma_8 = atof(arg4);
    double w0_de   = atof(arg5);
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, w0_de, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Precalculate some parameters, always put this after reinitialization
  set_peak_param(cmhm, peak, err);                         quitOnError(*err, __LINE__, stderr);
  peak->MPISize = MPISize;
  peak->MPIInd  = MPIInd;
    
  //-- Print, only for MPI root processor
  if (MPIInd == 0) {
    printf("\n");
    printf("                  *********************************\n");
    printf("                  **        Camelus v1.3         **\n");
    printf("                  **  Chieh-An Lin (CEA Saclay)  **\n");
    printf("                  *********************************\n");
    printf("\n");
    printf("Initialization done\n");
    printf("------------------------------------------------------------------------\n");
    if (task == 7) printParam_ABC(cmhm, peak);
    else           printParam(cmhm, peak);
  }
  
  //-- Synchronize all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
    
  //-- Below are different functions of Camelus
  
  //-- Sandbox
  if (task == -1) {sandbox(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);}
  
  //-- Print mass function
  else if (task == 1) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    double z = atof(arg2);
    char name[STRING_LENGTH_MAX];
    sprintf(name, "massFct_z%.3f", z);
    outputMassFct(name, cmhm, peak, z, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fast simulation
  else if (task == 2) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doFastSimulation(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Make maps
  else if (task == 3) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doKMap(NULL, cmhm, peak, 1, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list
  else if (task == 4) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doPeakList(NULL, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Multiscale data
  else if (task == 5) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doMultiscale(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Data matrix
  else if (task == 6) {
    if (!(argc == 3 || argc == 6)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    doDataMatrix(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 7) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  else {
    printInstructions(-1, 1);
    return 1;
  }
    
  free_peak_param(peak);
  free_parameters_hm(&cmhm);
  
  //-- Stop stopwatch
  clock_t finish = clock();
  if (MPIInd == 0) {
    printTime(start, finish);
    printf("\n");
  }
  
  MPI_Finalize();
  return 0;
}

void printInstructions(int task, int printHeader)
{
  if (printHeader) {
    printf("Use Camelus with:\n");
    printInstructions(1, 0); printInstructions(2, 0); printInstructions(3, 0); printInstructions(4, 0); printInstructions(5, 0);
    printInstructions(6, 0); printInstructions(7, 0);
    printf("\n");
  }
  
  else {
    switch (task) {
      case 1:
	printf("  ./camelus 1 z                        # Print halo mass function at z\n");
	break;
      case 2:
	printf("  ./camelus 2                          # Fast simulation\n");
	break;
      case 3:
	printf("  ./camelus 3                          # Lensing map and intermediate products for the first filter\n");
	break;
      case 4:
	printf("  ./camelus 4                          # Peak list and histogram from fast simulation\n");
	break;
      case 5:
	printf("  ./camelus 5                          # A realization of multiscale data\n");
	break;
      case 6:
	printf("  ./camelus 6 N                        # Data matrix of N realizations\n");
	printf("  ./camelus 6 N Omega_m sigma_8 w0_de  # Data matrix of N realizations for a given (Omega_m, sigma_8, w0_de)\n");
	break;
      case 7:
	printf("  ./camelus 7                          # ABC analysis\n");
	break;
    }
  }
  
  return;
}

#else
int main(int argc, char *argv[])
{
  int MPISize, MPIInd;
  
  //-- MPI initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIInd);
  
  //-- Start stopwatch
  clock_t start = clock();
  
  //-- Read inputs
  char *simulName = (argc >= 2) ? argv[1]       : "default";
  int task        = (argc >= 3) ? atoi(argv[2]) : 0;
  char *arg3      = (argc >= 4) ? argv[3]       : NULL;
  char *arg4      = (argc >= 5) ? argv[4]       : NULL;
  char *arg5      = (argc >= 6) ? argv[5]       : NULL;
  char *arg6      = (argc >= 7) ? argv[6]       : NULL;
  
  //-- Stopping lock for no MPI task
  if (MPISize > 1 && !(task == -1 || task == 35)) {
    printf("Cannot use multiple processors for task %d\n", task);
    MPI_Finalize();
    return 1;
  }
  
  error *myerr = NULL, **err = &myerr;
  cosmo_hm *cmhm;
  peak_param *peak;
  
  //-- Initialize cosmology and halo model
  if (!strcmp(simulName, "default")) {
    read_cosmo_hm("../param/hmParam.par", &cmhm, err);              quitOnError(*err, __LINE__, stderr);
    peak = malloc_err(sizeof(peak_param), err);                     quitOnError(*err, __LINE__, stderr);
    read_peak_param("../param/peakParam.par", peak, err);           quitOnError(*err, __LINE__, stderr);
  }
  else if (strstr(simulName, "aardvark") != NULL) {
    _paperI__read_cosmo_hm_aardvark(simulName, &cmhm, err);         quitOnError(*err, __LINE__, stderr);
    peak = _paperI__initialize_peak_param_aardvark(simulName, err); quitOnError(*err, __LINE__, stderr);
  }
  else if (!strcmp(simulName, "paperII")) {
    read_cosmo_hm("../param/hmParam_paperII.par", &cmhm, err);      quitOnError(*err, __LINE__, stderr);
    peak = _paperII__initialize_peak_param(simulName, err);         quitOnError(*err, __LINE__, stderr);
  }
  else if (!strcmp(simulName, "paperIII")) {
    read_cosmo_hm("../param/hmParam_paperIII.par", &cmhm, err);     quitOnError(*err, __LINE__, stderr);
    peak = _paperIII__initialize_peak_param(simulName, err);        quitOnError(*err, __LINE__, stderr);
  }
  else if (!strcmp(simulName, "paperIII_ABC_gauss") || !strcmp(simulName, "paperIII_ABC_star") || !strcmp(simulName, "paperIII_ABC_mrlens")) {
    read_cosmo_hm("../param/hmParam_paperIII.par", &cmhm, err);     quitOnError(*err, __LINE__, stderr);
    peak = _paperIII__initialize_peak_param_ABC(simulName, err);    quitOnError(*err, __LINE__, stderr);
  }
  else {
    *err = addError(peak_unknown, "Unknown simulName", *err, __LINE__);
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Reinitialization, if necessary
  if (argc == 7 && (task == 7 || task == 8 || task == 31)) {
    double Omega_m = atof(arg4); 
    double sigma_8 = atof(arg5);
    double w0_de   = atof(arg6);
    char newSimulName[STRING_LENGTH_MAX];
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, w0_de, err); quitOnError(*err, __LINE__, stderr);
    if (task == 8) {
      sprintf(newSimulName, "OmegaM%g_sigma%g_w0de%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8, cmhm->cosmo->w0_de);
      peak->simulName = newSimulName;
    }
    else if (task == 31) {
      sprintf(newSimulName, "paperIII_OmegaM%g_sigma%g_w0de%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8, cmhm->cosmo->w0_de);
      peak->simulName = newSimulName;
    }
  }
  else if (argc == 5 && task == 21) {
    int length = 10000;
    double_arr *Omega_m_arr = initialize_double_arr(length);
    double_arr *sigma_8_arr = initialize_double_arr(length);
    int begin = 0;
    begin = _paperII__generatePoints_zone1(Omega_m_arr, sigma_8_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    begin = _paperII__generatePoints_zone2(Omega_m_arr, sigma_8_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    begin = _paperII__generatePoints_zone3(Omega_m_arr, sigma_8_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    
    int index = atoi(arg4);
    peak->cosmoInd = index;
    double Omega_m = Omega_m_arr->array[index]; 
    double sigma_8 = sigma_8_arr->array[index];
    free_double_arr(Omega_m_arr);
    free_double_arr(sigma_8_arr);
    
    char newSimulName[STRING_LENGTH_MAX];
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, -1, err); quitOnError(*err, __LINE__, stderr);
    sprintf(newSimulName, "paperII_OmegaM%g_sigma%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8);
    peak->simulName = newSimulName;
  }
  else if ((argc == 5 || argc == 6) && task == 31) {
    int length = 40000;
    double_arr *Omega_m_arr = initialize_double_arr(length);
    double_arr *sigma_8_arr = initialize_double_arr(length);
    double_arr *w0_de_arr   = initialize_double_arr(length);
    int begin = 0;
    begin = _paperIII__generatePoints_series1(Omega_m_arr, sigma_8_arr, w0_de_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    begin = _paperIII__generatePoints_series2(Omega_m_arr, sigma_8_arr, w0_de_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    begin = _paperIII__generatePoints_series3(Omega_m_arr, sigma_8_arr, w0_de_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    begin = _paperIII__generatePoints_series4(Omega_m_arr, sigma_8_arr, w0_de_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    begin = _paperIII__generatePoints_series5(Omega_m_arr, sigma_8_arr, w0_de_arr, begin, err); quitOnError(*err, __LINE__, stderr);
    
    int index = atoi(arg4);
    peak->cosmoInd = index;
    double Omega_m = Omega_m_arr->array[index];
    double sigma_8 = sigma_8_arr->array[index];
    double w0_de   = w0_de_arr->array[index];
    free_double_arr(Omega_m_arr);
    free_double_arr(sigma_8_arr);
    free_double_arr(w0_de_arr);
    
    char newSimulName[STRING_LENGTH_MAX];
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, w0_de, err); quitOnError(*err, __LINE__, stderr);
    sprintf(newSimulName, "paperIII_OmegaM%g_sigma%g_w0de%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8, cmhm->cosmo->w0_de);
    peak->simulName = newSimulName;
  }
  
  //-- Precalculate some parameters, always put this after reinitialization
  set_peak_param(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  peak->MPISize = MPISize;
  peak->MPIInd  = MPIInd;
    
  //-- Print, only for MPI root processor
  if (MPIInd == 0) {
    printf("\n");
    printf("                  *********************************\n");
    printf("                  **      Camelus v1.3 beta      **\n");
    printf("                  **  Chieh-An Lin (CEA Saclay)  **\n");
    printf("                  *********************************\n");
    printf("\n");
    printf("Initialization done\n");
    printf("------------------------------------------------------------------------\n");
    
    if (task == 9 || task == 24 || task == 35)                                                  printParam_ABC(cmhm, peak);
    else if (!strcmp(simulName, "default")  || !strcmp(simulName, "test"))                      printParam(cmhm, peak);
    else if (strstr(simulName, "aardvark") != NULL)                                    _paperI__printParam_aardvark(cmhm, peak);
    else if (!strcmp(simulName, "paperII")  || strstr(simulName, "paperII_") != NULL) _paperII__printParam(cmhm, peak);
    else if (!strcmp(simulName, "paperIII"))                                         _paperIII__printParam(cmhm, peak);
    else if (strstr(simulName, "paperIII_ABC_") != NULL)                                        printParam_ABC(cmhm, peak);
    //else if (!strcmp(simulName, "paperIII"))                                                    printParam_complete(cmhm, peak);
  }
  
  //-- Synchronize all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
    
  //-- Below are different functions of Camelus
  
  //-- Sandbox
  if (task == -1) {sandbox(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);}
  
  //-- Print mass function
  else if (task == 1) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    double z = atof(arg3);
    char *fileName = arg4;
    peak->M_min = 1e+09;
    peak->M_max = 1e+17;
    outputMassFct(fileName, cmhm, peak, z, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Output mass functions from z = 0 to 1
  else if (task == 2) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    peak->M_min = 1e+09;
    peak->M_max = 1e+17;
    doMassFct(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fast simulation
  else if (task == 3) {
    if (!(argc == 3 || argc == 5)) {printInstructions(task, 1); return 1;}
    if (argc == 5) {
      peak->M_min = atof(arg3);
      peak->M_max = atof(arg4);
    }
    doFastSimulation(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Ray-tracing
  else if (task == 4) {
    if (!(argc == 4 || argc == 5)) {printInstructions(task, 1); return 1;}
    int doNoise = atoi(arg3);
    char *fileName = arg4;
    if (argc == 4) fileName = NULL;
    doRayTracing(fileName, cmhm, peak, doNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Compute K map
  else if (task == 5) {
    if (!(argc == 4 || argc == 5)) {printInstructions(task, 1); return 1;}
    int doNoise    = atoi(arg3);
    char *fileName = arg4;
    if (argc == 4) fileName = NULL;
    doKMap(fileName, cmhm, peak, doNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list
  else if (task == 6) {
    if (!(argc == 3 || argc == 4)) {printInstructions(task, 1); return 1;}
    char *fileName = arg3;
    if (argc == 3) fileName = NULL;
    doPeakList(fileName, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Repeat peak list
  else if (task == 7) {
    if (!(argc == 4 || argc == 7)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg3);
    doPeakList_repeat(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Data matrix
  else if (task == 8) {
    if (!(argc == 4 || argc == 7)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg3);
    doDataMatrix(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 9) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Mass sheet
  else if (task == 10) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    double z_halo_max = atof(arg3);
    double M_min = atof(arg4);
    doMassSheet(cmhm, peak, z_halo_max, M_min, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Mass function for Aardvark
  else if (task == 11) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    _paperI__doMassFct_aardvark(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Halo list for Aardvark
  else if (task == 12) {
    if (argc != 4) {printInstructions(task, 1); return 1;}
    int fastNb = atoi(arg3);
    _paperI__doHaloList_aardvark(cmhm, peak, fastNb, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Noise list for Aardvark
  else if (task == 13) {
    if (argc != 4) {printInstructions(task, 1); return 1;}
    int noiseNb = atoi(arg3);
    _paperI__doNoiseList_aardvark(cmhm, peak, noiseNb, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- kappa/K map for Aardvark
  else if (task == 14) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int fastNb      = atoi(arg3);
    int doSmoothing = atoi(arg4);
    _paperI__doKMap_aardvark(cmhm, peak, fastNb, doSmoothing, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- n/N map for Aardvark
  else if (task == 15) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int noiseNb     = atoi(arg3);
    int doSmoothing = atoi(arg4);
    _paperI__doNMap_aardvark(cmhm, peak, noiseNb, doSmoothing, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- kn/KN map for Aardvark
  else if (task == 16) {
    if (argc != 6) {printInstructions(task, 1); return 1;}
    int fastNb      = atoi(arg3);
    int noiseNb     = atoi(arg4);
    int doSmoothing = atoi(arg5);
    _paperI__doKNMap_aardvark(cmhm, peak, fastNb, noiseNb, doSmoothing, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list for Aardvark
  else if (task == 17) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int fastNb      = atoi(arg3);
    int noiseNb     = atoi(arg4);
    _paperI__doPeakList_aardvark(cmhm, peak, fastNb, noiseNb, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Pipeline for Paper I
  else if (task == 18) {
    if (!(argc == 5 || argc == 6)) {printInstructions(task, 1); return 1;}
    int fastNb, nbNoise, pipelineType = atoi(arg3);
    if (pipelineType == 2) {
      fastNb  = atoi(arg4);
      nbNoise = atoi(arg5);
    }
    else {
      fastNb  = 0;
      nbNoise = atoi(arg4);
    }
    _paperI__doPipeline_aardvark(cmhm, peak, pipelineType, fastNb, nbNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Full analysis for Paper I
  else if (task == 19) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int nbFast  = atoi(arg3);
    int nbNoise = atoi(arg4);
    _paperI__doFullAnalysis_aardvark(cmhm, peak, nbFast, nbNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- FSL model
  else if (task == 20) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    _FSL__doPeakPrediction(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list for Paper II
  else if (task == 21) {
    if (!(argc == 4 || argc == 5)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg3);
    _paperII__doPeakList_repeat(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fiducial peak list for Paper II
  else if (task == 22) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    _paperII__doPeakList_fidu(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Grid points for Paper II
  else if (task == 23) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    _paperII__doPrior(err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 24) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    _paperII__doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list for Paper III
  else if (task == 31) {
    if (!(argc == 4 || argc == 5 || argc == 6 || argc == 7)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg3);
    char *path;
    if (argc == 6) path = arg5;
    else           path = "../data/kappaMap_paperIII";
    _paperIII__doDataMatrix(cmhm, peak, N, path, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Grid points for Paper III
  else if (task == 32) {
    if (!(argc == 3 || argc == 4)) {printInstructions(task, 1); return 1;}
    int ind = 7941;
    if (argc == 4) ind = atof(arg3);
    _paperIII__doPrior(ind, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Unsmoothed maps for nonlinear filtering for Paper III
  else if (task == 33) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    
    int begin = atoi(arg3);
    int end   = atoi(arg4);
    _paperIII__doAllFilters(cmhm, peak, begin, end, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Unsmoothed maps for nonlinear filtering for Paper III
  else if (task == 34) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    _paperIII__doMapProcessing(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 35) {
    if (argc!=3) {printInstructions(task, 1); return 1;}
    _paperIII__doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  else {
    printInstructions(-2, 1);
    return 1;
  }
    
  free_peak_param(peak);
  free_parameters_hm(&cmhm);
  
  //-- Stop stopwatch
  clock_t finish = clock();
  if (MPIInd == 0) {
    printTime(start, finish);
    printf("\n");
  }
  
  MPI_Finalize();
  return 0;
}

void printInstructions(int task, int printHeader)
{
  if (printHeader) {
    if (task == -1) { //-- Complete menu
      printf("Use Camelus with:\n");
      printInstructions(1, 0); printInstructions(2, 0); printInstructions(3, 0); printInstructions(4, 0); printInstructions(5, 0); 
      printInstructions(6, 0); printInstructions(7, 0); printInstructions(8, 0); printInstructions(9, 0);
      printf("\n");
      printInstructions(31, 1); printInstructions(21, 1); printInstructions(11, 1);
    }
    else if (task == -2) { //-- Normal menu
      printf("Use Camelus with:\n");
      printInstructions(5, 0); printInstructions(8, 0); printInstructions(9, 0);
      printf("\n");
      printInstructions(31, 1);
    }
    else if (task <= 10) {
      printf("Use Camelus with:\n");
      printInstructions(task, 0);
      printf("\n");
    }
    else if (task <= 20) {
      printf("For Paper I:\n");
      printInstructions(11, 0);
      printf("\n");
    }
    else if (task <= 30) {
      printf("For paper II:\n");
      printInstructions(21, 0);
      printf("\n");
    }
    else if (task <= 35) {
      printf("For paper III:\n");
      printInstructions(31, 0);
      printf("\n");
    }
  }
  
  else {
    switch (task) {
      case 1:
	printf("  ./camelus default   1 z fileName                  # Print mass function at z >= 0\n");
	break;
      case 2:
	printf("  ./camelus default   2                             # Output mass functions from z = 0 to 1\n");
	break;
      case 3:
	printf("  ./camelus default   3                             # Fast simulation\n");
	printf("  ./camelus default   3 M_min M_max                 # Fast simulation a given [M_min, M_max]\n");
	break;
      case 4:
	printf("  ./camelus default   4 doNoise                     # Ray-tracing from mass function\n");
	printf("  ./camelus default   4 doNoise fileName            # Ray-tracing with a given halo catalogue\n");
	break;
      case 5:
	printf("  ./camelus default   5 doNoise                     # Lensing map and intermediate products for the first filter\n");
	printf("  ./camelus default   5 doNoise fileName            # Lensing map from a given halo catalogue\n");
	break;
      case 6:
	printf("  ./camelus default   6                             # Peak list and histogram from mass function\n");
	printf("  ./camelus default   6 fileName                    # Peak list and histogram from a kappa map\n");
	break;
      case 7:
	printf("  ./camelus default   7 N                           # N independent peak lists\n");
	printf("  ./camelus default   7 N Omega_m sigma_8 w0_de     # N independent peak lists for a given (Omega_m, sigma_8, w0_de)\n");
	break;
      case 8:
	printf("  ./camelus default   8 N                           # Data matrix of N realizations\n");
	printf("  ./camelus default   8 N Omega_m sigma_8 w0_de     # Data matrix of N realizations for a given (Omega_m, sigma_8, w0_de)\n");
	break;
      case 9:
	printf("  ./camelus default   9                             # ABC analysis\n");
	break;
      case 10:
	printf("  ./camelus default  10 z_halo_max M_min            # Compute mass sheet between definitions of kappa\n");
	break;
      case 11:
	printf("  ./camelus paperI   11                             # Mass function\n");
	printf("  ./camelus paperI   12 fastNb                      # Halo list\n");
	printf("  ./camelus paperI   13 noiseNb                     # Noise list\n");
	printf("  ./camelus paperI   14 fastNb doSmoothing          # kappa or K maps\n");
	printf("  ./camelus paperI   15 noiseNb doSmoothing         # n or N maps\n");
	printf("  ./camelus paperI   16 fastNb noiseNb doSmoothing  # k+n or K+N maps\n");
	printf("  ./camelus paperI   17 fastNb noiseNb              # Peaks\n");
	printf("  ./camelus paperI   18 1 noiseNb                   # Pipeline for noise\n");
	printf("  ./camelus paperI   18 2 fastNb nbNoise            # Pipeline for peaks\n");
	printf("  ./camelus paperI   19 nbFast nbNoise              # Full analysis for Paper I\n");
	printf("  ./camelus paperI   20                             # FSL model\n");
	break;
      case 21:
	printf("  ./camelus paperII  21 N                           # N independent peak lists\n");
	printf("  ./camelus paperII  21 N ind                       # N independent peak lists for a given index\n");
	printf("  ./camelus paperII  22                             # Create black box peak list\n");
	printf("  ./camelus paperII  23                             # Print grid points\n");
	printf("  ./camelus paperII  24                             # ABC analysis\n");
	break;
      case 31:
	printf("  ./camelus paperIII 31 1                           # Output observation\n");
	printf("  ./camelus paperIII 31 N                           # Data matrix of N realizations\n");
	printf("  ./camelus paperIII 31 N Omega_m sigma_8 w0_de     # Data matrix of N realizations for a given (Omega_m, sigma_8, w0_de)\n");
	printf("  ./camelus paperIII 31 N ind                       # Data matrix of N realizations for a given index\n");
	printf("  ./camelus paperIII 32                             # Print grid points\n");
	printf("  ./camelus paperIII 32 ind                         # Print grid points for a given index\n");
	printf("  ./camelus paperIII 33 begin end                   # Output all kappa maps in FITS, index = [begin, end[\n");
	printf("  ./camelus paperIII 34                             # Output map processing with mask\n");
	printf("  ./camelus paperIII_ABC_gauss  35                  # ABC analysis\n");
	printf("  ./camelus paperIII_ABC_star   35                  # ABC analysis\n");
	printf("  ./camelus paperIII_ABC_mrlens 35                  # ABC analysis\n");
	break;
    }
  }
  
  return;
}
#endif

void MPI_finalize(int MPISize, int MPIInd)
{
  MPI_Barrier(MPI_COMM_WORLD);
  if (MPIInd == 0) {
    sleep(1.0);
    printf("------------------------------------------------------------------------\n");
  }
  return;
}

void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  //printParam_complete(cmhm, peak);
  
  double z_l = 0.35;
  double M   = 1.2e+15;
  double z_s = 1.0;
  //doProfile("test", cmhm, peak, z_l, M, z_s, err); forwardError(*err, __LINE__,);
  
    
  
  //mpirun -n 3 ./camelus paperIII_ABC_gauss 35
  //-- Hold root processor if others have not finished
  MPI_finalize(peak->MPISize, peak->MPIInd);
  return;
}

