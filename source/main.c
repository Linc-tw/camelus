

  /*******************************
   **  main.c			**
   **  Chieh-An Lin		**
   **  Version 2015.04.05	**
   *******************************/


#include "main.h"


#ifdef __releaseMenu__
int main(int argc, char *argv[])
{
  clock_t start = clock();
  
  //-- Read inputs
  char *arg2, *arg3, *arg4;
  int task;
  if (argc>=2) task = atoi(argv[1]);
  if (argc>=3) arg2 = argv[2];
  if (argc>=4) arg3 = argv[3];
  if (argc>=5) arg4 = argv[4];
  
  //-- Initialize cosmology and halo model
  error *myerr = NULL, **err = &myerr;
  cosmo_hm *cmhm;
  peak_param *peak;
  read_cosmo_hm("../param/hmParam.par", &cmhm, err);    quitOnError(*err, __LINE__, stderr);
  peak = malloc_err(sizeof(peak_param), err);           quitOnError(*err, __LINE__, stderr);
  read_peak_param("../param/peakParam.par", peak, err); quitOnError(*err, __LINE__, stderr);
  
  //-- Reinitialization, if necessary
  if (argc == 5 && task == 4) {
    double Omega_m = atof(arg3); 
    double sigma_8 = atof(arg4);
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, err);     quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Precalculate some parameters, always put this after reinitialization
  set_peak_param(cmhm, peak, err);                      quitOnError(*err, __LINE__, stderr);
  
  //-- Print
  printf("\n                  ---   Camelus v1.2   ---\n");
  printf("\nInitialization done\n");
  printf("------------------------------------------------------------------------\n");
  if (task == 5) printParam_ABC(cmhm, peak);
  else           printParam(cmhm, peak);
  
  if (task == -1) {sandbox(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);}
  
  //-- Print mass function
  else if (task == 1) {
    if (argc != 3) {
      printf("Usage:\n");
      printf("  ./camelus 1 z                  # Print halo mass function at z\n");
      printf("\n");
      return 1;
    }
    
    double z = atof(arg2);
    peak->M_min = 1e+09;
    peak->M_max = 1e+17;
    
    char name[STRING_LENGTH_MAX];
    sprintf(name, "massFct_z%.3f", z);
    outputMassFct(name, cmhm, peak, z, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fast simulation
  else if (task == 2) {
    if (argc != 2) {
      printf("Usage:\n");
      printf("  ./camelus 2                    # Fast simulation\n");
      printf("\n");
      return 1;
    }
    
    doFastSimulation(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list
  else if (task == 3) {
    if (argc != 2) {
      printf("Usage:\n");
      printf("  ./camelus 3                    # Peak list and histogram from fast simulation\n");
      printf("\n");
      return 1;
    }
    
    doPeakList(NULL, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Repeat peak list
  else if (task == 4) {
    if (!(argc==3 || argc==5)) {
      printf("Usage:\n");
      printf("  ./camelus 4 N                  # N independent peak lists\n");
      printf("  ./camelus 4 N Omega_m sigma_8  # N independent peak lists for given (Omega_m, sigma_8)\n");
      printf("\n");
      return 1;
    }
    
    int N = atoi(arg2);
    doPeakList_repeat(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 5) {
    if (argc!=2) {
      printf("Usage:\n");
      printf("  ./camelus 5                    # ABC analysis\n");
      printf("\n");
      return 1;
    }
    
    doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  else {
      printf("Usage:\n");
      printf("  ./camelus 1 z                  # Print halo mass function at z\n");
      printf("  ./camelus 2                    # Fast simulation\n");
      printf("  ./camelus 3                    # Peak list and histogram from fast simulation\n");
      printf("  ./camelus 4 N                  # N independent peak lists\n");
      printf("  ./camelus 4 N Omega_m sigma_8  # N independent peak lists for given (Omega_m, sigma_8)\n");
      printf("  ./camelus 5                    # ABC analysis\n");
      printf("\n");
      return 1;
  }
  
  free_peak_param(peak);
  free_parameters_hm(&cmhm);
  
  clock_t finish = clock();
  printTime(start, finish);
  printf("\n");
  return 0;
}

#else
int main(int argc, char *argv[])
{
  clock_t start = clock();
  
  //-- Read inputs
  char *simulName, *arg3, *arg4, *arg5, *arg6;
  int task;
  if (argc==1) simulName = "default";
  if (argc>=2) simulName = argv[1];
  if (argc>=3) task = atoi(argv[2]);
  if (argc>=4) arg3 = argv[3];
  if (argc>=5) arg4 = argv[4];
  if (argc>=6) arg5 = argv[5];
  if (argc>=7) arg6 = argv[6];
  if (argc<=2) task = 0;
  
  //-- Initialize cosmology and halo model
  error *myerr = NULL, **err = &myerr;
  cosmo_hm *cmhm;
  peak_param *peak;
  if (!strcmp(simulName, "default")) {
    read_cosmo_hm("../param/hmParam.par", &cmhm, err);    quitOnError(*err, __LINE__, stderr);
    peak = malloc_err(sizeof(peak_param), err);           quitOnError(*err, __LINE__, stderr);
    read_peak_param("../param/peakParam.par", peak, err); quitOnError(*err, __LINE__, stderr);
  }
  else if (strstr(simulName, "aardvark") != NULL) {
    _paperI__read_cosmo_hm_aardvark(simulName, &cmhm, err);         quitOnError(*err, __LINE__, stderr);
    peak = _paperI__initialize_peak_param_aardvark(simulName, err); quitOnError(*err, __LINE__, stderr);
  }
  else if (strstr(simulName, "paperII") != NULL) {
    read_cosmo_hm("../param/hmParam_paperII.par", &cmhm, err);      quitOnError(*err, __LINE__, stderr);
    peak = _paperII__initialize_peak_param_paperII(simulName, err); quitOnError(*err, __LINE__, stderr);
  }
  else {
    *err = addError(peak_unknown, "Unknown simulName", *err, __LINE__);
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Reinitialization, if necessary
  if (argc == 6 && (task == 9 || task == 21)) {
    double Omega_m = atof(arg4); 
    double sigma_8 = atof(arg5);
    char newSimulName[STRING_LENGTH_MAX];
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, err); quitOnError(*err, __LINE__, stderr);
    sprintf(newSimulName, "paperII_OmegaM%g_sigma%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8);
    peak->simulName = newSimulName;
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
    peak->paperII_index = index;
    double Omega_m = Omega_m_arr->array[index]; 
    double sigma_8 = sigma_8_arr->array[index];
    free_double_arr(Omega_m_arr);
    free_double_arr(sigma_8_arr);
    
    char newSimulName[STRING_LENGTH_MAX];
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, err); quitOnError(*err, __LINE__, stderr);
    sprintf(newSimulName, "paperII_OmegaM%g_sigma%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8);
    peak->simulName = newSimulName;
  }
  
  //-- Precalculate some parameters, always put this after reinitialization
  set_peak_param(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  
  //-- Print
  printf("\n                  ---   Camelus v1.2   ---\n");
  printf("\nInitialization done\n");
  printf("------------------------------------------------------------------------\n");
  //printParam_complete(cmhm, peak);
  if (task == 24)                                                       printParam_ABC(cmhm, peak);
  else if (!strcmp(simulName, "default") || !strcmp(simulName, "test")) printParam(cmhm, peak);
  else if (strstr(simulName, "aardvark") != NULL)              _paperI__printParam_aardvark(cmhm, peak);
  else if (strstr(simulName, "paperII") != NULL)              _paperII__printParam(cmhm, peak);
  
  //-- Below are different functions of Camelus
  
  //-- Sandbox
  if (task == -1) {sandbox(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);}
  
  //-- Print mass function
  else if (task == 1) {
    if (!(argc==4 || argc==6)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 1 z                            # Print mass function at z >= 0\n");
      printf("  ./camelus simulName 1 z M_min M_max                # Print mass function at z >= 0 given [M_min, M_max]\n");
      printf("\n");
      return 1;
    }
    
    double z = atof(arg3);
    if (argc == 6) {
      peak->M_min = atof(arg4);
      peak->M_max = atof(arg5);
    }
    else {
      peak->M_min = 1e+09;
      peak->M_max = 1e+17;
    }
    outputMassFct(NULL, cmhm, peak, z, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Output mass functions from z = 0 to 1
  else if (task == 2) {
    if (argc!=3) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 2                              # Output mass functions from z = 0 to 1\n");
      printf("\n");
      return 1;
    }
    
    doMassFct(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fast simulation
  else if (task == 3) {
    if (!(argc==3 || argc==5)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 3                              # Fast simulation\n");
      printf("  ./camelus simulName 3 M_min M_max                  # Fast simulation given [M_min, M_max]\n");
      printf("\n");
      return 1;
    }
    
    if (argc == 5) {
      peak->M_min = atof(arg3);
      peak->M_max = atof(arg4);
    }
    doFastSimulation(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Ray-tracing
  else if (task == 4) {
    if (!(argc==3 || argc==4)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 4                              # Ray-tracing from mass function\n");
      printf("  ./camelus simulName 4 haloFileName                 # Ray-tracing with a given halo list\n");
      printf("\n");
      return 1;
    }
    
    char *haloFileName = arg3;
    if (argc == 3) haloFileName = NULL;
    doRayTracing(haloFileName, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Compute K map
  else if (task == 5) {
    if (!(argc==3 || argc==4)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 5                              # Smoothed noiseless map from mass function\n");
      printf("  ./camelus simulName 5 haloFileName                 # Smoothed noiseless map from a given halo list\n");
      printf("\n");
      return 1;
    }
    
    char *haloFileName = arg3;
    if (argc == 3) haloFileName = NULL;
    doKMap(haloFileName, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Compute N map
  else if (task == 6) {
    if (argc!=3) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 6                              # Smoothed noise map\n");
      printf("\n");
      return 1;
    }
    
    doNMap(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Compute K+N map
  else if (task == 7) {
    if (!(argc==3 || argc==4 || argc==5)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 7                              # Smoothed noisy map from mass function\n");
      printf("  ./camelus simulName 7 KName                        # Smoothed noisy map from a noiseless map\n");
      printf("  ./camelus simulName 7 KName NName                  # Smoothed noisy map from a noiseless and a noise map\n");
      printf("\n");
      return 1;
    }
    
    char *KName = arg3;
    char *NName = arg4;
    if (argc <= 4) NName = NULL;
    if (argc == 3) KName = NULL;
    doKNMap(KName, NName, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list
  else if (task == 8) {
    if (!(argc==3 || argc==4)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 8                              # Peak list and histogram from mass function\n");
      printf("  ./camelus simulName 8 KNName                       # Peak list and histogram from a noisy map\n");
      printf("\n");
      return 1;
    }
    
    char *KNName = arg3;
    if (argc == 3) KNName = NULL;
    doPeakList(KNName, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Repeat peak list
  else if (task == 9) {
    if (!(argc==4 || argc==6)) {
      printf("Usage: (simulName can be default)\n");
      printf("  ./camelus simulName 9 N                            # N independent peak lists\n");
      printf("  ./camelus simulName 9 N Omega_m sigma_8            # N independent peak lists for given (Omega_m, sigma_8)\n");
      printf("\n");
      return 1;
    }
    
    int N = atoi(arg3);
    doPeakList_repeat(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Mass function for Aardvark
  else if (task == 11) {
    if (argc!=3) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 11                             # Mass function\n");
      printf("\n");
      return 1;
    }
    
    _paperI__doMassFct_aardvark(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Halo list for Aardvark
  else if (task == 12) {
    if (argc!=4) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 12 fastNb                      # Halo list\n");
      printf("\n");
      return 1;
    }
    
    int fastNb = atoi(arg3);
    _paperI__doHaloList_aardvark(cmhm, peak, fastNb, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Noise list for Aardvark
  else if (task == 13) {
    if (argc!=4) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 13 noiseNb                     # Noise list\n");
      printf("\n");
      return 1;
    }
    
    int noiseNb = atoi(arg3);
    _paperI__doNoiseList_aardvark(cmhm, peak, noiseNb, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- kappa/K map for Aardvark
  else if (task == 14) {
    if (argc!=5) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 14 fastNb doSmoothing          # kappa or K maps\n");
      printf("\n");
      return 1;
    }
    
    int fastNb      = atoi(arg3);
    int doSmoothing = atoi(arg4);
    _paperI__doKMap_aardvark(cmhm, peak, fastNb, doSmoothing, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- n/N map for Aardvark
  else if (task == 15) {
    if (argc!=5) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 15 noiseNb doSmoothing         # n or N maps\n");
      printf("\n");
      return 1;
    }
    
    int noiseNb     = atoi(arg3);
    int doSmoothing = atoi(arg4);
    _paperI__doNMap_aardvark(cmhm, peak, noiseNb, doSmoothing, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- kn/KN map for Aardvark
  else if (task == 16) {
    if (argc!=6) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 16 fastNb noiseNb doSmoothing  # k+n or K+N maps\n");
      printf("\n");
      return 1;
    }
    
    int fastNb      = atoi(arg3);
    int noiseNb     = atoi(arg4);
    int doSmoothing = atoi(arg5);
    _paperI__doKNMap_aardvark(cmhm, peak, fastNb, noiseNb, doSmoothing, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list for Aardvark
  else if (task == 17) {
    if (argc!=5) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 17 fastNb noiseNb              # Peaks\n");
      printf("\n");
      return 1;
    }
    
    int fastNb      = atoi(arg3);
    int noiseNb     = atoi(arg4);
    _paperI__doPeakList_aardvark(cmhm, peak, fastNb, noiseNb, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Pipeline for Paper I
  else if (task == 18) {
    if (!(argc==5 || argc==6)) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 18 1 noiseNb                   # Pipeline for noise\n");
      printf("  ./camelus simulName 18 2 fastNb nbNoise            # Pipeline for peaks\n");
      printf("\n");
      return 1;
    }
    
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
    if (argc!=5) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 19 nbFast nbNoise              # Full analysis for Paper I\n");
      printf("\n");
      return 1;
    }
    
    int nbFast  = atoi(arg3);
    int nbNoise = atoi(arg4);
    _paperI__doFullAnalysis_aardvark(cmhm, peak, nbFast, nbNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- FSL model
  else if (task == 20) {
    if (argc!=3) {
      printf("For Paper I:\n");
      printf("  ./camelus simulName 20                             # FSL model\n");
      printf("\n");
      return 1;
    }
    
    _FSL__doPeakPrediction(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list for Paper II
  else if (task == 21) {
    if (!(argc==4 || argc==5 || argc==6)) {
      printf("For paper II:\n");
      printf("  ./camelus simulName 21 N                           # N independent peak lists\n");
      printf("  ./camelus simulName 21 N Omega_m sigma_8           # N independent peak lists for given (Omega_m, sigma_8)\n");
      printf("  ./camelus simulName 21 N index                     # N independent peak lists for a given index\n");
      printf("\n");
      return 1;
    }
    
    int N = atoi(arg3);
    _paperII__doPeakList_repeat(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fiducial peak list for Paper II
  else if (task == 22) {
    if (argc!=3) {
      printf("For paper II:\n");
      printf("  ./camelus simulName 22                             # Create black box peak list\n");
      printf("\n");
      return 1;
    }
    
    _paperII__doPeakList_fidu(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Prior for Paper II
  else if (task == 23) {
    if (argc!=3) {
      printf("For paper II:\n");
      printf("  ./camelus simulName 23                             # Print prior information\n");
      printf("\n");
      return 1;
    }
    
    _paperII__printPrior(err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 24) {
    if (argc!=3) {
      printf("For paper II:\n");
      printf("  ./camelus simulName 24                             # ABC analysis\n");
      printf("\n");
      return 1;
    }
    
    _paperII__doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  else {
      printf("Usage: (simulName can be default)\n");
#ifdef __completeMenu__
      printf("  ./camelus simulName 1 z                            # Print mass function at z >= 0\n");
      printf("  ./camelus simulName 1 z M_min M_max                # Print mass function at z >= 0 given [M_min, M_max]\n");
      printf("  ./camelus simulName 2                              # Output mass functions from z = 0 to 1\n");
      printf("  ./camelus simulName 3                              # Fast simulation\n");
      printf("  ./camelus simulName 3 M_min M_max                  # Fast simulation given [M_min, M_max]\n");
      printf("  ./camelus simulName 4                              # Ray-tracing from mass function\n");
      printf("  ./camelus simulName 4 haloFileName                 # Ray-tracing with a given halo list\n");
      printf("  ./camelus simulName 5                              # Smoothed noiseless map from mass function\n");
      printf("  ./camelus simulName 5 haloFileName                 # Smoothed noiseless map from a given halo list\n");
      printf("  ./camelus simulName 6                              # Smoothed noise map\n");
      printf("  ./camelus simulName 7                              # Smoothed noisy map from mass function\n");
      printf("  ./camelus simulName 7 KName                        # Smoothed noisy map from a noiseless map\n");
      printf("  ./camelus simulName 7 KName NName                  # Smoothed noisy map from a noiseless and a noise map\n");
#endif
      printf("  ./camelus simulName 8                              # Peak list and histogram from mass function\n");
#ifdef __completeMenu__
      printf("  ./camelus simulName 8 KNName                       # Peak list and histogram from a noisy map\n");
#endif
      printf("  ./camelus simulName 9 N                            # N independent peak lists\n");
      printf("  ./camelus simulName 9 N Omega_m sigma_8            # N independent peak lists for given (Omega_m, sigma_8)\n");
      printf("\n");
#ifdef __completeMenu__
      printf("For Paper I:\n");
      printf("  ./camelus simulName 11                             # Mass function\n");
      printf("  ./camelus simulName 12 fastNb                      # Halo list\n");
      printf("  ./camelus simulName 13 noiseNb                     # Noise list\n");
      printf("  ./camelus simulName 14 fastNb doSmoothing          # kappa or K maps\n");
      printf("  ./camelus simulName 15 noiseNb doSmoothing         # n or N maps\n");
      printf("  ./camelus simulName 16 fastNb noiseNb doSmoothing  # k+n or K+N maps\n");
      printf("  ./camelus simulName 17 fastNb noiseNb              # Peaks\n");
      printf("  ./camelus simulName 18 1 noiseNb                   # Pipeline for noise\n");
      printf("  ./camelus simulName 18 2 fastNb nbNoise            # Pipeline for peaks\n");
      printf("  ./camelus simulName 19 nbFast nbNoise              # Full analysis for Paper I\n");
      printf("  ./camelus simulName 20                             # FSL model\n");
      printf("\n");
#endif
      printf("For paper II:\n");
      printf("  ./camelus simulName 21 N                           # N independent peak lists\n");
      printf("  ./camelus simulName 21 N Omega_m sigma_8           # N independent peak lists for given (Omega_m, sigma_8)\n");
      printf("  ./camelus simulName 21 N index                     # N independent peak lists for a given index\n");
      printf("  ./camelus simulName 22                             # Create black box peak list\n");
      printf("  ./camelus simulName 23                             # Print prior information\n");
      printf("  ./camelus simulName 24                             # ABC analysis\n");
      printf("\n");
      return 1;
  }
  
  free_peak_param(peak);
  free_parameters_hm(&cmhm);
  
  clock_t finish = clock();
  printTime(start, finish);
  printf("\n");
  return 0;
}
#endif


void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  
  printf("This is a sandbox. Test me!");
  
  printf("------------------------------------------------------------------------\n");
  return;
}

