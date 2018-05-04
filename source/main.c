

  /*******************************
   **  main.c			**
   **  Chieh-An Lin		**
   **  Version 2016.03.22	**
   *******************************/


#include "main.h"


#ifdef __releaseMenu__

// comment
int main(int argc, char *argv[])
{
  //int MPISize, MPIInd;
  
  //-- MPI initialization
  //MPI_Init(&argc, &argv);
  //MPI_Comm_size(MPI_COMM_WORLD, &MPISize);
  //MPI_Comm_rank(MPI_COMM_WORLD, &MPIInd);
  
  //-- Start stopwatch
  clock_t start = clock();
  
  //-- Read inputs
  int task   = (argc >= 2) ? atoi(argv[1]) : 0;
  char *arg2 = (argc >= 3) ? argv[2]       : NULL;
  char *arg3 = (argc >= 4) ? argv[3]       : NULL;
  char *arg4 = (argc >= 5) ? argv[4]       : NULL;
  char *arg5 = (argc >= 6) ? argv[5]       : NULL;
  char *arg6 = (argc >= 7) ? argv[6]       : NULL;
  
  //-- Stopping lock for no MPI task
  //if (MPISize > 1 && !(task == -1 || task == 7)) {
  //  printf("Cannot use multiple processors for task %d\n", task);
  //  MPI_Finalize();
  //  return 1;
  //}
  
  error *myerr = NULL, **err = &myerr;
  cosmo_hm *cmhm;
  peak_param *peak;
  
  //-- Initialize cosmology and halo model
  read_cosmo_hm("../param/hmParam.par", &cmhm, err);       quitOnError(*err, __LINE__, stderr);
  peak = malloc_err(sizeof(peak_param), err);              quitOnError(*err, __LINE__, stderr);
  read_peak_param("../param/peakParam.par", peak, err);    quitOnError(*err, __LINE__, stderr);
  
  //-- Reinitialization, if necessary
  if (argc == 6 && task == 7) {
    double Omega_m = atof(arg3); 
    double sigma_8 = atof(arg4);
    double w0_de   = atof(arg5);
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, w0_de, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Precalculate some parameters, always put this after reinitialization
  set_peak_param(cmhm, peak, err);                         quitOnError(*err, __LINE__, stderr);
  peak->MPISize = 1;
  peak->MPIInd  = 0;
  //peak->MPISize = MPISize;
  //peak->MPIInd  = MPIInd;
  
  //-- Print, only for MPI root processor
  //if (MPIInd == 0) {
    printf("\n");
    printf("                  *********************************\n");
    printf("                  **        Camelus v1.31        **\n");
    printf("                  **  Chieh-An Lin (CEA Saclay)  **\n");
    printf("                  *********************************\n");
    printf("\n");
    printf("Initialization done\n");
    printf("------------------------------------------------------------------------\n");
    printParam(cmhm, peak);
  //}
  
  //-- Synchronize all MPI processes
  //MPI_Barrier(MPI_COMM_WORLD);
    
  //-- Below are different functions of Camelus
  
  //-- Sandbox
  if (task == -1) {sandbox(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);}
  
  //-- Mass function
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
  
  //-- Profile
  else if (task == 3) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    double z_l     = atof(arg2);
    double M       = atof(arg3);
    double z_s     = atof(arg4);
    char name[STRING_LENGTH_MAX];
    sprintf(name, "profile_zL%.3f_M%.1e_zS%.3f", z_l, M, z_s);
    doProfile(name, cmhm, peak, z_l, M, z_s, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Compute lensing maps
  else if (task == 4) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doKMap(NULL, cmhm, peak, 1, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list
  else if (task == 5) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doPeakList(NULL, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Multiscale data
  else if (task == 6) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doMultiscale(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Data matrix
  else if (task == 7) {
    if (!(argc == 3 || argc == 6)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    doDataMatrix(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 8) {
    if (argc != 2) {printInstructions(task, 1); return 1;}
    doABC(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
    //-- only catalogue of galaxies and haloes are produced 
  else if (task == 15) {
    if (argc != 4) {printInstructions(task, 1); return 1;}
    char *input_name = arg2;
    char *input_name2 = arg3;
    doProduce_Catalog(input_name,input_name2, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
    //-- read catalogue of halo and galaxy then compute peak count / list / histogramm
  else if (task == 16) {
    if (argc != 4) {printInstructions(task, 1); return 1;}
    char *input_name = arg2;
	char *opt = arg3;
 	doPeakList_withInputs(input_name,opt, cmhm, peak, err);
	quitOnError(*err, __LINE__, stderr);
  }
  else if (task == 151) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    char *input_catHal = arg3;
    char *input_catGal = arg4;
    doProduce_Catalog_N(N,input_catHal,input_catGal, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  else if (task == 161) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    char *input_name = arg3;
	char *opt = arg4;
 	doPeakList_withInputs_N(N,input_name,opt, cmhm, peak, err);
	quitOnError(*err, __LINE__, stderr);
  }
  else if (task == 171) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    char *input_hal = arg2;
    char *input_gal = arg3;
	char *opt = arg4;
 	doPeakList_withInputs_hod(input_hal,input_gal,opt, cmhm, peak, err); 
	quitOnError(*err, __LINE__, stderr);
  }
  else if (task == 999) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    char *input_name = arg3;
    char *input_name2 = arg4;

	printf("Nb realisation : %i \n",N);
	printf("Input param : %s \n",input_name );
	printf("Output CatHalo : %s \n",input_name2 );

  	read_cosmo_hm(input_name, &cmhm, err);       
	quitOnError(*err, __LINE__, stderr);
    doProduce_Catalog_DM_HOD(N,input_name,input_name2, cmhm, peak, err);
	quitOnError(*err, __LINE__, stderr);
  }
  else if (task == 900){
    if (argc != 6) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    char *input_name = arg3;
    char *input_name2 = arg4;
    char *input_name3 = arg5;

	printf("Nb realisation : %i \n",N);
	printf("Input param : %s \n",input_name );
	printf("Output CatHalo : %s \n",input_name2 );
	printf("Output CatGal : %s \n",input_name3 );

  	read_cosmo_hm(input_name, &cmhm, err);       
	quitOnError(*err, __LINE__, stderr);
	doProduce_Catalog_DM_galaxies(N,input_name,input_name2,input_name3, cmhm, peak, err);
	quitOnError(*err, __LINE__, stderr);
  }

   else if (task == 971){
    if (argc != 6) {printInstructions(task, 1); return 1;}
    int N = atoi(arg2);
    char *input_name = arg3;
    char *input_name2 = arg4;
    char *input_name3 = arg5;

	printf("Nb realisation : %i \n",N);
	printf("Input param : %s \n",input_name );
	printf("Output CatHalo : %s \n",input_name2 );
	printf("Output CatGal_lensed : %s \n",input_name3 );

  	read_cosmo_hm(input_name, &cmhm, err);       
	quitOnError(*err, __LINE__, stderr);
	doProduce_Catalog_DM_galaxies(N,input_name,input_name2,input_name3, cmhm, peak, err);
	quitOnError(*err, __LINE__, stderr);
  }
    
  else {
    printInstructions(-1, 1);
    return 1;
  }
    
  free_peak_param(peak);
  free_parameters_hm(&cmhm);
  
  //-- Stop stopwatch
  clock_t finish = clock();
  //if (MPIInd == 0) {
    printTime(start, finish);
    printf("\n");
  //}
  
  //MPI_Finalize();
  return 0;
}

void printInstructions(int task, int printHeader)
{
  if (printHeader) {
    printf("Use Camelus with:\n");
    printInstructions(1, 0); printInstructions(2, 0); printInstructions(3, 0); printInstructions(4, 0); printInstructions(5, 0);
    printInstructions(6, 0); printInstructions(7, 0); printInstructions(8, 0);
    printInstructions(15, 0); printInstructions(16, 0); printInstructions(900, 0);
    printf("\n");
  }
  
  else {
     switch (task) {
        case 1:
           printf("  ./camelus 1 z                	 # Output mass function at z\n");
           break;
        case 2:
           printf("  ./camelus 2                      # Fast simulation\n");
           break;
        case 3:
           printf("  ./camelus 3 z_l M z_s            # Halo lensing profile with 1-h and 2-h terms\n");
           break;
        case 4:
           printf("  ./camelus 4                      # Lensing map and intermediate products for the first filter\n");
           break;
        case 5:
           printf("  ./camelus 5                      # Peak list and histogram from fast simulation\n");
           break;
        case 6:
           printf("  ./camelus 6                     # A realization of multiscale data\n");
           break;
        case 7:
           printf("  ./camelus 7 N                        # Data matrix of N realizations\n");
           printf("  ./camelus 7 N Omega_m sigma_8 w0_de  # Data matrix of N realizations for a given (Omega_m, sigma_8, w0_de)\n");
           break;
        case 8:
           printf("  ./camelus 8                          # ABC analysis\n");
           break;
        case 15:
           printf("  ./camelus 15 halocat galcat          # Creates halo and galaxy catalogues\n");
           
        case 151:
           printf("  ./camelus 151 N halocat galcat        # Creates N halo and galaxy catalogues\n");
           break;
        case 16:
           printf("  ./camelus 16 galcat   end     # Reads galaxy catalogues and creates peak histogram // end name files \n");
        case 161:
           printf("  ./camelus 161 N galcat  end   # Reads N galaxy catalogues and creates peak histogram // end name files \n");
        case 171:
           printf("  ./camelus 171  halocat galcat_nolensed  end   # Reads halo/galaxy catalogues and compute lensing quantities // end name files \n");
	   break;
        case 900:
           printf("  ./camelus 900 N paramhm halocat galcat   # create catalog haloes with Ngal and galaxy catalog \n");
        case 999:
           printf("  ./camelus 999 N paramhm halocat   # create catalog haloes with Ngal \n");
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
  if (MPISize > 1 && !(task == -1 || task == 9 || task == 35)) {
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
  else if (!strcmp(simulName, "paperIII")) {
    read_cosmo_hm("../param/hmParam_paperIII.par", &cmhm, err);     quitOnError(*err, __LINE__, stderr);
    peak = _paperIII__initialize_peak_param(simulName, err);        quitOnError(*err, __LINE__, stderr);
  }
  else if (!strcmp(simulName, "paperIII2") || strstr(simulName, "paperIII_ABC_") != NULL) {
    read_cosmo_hm("../param/hmParam_paperIII.par", &cmhm, err);     quitOnError(*err, __LINE__, stderr);
    peak = _paperIII__initialize_peak_param_ABC(simulName, err);    quitOnError(*err, __LINE__, stderr);
  }
  else {
    *err = addError(peak_unknown, "Unknown simulName", *err, __LINE__);
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Reinitialization, if necessary
  if (argc == 7 && (task == 8 || task == 31)) {
    double Omega_m = atof(arg4);
    double sigma_8 = atof(arg5);
    double w0_de   = atof(arg6);
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, w0_de, err); quitOnError(*err, __LINE__, stderr);
    if (task == 8)       sprintf(peak->simulName, "OmegaM%g_sigEig%g_w0de%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8, cmhm->cosmo->w0_de);
    else if (task == 31) sprintf(peak->simulName, "paperIII_OmegaM%g_sigEig%g_w0de%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8, cmhm->cosmo->w0_de);
  }
  else if ((argc == 5 || argc == 6) && (task == 31 || task == 36)) {
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
    
    cmhm = updateCmhm(cmhm, Omega_m, sigma_8, w0_de, err); quitOnError(*err, __LINE__, stderr);
    sprintf(peak->simulName, "paperIII_OmegaM%g_sigEig%g_w0de%g", cmhm->cosmo->Omega_m, cmhm->cosmo->sigma_8, cmhm->cosmo->w0_de);
  }
  
  //-- Precalculate some parameters, always put this after reinitialization
  set_peak_param(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  peak->MPISize = MPISize;
  peak->MPIInd  = MPIInd;
    
  //-- Print, only for MPI root processor
  if (MPIInd == 0) {
    printf("\n");
    printf("                  *********************************\n");
    printf("                  **      Camelus v1.31 beta     **\n");
    printf("                  **  Chieh-An Lin (CEA Saclay)  **\n");
    printf("                  *********************************\n");
    printf("\n");
    printf("Initialization done\n");
    printf("------------------------------------------------------------------------\n");
    
    if (!strcmp(simulName, "default") || !strcmp(simulName, "test"))                                  printParam(cmhm, peak);
    else if (!strcmp(simulName, "paperIII") || !strcmp(simulName, "paperIII2") || strstr(simulName, "paperIII_ABC_") != NULL) {
      _paperIII__printParam(cmhm, peak);
    }
  }
  
  //-- Synchronize all MPI processes
  MPI_Barrier(MPI_COMM_WORLD);
    
  //-- Below are different functions of Camelus
  
  //-- Sandbox
  if (task == -1) {sandbox(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);}
  else if (task == 0)  {printInstructions(-2, 1); return 1;}
  
  //-- Mass function
  else if (task == 1) {
    if (!(argc == 3 || argc == 5)) {printInstructions(task, 1); return 1;}
    if (argc == 3) doMassFct(cmhm, peak, err);
    else {
      double z = atof(arg3);
      char *fileName = arg4;
      peak->M_min = 1e+09;
      peak->M_max = 1e+17;
      outputMassFct(fileName, cmhm, peak, z, err);
    }
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Fast simulation
  else if (task == 2) {
    if (!(argc == 3 || argc == 5)) {printInstructions(task, 1); return 1;}
    if (argc == 5) {
      peak->M_min = atof(arg3);
      peak->M_max = atof(arg4);
    }
    doFastSimulation(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Profile
  else if (task == 3) {
    if (argc != 7) {printInstructions(task, 1); return 1;}
    double z_l     = atof(arg3);
    double M       = atof(arg4);
    double z_s     = atof(arg5);
    char *fileName = arg6;
    doProfile(fileName, cmhm, peak, z_l, M, z_s, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Ray-tracing
  else if (task == 4) {
    if (!(argc == 4 || argc == 5)) {printInstructions(task, 1); return 1;}
    int doNoise = atoi(arg3);
    char *fileName = (argc == 5) ? arg4 : NULL;
    doRayTracing(fileName, cmhm, peak, doNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Compute lensing maps
  else if (task == 5) {
    if (!(argc == 4 || argc == 5)) {printInstructions(task, 1); return 1;}
    int doNoise    = atoi(arg3);
    char *fileName = (argc == 5) ? arg4 : NULL;
    doKMap(fileName, cmhm, peak, doNoise, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Peak list
  else if (task == 6) {
    if (!(argc == 3 || argc == 4)) {printInstructions(task, 1); return 1;}
    char *fileName = (argc == 4) ? arg3 : NULL;
    doPeakList(fileName, cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Multiscale data
  else if (task == 7) {
    if (argc != 3) {printInstructions(task, 1); return 1;}
    doMultiscale(cmhm, peak, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Data matrix
  else if (task == 8) {
    if (!(argc == 4 || argc == 7)) {printInstructions(task, 1); return 1;}
    int N = atoi(arg3);
    doDataMatrix(cmhm, peak, N, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- ABC
  else if (task == 9) {
    if (!(argc == 3 || argc == 5 || argc == 6)) {printInstructions(task, 1); return 1;}
    int Q_real  = (argc >= 4) ? atoi(arg3) : 0;
    int t       = (argc >= 5) ? atoi(arg4) : 0;
    int procInd = (argc >= 6) ? atoi(arg5) : 0;
    if      (argc == 3) doABC(cmhm, peak, err); 
    else if (argc == 6) doABC_subset(cmhm, peak, Q_real, t, procInd, err);
    else                doABC_gather(cmhm, peak, Q_real, t, err);
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Mass sheet
  else if (task == 10) {
    if (argc != 5) {printInstructions(task, 1); return 1;}
    double z_halo_max = atof(arg3);
    double M_min = atof(arg4);
    doMassSheet(cmhm, peak, z_halo_max, M_min, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Data matrix for Paper III
  else if (task == 31) {
    if (!(argc == 4 || argc == 5 || argc == 6 || argc == 7)) {printInstructions(task, 1); return 1;}
    int N      = atoi(arg3);
    char *path = (argc == 6) ? arg5 : ".";
    _paperIII__doDataMatrix(cmhm, peak, N, path, err); quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Grid points for Paper III
  else if (task == 32) {
    if (!(argc == 3 || argc == 4)) {printInstructions(task, 1); return 1;}
    int ind = (argc == 4) ? atof(arg3) : 7941;
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
    if (!(argc == 3 || argc == 5 || argc == 6 || argc == 7)) {printInstructions(task, 1); return 1;}
    int Q_real     = (argc >= 4) ? atoi(arg3) : 0;
    int t          = (argc >= 5) ? atoi(arg4) : 0;
    int procInd    = (argc >= 6) ? atoi(arg5) : 0;
    char *tempPath = (argc >= 7) ? arg6 : ".";
    if      (argc == 3) _paperIII__doABC_subset(cmhm, peak, peak->ABC_Q, 0, procInd, tempPath, err);
    else if (argc >= 6) _paperIII__doABC_subset(cmhm, peak, Q_real, t, procInd, tempPath, err);
    else                _paperIII__doABC_gather(cmhm, peak, Q_real, t, err);
    quitOnError(*err, __LINE__, stderr);
  }
  
  //-- Data matrix for Paper III
  else if (task == 36) {
    if (!(argc == 4 || argc == 5 || argc == 6)) {printInstructions(task, 1); return 1;}
    int N      = atoi(arg3);
    char *path = (argc == 6) ? arg5 : ".";
    _paperIII__doDataMatrix2(cmhm, peak, N, path, err); quitOnError(*err, __LINE__, stderr);
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
    if (task == -1) { //-- Complete menu
      printf("Use Camelus with:\n");
      printInstructions(1,  0); printInstructions(2, 0); printInstructions(3, 0); printInstructions(4, 0); printInstructions(5, 0); 
      printInstructions(6,  0); printInstructions(7, 0); printInstructions(8, 0); printInstructions(9, 0); printInstructions(10, 0);
      printf("\n");
      printInstructions(31, 1);
    }
    else if (task == -2) { //-- Normal menu
      printf("Use Camelus with:\n");
      printInstructions(5, 0); printInstructions(7, 0); printInstructions(8, 0); printInstructions(9, 0);
      printf("\n");
      printInstructions(31, 1);
    }
    else if (task <= 11) {
      printf("Use Camelus with:\n");
      printInstructions(task, 0);
      printf("\n");
    }
    else if (task >= 31 && task <= 36) {
      printf("For paper III:\n");
      printInstructions(31, 0);
      printf("\n");
    }
  }
  
  else {
    switch (task) {
      case 1:
	printf("  ./camelus default    1                             # Output mass functions from z = 0 to 1\n");
	printf("  ./camelus default    1 z fileName                  # Output mass function at z\n");
	break;
      case 2:
	printf("  ./camelus default    2                             # Fast simulation\n");
	printf("  ./camelus default    2 M_min M_max                 # Fast simulation a given [M_min, M_max]\n");
	break;
      case 3:
	printf("  ./camelus default    3 z_l M z_s fileName          # Halo lensing profile with 1-h and 2-h terms\n");
	break;
      case 4:
	printf("  ./camelus default    4 doNoise                     # Ray-tracing from mass function\n");
	printf("  ./camelus default    4 doNoise fileName            # Ray-tracing with a given halo catalogue\n");
	break;
      case 5:
	printf("  ./camelus default    5 doNoise                     # Lensing map and intermediate products for the first filter\n");
	printf("  ./camelus default    5 doNoise fileName            # Lensing map from a given halo catalogue\n");
	break;
      case 6:
	printf("  ./camelus default    6                             # Peak list and histogram from mass function\n");
	printf("  ./camelus default    6 fileName                    # Peak list and histogram from a kappa map\n");
	break;
      case 7:
	printf("  ./camelus default    7                             # A realization of multiscale data\n");
	break;
      case 8:
	printf("  ./camelus default    8 N                           # Data matrix of N realizations\n");
	printf("  ./camelus default    8 N Omega_m sigma_8 w0_de     # Data matrix of N realizations for a given (Omega_m, sigma_8, w0_de)\n");
	break;
      case 9:
	printf("  ./camelus default    9                             # ABC analysis\n");
	printf("  ./camelus default    9 Q_real t procInd            # Do a subset for ABC analysis\n");
	printf("  ./camelus default    9 Q_real t                    # Gathering subsets of ABC analysis\n");
	break;
      case 10:
	printf("  ./camelus default   10 z_halo_max M_min            # Compute mass sheet between definitions of kappa\n");
	break;
      case 31:
	printf("  ./camelus paperIII  31 1                          # Output observation\n");
	printf("  ./camelus paperIII  31 N                          # Data matrix of N realizations\n");
	printf("  ./camelus paperIII  31 N ind                      # Data matrix of N realizations for a given index\n");
	printf("  ./camelus paperIII  31 N ind tempPath             # Data matrix of N realizations for a given index in tempPath\n");
	printf("  ./camelus paperIII  31 N Omega_m sigma_8 w0_de    # Data matrix of N realizations for a given (Omega_m, sigma_8, w0_de)\n");
	printf("  ./camelus paperIII  32                            # Print grid points\n");
	printf("  ./camelus paperIII  32 ind                        # Print grid points for a given index\n");
	printf("  ./camelus paperIII  33 begin end                  # Output all kappa maps in FITS, index = [begin, end[\n");
	printf("  ./camelus paperIII  34                            # Output map processing with mask\n");
	printf("  ./camelus paperIII_ABC_* 35                           # ABC analysis\n");
	printf("  ./camelus paperIII_ABC_* 35 Q_real t procInd          # Do a subset for ABC analysis\n");
	printf("  ./camelus paperIII_ABC_* 35 Q_real t procInd tempPath # Do a subset for ABC analysis, put mrlens products in tempPath\n");
	printf("  ./camelus paperIII_ABC_* 35 Q_real t                  # Gathering subsets of ABC analysis\n");
	printf("  ./camelus paperIII2 36 1                          # Output observation\n");
	printf("  ./camelus paperIII2 36 N                          # Data matrix of N realizations\n");
	printf("  ./camelus paperIII2 36 N ind tempPath             # Data matrix of N realizations for a given index in tempPath\n");
	break;
    }
  }
  
  return;
}
#endif

void MPI_terminate(int MPISize, int MPIInd)
{
#ifndef __releaseMenu__
  MPI_Barrier(MPI_COMM_WORLD);
  if (MPIInd == 0) {
    sleep(1.0);
    printf("------------------------------------------------------------------------\n");
  }
#endif
  return;
}

void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  char name[STRING_LENGTH_MAX];
  char name2[STRING_LENGTH_MAX];
  //printParam_complete(cmhm, peak);
  
  /*
  double z_l   = 0.35;
  double M     = 2e+14;
  halo_t *halo = malloc(sizeof(halo_t));
  set_halo_t(cmhm, halo, z_l, M, NULL, err);
  
  double R_vir = halo->r_vir / halo->a;
  double w_l   = halo->w;
  int N_w      = 5000;
  double dw    = 2.0 * R_vir / (double)N_w;
  
  double z_s = 1.0;
  interpolator_t *a_inter = initialize_interpolator_t(N_w+1);
  makeAInterpolator(cmhm, a_inter, z_s, err);
  
  double a_s = 1.0 / (1.0 + z_s);
  double w_s = w(cmhm->cosmo, a_s, 0, err);
  double sum = 0.0;
  double w_curr, a_curr;
  int i;
  
  for (i=0; i<N_w; i++) {
    w_curr = w_l - R_vir + (i + 0.5) * dw;
    a_curr = execute_interpolator_t(a_inter, w_curr);
    sum   += f_K(cmhm->cosmo, w_s - w_curr, err) * f_K(cmhm->cosmo, w_curr, err) / a_curr; forwardError(*err, __LINE__,);
  }
  
  double factor = (FOUR_PI_G_OVER_C2 * CRITICAL_DENSITY * cmhm->cosmo->Omega_m * dw) / f_K(cmhm->cosmo, w_s, err);
  printf("delta kappa = %10f\n", sum * factor);
  */
  
  /*
  cosmo *cm = cmhm->cosmo;
  double logkmin = log(k_min);
  double logkmax = log(k_max);
  double dk = (logkmax - logkmin)/(N_k-1.0);
  
  P_NL_fitting(cm, 1.0, 0.0001, err);
  int i, j;
  double PLk, k;
  
  FILE *file = fopen_err("P_L", "w", err); forwardError(*err, __LINE__,);
  for (j=0; j<N_k; j++) {
    k = exp(logkmin+j*dk);
    PLk = P_L(cm, 1.0, k, err);
    fprintf(file, "%e %e %e\n", exp(logkmin+j*dk), exp(cm->P_NL->table[cm->N_a-1][j]), PLk);
  }
  fclose(file);
  */
  
  double nbHalos[9] = {         839,         1346,          851,           1052,         1681,         1010,            812,         1158,          665};
  double zLArr[9]   = {    0.195620,     0.305570,     0.374820,       0.192260,     0.303390,     0.374685,       0.187840,     0.302895,     0.376040};
  double logMArr[9] = {   14.046000,    14.047000,    14.043000,      14.185000,    14.182000,    14.183000,      14.438500,    14.430000,    14.438000};
  double MArr[9]    = {1.111731e+14, 1.114294e+14, 1.104079e+14,   1.531089e+14, 1.520548e+14, 1.524051e+14,   2.744732e+14, 2.691537e+14, 2.741572e+14};
  double z_s = 1.1689;
  
  int ind = 6;
  for (ind=0; ind<9; ind++) {
    sprintf(name, "profile_z%.4f_hBin%d", z_s, ind);
    doProfile(name, cmhm, peak, zLArr[ind], MArr[ind], z_s, err); forwardError(*err, __LINE__,);
  }
  /*
  sprintf(name, "profile_BMO_z%.4f", z_s);
  doProfile(name, cmhm, peak, zLArr[ind], MArr[ind], z_s, err); forwardError(*err, __LINE__,);
  */
  
  //mpirun -n 3 ./camelus paperIII_ABC_gauss 35
  //-- Hold root processor if others have not finished
  MPI_terminate(peak->MPISize, peak->MPIInd);
  return;
}

