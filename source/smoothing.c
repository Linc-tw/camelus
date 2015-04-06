

  /*******************************
   **  smoothing.c		**
   **  Chieh-An Lin		**
   **  Version 2015.02.25	**
   *******************************/


#include "smoothing.h"


//----------------------------------------------------------------------
//-- Functions related to map_t

map_t *initialize_map_t(int N1, int N2, double theta_pix, error **err)
{
  map_t *kMap      = (map_t*)malloc_err(sizeof(map_t), err);        forwardError(*err, __LINE__,);
  kMap->N1         = N1;
  kMap->N2         = N2;
  kMap->length     = N1 * N2;
  kMap->theta_pix  = theta_pix;
  kMap->center[0]  = 0.0;
  kMap->center[1]  = 0.0;
  kMap->limits[0]  = -0.5 * N1 * theta_pix;
  kMap->limits[1]  = -kMap->limits[0];
  kMap->limits[2]  = -0.5 * N2 * theta_pix;
  kMap->limits[3]  = -kMap->limits[2];
  kMap->type       = kappa_map;
  kMap->kappa_mean = 0.0;
  kMap->kappa      = (double*)malloc_err(kMap->length * sizeof(double), err); forwardError(*err, __LINE__,);
  return kMap;
}

void free_map_t(map_t *kMap)
{
  if (kMap->kappa) {free(kMap->kappa); kMap->kappa = NULL;}
  free(kMap); kMap = NULL;
  return;
}

void getPixPos_map_t(map_t *kMap, double pos[2], int i, int j)
{
  //-- Compute the position of the pixel indexed (i, j) and stock this information in pos.
  pos[0] = kMap->limits[0] + (0.5 + i) * kMap->theta_pix;
  pos[1] = kMap->limits[2] + (0.5 + j) * kMap->theta_pix;
  return;
}

void subtractMean_map_t(peak_param *peak, map_t *kMap)
{
  //-- Get mean
  double mean = 0.0;
  int i;
  for (i=0; i<kMap->length; i++) mean += kMap->kappa[i];
  mean /= (double)kMap->length;
  kMap->kappa_mean = mean;
  
  //-- Subtract
  for (i=0; i<kMap->length; i++) kMap->kappa[i] -= mean;
  
  if (peak->printMode == 1) {
    printf("kappa mean = %.5f, ", mean);
    fflush(stdout);
  }
  else printf("kappa mean %g subtracted\n", mean);
  return;
}

void read_map_t(char *name, peak_param *peak, map_t *kMap, mapType_t type, error **err)
{
  testError(kMap->length!=peak->resol[0]*peak->resol[1], peak_match, "Resolution match error", *err, __LINE__); forwardError(*err, __LINE__,);
  FILE *file = fopen_err(name, "r", err); forwardError(*err, __LINE__,);
  
  char buffer[STRING_LENGTH_MAX], *buffer1;
  int buffer2, count = 0;
  double *kappa = kMap->kappa;
  
  //-- Read
  int c = fgetc(file);
  while (c != EOF) {
    if (c == (int)'#') buffer1 = fgets(buffer, STRING_LENGTH_MAX, file);
    else {
      testError(count>=kMap->length, peak_overflow, "Too many pixels", *err, __LINE__); forwardError(*err, __LINE__,);
      ungetc(c, file);
      buffer2 = fscanf(file, "%*f %*f %lf\n", &kappa[count]);
      count++;
    }
    c = fgetc(file);
  }
  fclose(file);
  
  //-- Length check
  testError(count!=kMap->length, peak_match, "Array length error", *err, __LINE__); forwardError(*err, __LINE__,);
  kMap->type = type;
  printf("\"%s\" read\n", name);
  return;
}

void output_map_t(FILE *file, map_t *kMap)
{
  mapType_t type = kMap->type;
  
  if (kMap->kappa_mean == 0) fprintf(file, "# Number of pixels = %d\n", kMap->length);
  else                       fprintf(file, "# Number of pixels = %d, kappa mean = %g\n", kMap->length, kMap->kappa_mean);
  fprintf(file, "#\n");
  if      (type == kappa_map) fprintf(file, "#  theta_x   theta_y     kappa\n");
  else if (type == K_map)     fprintf(file, "#  theta_x   theta_y        K \n");
  else if (type == N_map)     fprintf(file, "#  theta_x   theta_y        N \n");
  else if (type == kn_map)    fprintf(file, "#  theta_x   theta_y       k+n\n");
  else if (type == KN_map)    fprintf(file, "#  theta_x   theta_y       K+N\n");
  else                        fprintf(file, "#  theta_x   theta_y     noise\n"); //-- if (type == noise_map || type == noise_list)
  fprintf(file, "# [arcmin]  [arcmin]       [-]\n");
  
  double pos[2];
  int i, j;
  for (j=0; j<kMap->N2; j++) {
    for (i=0; i<kMap->N1; i++) {
      getPixPos_map_t(kMap, pos, i, j);
      fprintf(file, "  %8.3f  %8.3f  %8.5f\n", pos[0], pos[1], kMap->kappa[i+j*kMap->N1]);
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to smoothing and noise

void fillGaussianKernel(FFT_t *transformer, double s)
{
  //-- A kernel should be placed in this way
  //--  1.0  0.6   0   0   0  0.6
  //--  0.6  0.1   0   0   0  0.1
  //--    0    0   0   0   0    0
  //--    0    0   0   0   0    0
  //--    0    0   0   0   0    0
  //--  0.6  0.1   0   0   0  0.1
  //-- 
  //-- Gaussian kernel:
  //-- exp[-(x^2 + y^2)/s^2]
  //-- where s is theta_G / theta_pix
  
  double *kernel = transformer->kernel;
  double s_sq    = SQ(s);                    //-- [pix^2]
  double cut     = CUTOFF_FACTOR_FILTER * s; //-- CUTOFF_FACTOR_FILTER set to 2.2 in peakParameters.h, equivalent to a Gaussian trunked at 3 sigma
  double cut_sq  = SQ(cut);                  //-- [pix^2]
  
  int M1    = transformer->N1; //-- M1, M2 are resolution + size of zero padding
  int M2    = transformer->N2;
  int upper = (int)ceil(cut);
  
  double value, sum = 0.0;
  int i, j, r_sq;
  for (i=-upper; i<=upper; i++) {
    for (j=-upper; j<=upper; j++) {
      r_sq = SUM_SQ_2(i, j);
      if (r_sq > cut_sq) continue;
      value = exp(-r_sq / s_sq);
      sum += value;
      kernel[imod(M1, i) + imod(M2, j) * M1] = value;
    }
  }
  
  //-- Normalization
  for (i=0; i<transformer->length; i++) kernel[i] /= sum;
  
  //-- To Fourier space
  fftw_execute(transformer->kernel_p);
  return;
}

void doFFTSmoothing(peak_param *peak, map_t *kMap, FFT_t *transformer, error **err)
{
  testError((kMap->type!=kappa_map && kMap->type!=noise_map && kMap->type!=kn_map), peak_mapType, "Map type error", *err, __LINE__); forwardError(*err, __LINE__,);
  reset_FFT_t(transformer);
  
  int M1 = transformer->N1;
  int N1 = kMap->N1;
  int N2 = kMap->N2;
  
  //-- Fill transformer->before and add zero-padding
  int i, j;
  for (i=0; i<N1; i++) {
    for (j=0; j<N2; j++) {
      transformer->before[i+j*M1] = kMap->kappa[i+j*N1];
    }
  }
  
  execute_FFT_t(transformer);
  
  //-- Retrieve values from transformer->after, set to smoothed->kappa
  for (i=0; i<N1; i++) {
    for (j=0; j<N2; j++) {
      kMap->kappa[i+j*N1] = transformer->after[i+j*M1];
    }
  }
  
  if (kMap->type == kappa_map)      kMap->type = K_map;
  else if (kMap->type == noise_map) kMap->type = N_map;
  else                              kMap->type = KN_map;
  
  if (peak->printMode != 1) printf("Smoothing with FFT done\n");
  return;
}

void fillNoise(peak_param *peak, map_t *nMap, u_int64_t seed)
{
  gsl_rng *generator = peak->generator;
  double sigma_pix = peak->sigma_pix;
  gsl_rng_set(generator, seed);
  
  int i;
  for (i=0; i<nMap->length; i++) {
    nMap->kappa[i] = gsl_ran_gaussian(generator, sigma_pix);
  }
  nMap->type = noise_map;
  
  if (peak->printMode != 1) printf("%d samples of Gaussian noise generated\n", nMap->length);
  return;
}

void addNoise(map_t *kMap, map_t *nMap, error **err)
{
  testError(kMap->N1!=nMap->N1, peak_match, "Map size error", *err, __LINE__); forwardError(*err, __LINE__,);
  testError(kMap->N2!=nMap->N2, peak_match, "Map size error", *err, __LINE__); forwardError(*err, __LINE__,);
  testError(!((kMap->type==kappa_map && nMap->type==noise_map) || (kMap->type==K_map && nMap->type==N_map)), peak_mapType, "Map type error", *err, __LINE__); forwardError(*err, __LINE__,);
  
  int i;
  for (i=0; i<kMap->length; i++) kMap->kappa[i] += nMap->kappa[i];
  if (kMap->type==kappa_map) kMap->type = kn_map;
  else                       kMap->type = KN_map;
  return;
}

//----------------------------------------------------------------------
//-- Functions related to map making

void kappaMapFromKappa(gal_map *gMap, map_t *kMap, error **err)
{
  testError(gMap->length!=kMap->length, peak_match, "Resolution match error", *err, __LINE__); forwardError(*err, __LINE__,);
  int i;
  for (i=0; i<gMap->length; i++) kMap->kappa[i] = kappaMean_gal_list(gMap->map[i], err);       forwardError(*err, __LINE__,);
  kMap->type = kappa_map;
  return;
}

void outputMap(char name[], cosmo_hm *cmhm, peak_param *peak, map_t *kMap)
{
  mapType_t type = kMap->type;
  FILE *file     = fopen(name, "w");
  
  fprintf(file, "# %s\n", STR_MAPTYPE_T(type));
  fprintf(file, "# Field = %s, Omega = (%g, %g) [arcmin], theta_pix = %g [arcmin]\n", STR_FIELD_T(peak->field), peak->Omega[0], peak->Omega[1], peak->theta_pix);
  
  if (type != noise_map && type != N_map)     fprintf(file, "# n_gal = %g [arcmin^-2], z_s = %g\n", peak->n_gal, peak->z_s);
  if (type != kappa_map && type != noise_map) fprintf(file, "# Filter = %s, theta_G = %g [arcmin], s = %g [pix]\n", STR_FILTER_T(peak->filter), peak->theta_G, peak->s);
  if (type != kappa_map && type != K_map)     fprintf(file, "# sigma_eps = %g, sigma_pix = %g\n", peak->sigma_eps, peak->sigma_pix);
  fprintf(file, "#\n");
  
  if (type != noise_map && type != N_map) {
    outputCosmoParam(file, cmhm, peak);
    fprintf(file, "#\n");
  }
  
  output_map_t(file, kMap);
  
  fclose(file);
  printf("\"%s\" made\n", name);
  return;
}

void makeMap(peak_param *peak, gal_map *gMap, map_t *kMap, map_t *nMap, FFT_t *transformer, error **err)
{
  //-- Map making main function
  
  if (peak->doKappa == 1) {
    //-- Do convergence
    
    kappaMapFromKappa(gMap, kMap, err);           forwardError(*err, __LINE__,);
    subtractMean_map_t(peak, kMap);
    if (peak->doNewNoise == 1) {
      //-- Overwrite nMap
      u_int64_t seed = renewSeed();
      fillNoise(peak, nMap, seed);
    }
    addNoise(kMap, nMap, err);                    forwardError(*err, __LINE__,);
    doFFTSmoothing(peak, kMap, transformer, err); forwardError(*err, __LINE__,);
  }
  else {
    //-- TODO: Do shear
    
    if (peak->z_s > 0) {
      //-- Fixed redshift
    }
    else {
      //-- Real redshift
    }
  }
  return;
}

//----------------------------------------------------------------------
//-- Main functions

void doKMap(char haloFileName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  //-- K map for "smoothed noiseless map"
  halo_map *hMap     = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  gal_map *gMap      = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
  makeGalaxies(cmhm, peak, gMap, err);                                                            forwardError(*err, __LINE__,);
  map_t *kMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err);    forwardError(*err, __LINE__,);
  FFT_t *transformer = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
  fillGaussianKernel(transformer, peak->s);
  
  if (haloFileName == NULL) {
    //-- Carry out fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
    setMassSamplers(cmhm, peak, sampArr, err);     forwardError(*err, __LINE__,);
    makeFastSimul(cmhm, peak, sampArr, hMap, err); forwardError(*err, __LINE__,);
    outputFastSimul("haloList", cmhm, peak, hMap);
    free_sampler_arr(sampArr);
  }
  else {
    read_halo_map(haloFileName, cmhm, hMap, err);  forwardError(*err, __LINE__,);
  }
  
  lensingForMap(cmhm, peak, hMap, gMap, err);      forwardError(*err, __LINE__,);
  outputGalaxies("galList", cmhm, peak, gMap);
  
  kappaMapFromKappa(gMap, kMap, err);              forwardError(*err, __LINE__,);
  subtractMean_map_t(peak, kMap);
  doFFTSmoothing(peak, kMap, transformer, err);    forwardError(*err, __LINE__,);
  outputMap("KMap", cmhm, peak, kMap);
  
  free_halo_map(hMap);
  free_gal_map(gMap);
  free_map_t(kMap);
  free_FFT_t(transformer);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doNMap(cosmo_hm *cmhm, peak_param *peak, error **err)
{
  //-- N map for "smoothed noise map"
  map_t *nMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  FFT_t *transformer = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
  fillGaussianKernel(transformer, peak->s);
  
  u_int64_t seed = renewSeed();
  fillNoise(peak, nMap, seed);
  doFFTSmoothing(peak, nMap, transformer, err); forwardError(*err, __LINE__,);
  outputMap("NMap", cmhm, peak, nMap);
  
  free_map_t(nMap);
  printf("------------------------------------------------------------------------\n");
  return;
}

void doKNMap(char KName[], char NName[], cosmo_hm *cmhm, peak_param *peak, error **err)
{
  map_t *kMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  map_t *nMap        = initialize_map_t(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
  FFT_t *transformer = initialize_FFT_t(peak->FFTSize, peak->FFTSize);
  fillGaussianKernel(transformer, peak->s);
    
  if (KName == NULL) {
    //-- Carry out KMap from fast simulation
    sampler_arr *sampArr = initialize_sampler_arr(peak->N_z_halo, peak->nbMassBins);
    setMassSamplers(cmhm, peak, sampArr, err);                                                        forwardError(*err, __LINE__,);
    halo_map *hMap       = initialize_halo_map(peak->resol[0], peak->resol[1], peak->theta_pix, err); forwardError(*err, __LINE__,);
    gal_map *gMap        = initialize_gal_map(peak->resol[0], peak->resol[1], peak->theta_pix, err);  forwardError(*err, __LINE__,);
    makeGalaxies(cmhm, peak, gMap, err);                                                              forwardError(*err, __LINE__,);
    
    makeFastSimul(cmhm, peak, sampArr, hMap, err); forwardError(*err, __LINE__,);
    lensingForMap(cmhm, peak, hMap, gMap, err);    forwardError(*err, __LINE__,);
    kappaMapFromKappa(gMap, kMap, err);            forwardError(*err, __LINE__,);
    subtractMean_map_t(peak, kMap);
    doFFTSmoothing(peak, kMap, transformer, err);  forwardError(*err, __LINE__,);
    outputMap("KMap", cmhm, peak, kMap);
    
    free_sampler_arr(sampArr);
    free_halo_map(hMap);
    free_gal_map(gMap);
  }
  else {read_map_t(KName, peak, kMap, K_map, err); forwardError(*err, __LINE__,);}
  
  if (NName == NULL) {
    //-- Carry out KMap from fast simulation
    u_int64_t seed = renewSeed();
    fillNoise(peak, nMap, seed);
    doFFTSmoothing(peak, nMap, transformer, err);  forwardError(*err, __LINE__,);
    outputMap("NMap", cmhm, peak, nMap);
  }
  else {read_map_t(NName, peak, nMap, N_map, err); forwardError(*err, __LINE__,);}
  
  addNoise(kMap, nMap, err);                       forwardError(*err, __LINE__,);
  outputMap("KNMap", cmhm, peak, kMap);
  
  free_map_t(kMap);
  free_map_t(nMap);
  free_FFT_t(transformer);
  printf("------------------------------------------------------------------------\n");
  return;
}

//----------------------------------------------------------------------



