

  /*******************************************************
   **  wrapper.cpp					**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#ifdef __cplusplus
#include <stdio.h>
#include <healpix_cxx/alm.h>
#include <healpix_cxx/alm_healpix_tools.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/xcomplex.h>

extern "C" {

#include "FITSFunctions.h"
#include "wrapper.h"
#endif

  
//----------------------------------------------------------------------
//-- Functions related to a_lm

void readWeight(long nside, double *weight)
{
  char name[1024];
  sprintf(name, "/usr/lib/python3/dist-packages/healpy/data/weight_ring_n%05ld.fits", nside);
  FITS_t *fits = initializeTableReader_FITS_t(name);
  readTableColumn(fits, 0, weight);
  int i;
  for (i=0; i<2*nside; i++) weight[i] += 1.0;
  return;
}

void mapToAlm(long nside, float *map, int l_maxx, double *alm, double *weight)
{
  arr<float> fltArr_cxx(map, 12*nside*nside);
  Healpix_Map<float> map_cxx(fltArr_cxx, RING);
  Alm< xcomplex<float> > alm_cxx(l_maxx, l_maxx);
  
  double *weight2;
  if (weight == NULL) {
    weight2 = (double*)malloc(2*nside * sizeof(double));
    readWeight(nside, weight2);
  }
  else weight2 = weight;
  arr<double> weight_cxx(weight2, 2*nside);
  
  int l, m;
  
  //-- Set to 0
  for (l=0; l<=l_maxx; l++) {
    for (m=0; m<=l; m++) alm_cxx(l, m).Set(0.0, 0.0);
  }
  
  map2alm(map_cxx, alm_cxx, weight_cxx);
  
  int index;
  
  //-- Retrieve
  for (m=0; m<=l; m++) {
    for (l=m; l<=l_maxx; l++) {
      index = m*(2*l_maxx+1-m)/2+l;
      alm[0+2*index] = alm_cxx(l, m).real();
      alm[1+2*index] = alm_cxx(l, m).imag();
    }
  }
  
  if (weight == NULL) free(weight2);
  return;
}

void almToMap(int l_maxx, double *alm, long nside, float *map)
{
  Alm< xcomplex<float> > alm_cxx(l_maxx, l_maxx);
  Healpix_Map<float> map_cxx((int)nside, RING, SET_NSIDE);
  
  long l, m, index;
  
  //-- Fill
  for (m=0; m<=l; m++) {
    for (l=m; l<=l_maxx; l++) {
      index = m*(2*l_maxx+1-m)/2+l;
      alm_cxx(l, m).Set(alm[0+2*index], alm[1+2*index]);
    }
  }
  
  alm2map(alm_cxx, map_cxx);
  
  int64 i;
  //-- Retrieve
  for (i=0; i<12*nside*nside; i++) map[i] = (float)map_cxx[i];
  return;
}

void almToMap_spin2(int l_maxx, double *alm, long nside, float *map1, float *map2)
{
  Alm< xcomplex<float> > alm1_cxx(l_maxx, l_maxx);
  Alm< xcomplex<float> > alm2_cxx(l_maxx, l_maxx);
  Healpix_Map<float> map1_cxx((int)nside, RING, SET_NSIDE);
  Healpix_Map<float> map2_cxx((int)nside, RING, SET_NSIDE);
  
  long l, m, index;
  
  //-- Fill
  for (m=0; m<=l; m++) {
    for (l=m; l<=l_maxx; l++) {
      index = m*(2*l_maxx+1-m)/2+l;
      alm1_cxx(l, m).Set(alm[0+2*index], alm[1+2*index]);
      alm2_cxx(l, m).Set(0.0, 0.0); //-- No B modes
    }
  }
  
  alm2map_spin(alm1_cxx, alm2_cxx, map1_cxx, map2_cxx, 2); //-- spin = 2
  
  int64 i;
  //-- Retrieve
  for (i=0; i<12*nside*nside; i++) {
    map1[i] = (float)map1_cxx[i];
    map2[i] = (float)map2_cxx[i];
  }
  return;
}

//----------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

