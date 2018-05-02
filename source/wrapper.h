

  /*******************************************************
   **  wrapper.h					**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#ifdef __cplusplus
extern "C" {
#endif


//-- Functions related to a_lm
void readWeight(long nside, double *weight);
void mapToAlm(long nside, float *map, int l_maxx, double *alm, double *weight);
void almToMap(int l_maxx, double *alm, long nside, float *map);
void almToMap_spin2(int l_maxx, double *alm, long nside, float *map1, float *map2);

#ifdef __cplusplus
}
#endif

