

  /*******************************************************
   **  HEALPixFunctions.h				**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifdef __CAMELUS_USE_HEALPIX__
#ifndef __CAMELUS_HEALPIX_FUNCTIONS__
#define __CAMELUS_HEALPIX_FUNCTIONS__

#include <chealpix.h>

#ifdef __CAMELUS_USE_HEALPIX_CXX__
#include "wrapper.h"
#endif


//-- Functions related to HEALPix conversion
int compare(const void *a, const void *b);
void pixelToPatch(long nsidePat, long nsidePix, long pix, long *patch, int doNest);
void patchToPixels(long nsidePat, long nsidePix, long patch, hpint64 *pix, int doNest, int doSort);
void nsideToLevels(long nside, long *lenArr, long *cumLenArr);
void patchToCap(long nside, long patch, int doNest, int *cap);

//-- Functions related to Cartesian ordering
int log2_int(int a);
void ijPixToLocalNest(int resol, int i_pix, int j_pix, int *localNest);
void localNestToIJPix(int resol, int localNest, int *i_pix, int *j_pix);
void CartesianToLocalNest(int resol, int carte, int *localNest);
void localNestToCartesian(int resol, int localNest, int *carte);
void resolToIJPix(int resol, int *i_pix, int *j_pix);
void resolToLocalNest(int resol, int *localNest);
void CartesianToRingInPatch(long nside, long patch, int resol, long *pix);

//-- Functions related to neighbors
void decompose(long nside, long patch, int doNest, int *cap, int *level, int *length, int *off, int *j);
void decompose2(long nside, long patch, int doNest, long lenArr[], long cumLenArr[], long v1, long v2, int *cap, int *level, int *length, int *off, int *j);
void recombine(long nside, int level, int off, int j, long *patch);
void findNeighborsForCap(long nside, int level, int off, int j, long neiInfo[3][8]);
void findNeighborsForBelt(long nside, int level, int off, int j, long neiInfo[3][8]);
void findNeighbors(long nside, long patch, int doNest, long neighbors[8]);

//-- Functions related to boundary
void getLimits(long nside, long patch, int doNest, double limits[4]);
void adjustRA(double *RA, int N, int mode, int sign);
double HEALPixCapLine(double psi, double l);
double HEALPixBeltLine(double psi, double l);
double HEALPixCapArea(double psi_0, double psi, double l, double z_0);
double HEALPixCapInverseArea_plus(double A, double cst1, double cst2, double cst3);
double HEALPixCapInverseArea_minus(double A, double cst1, double cst2, double cst3);
double HEALPixBeltInverseArea(double A, double cst4, double cst5);
void capSampling(gsl_rng *generator, long nside, int N, int level, int length, int j, double z_0, double *pos1, double *pos2);
void beltSampling(gsl_rng *generator, long nside, int N, double *pos1, double *pos2);
void patchSampling(gsl_rng *generator, long nside, long patch, int N, int doNest, double *pos1, double *pos2);
void patchSampling2(gsl_rng *generator, long nside, long patch, int N, int doNest, double *pos1, double *pos2);
void patchSampling3(gsl_rng *generator, long nside, long patch, int N, int doNest, long lenArr[], long cumLenArr[], long v1, long v2, double *pos1, double *pos2);
void patchSampling4(gsl_rng *generator, long nside, int cap, int level, int length, int off, int j, double z0, double ctrPhi, double *pos1, double *pos2);

#endif
#endif

