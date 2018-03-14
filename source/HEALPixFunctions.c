

  /*******************************************************
   **  HEALPixFunctions.c				**
   **  Version 2018.03.13				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "HEALPixFunctions.h"


//----------------------------------------------------------------------
//-- Functions related to HEALPix conversion

#ifdef __CAMELUS_USE_HEALPIX__
int compare(const void *a, const void *b)
{
  hpint64 a2 = *(hpint64*)a;
  hpint64 b2 = *(hpint64*)b;
  return (a2 > b2) - (a2 < b2);
}

void pixelToPatch(long nsidePat, long nsidePix, long pix, long *patch, int doNest)
{
  long pixNest = pix;
  if (doNest == 0) ring2nest(nsidePix, pix, &pixNest); //-- Require nest
  
  int resol  = nsidePix / nsidePat;
  int length = resol * resol;
  
  patch[0]  = pixNest / length;
  if (doNest == 0) nest2ring(nsidePat, patch[0], patch);
  return;
}

void patchToPixels(long nsidePat, long nsidePix, long patch, hpint64 *pix, int doNest, int doSort)
{
  //-- *pix should be initialized to have length = (nsidePix / nsidePat)^2
  
  long patNest = patch;
  if (doNest == 0) ring2nest(nsidePat, patch, &patNest); //-- Require nest
  
  int resol  = nsidePix / nsidePat;
  int length = resol * resol;
  
  int i;
  for (i=0; i<length; i++) pix[i] = patNest * length + i;
  
  if (doNest == 0) for (i=0; i<length; i++) nest2ring64(nsidePix, pix[i], &pix[i]); //-- Require nest
  if (doSort == 1) qsort(pix, length, sizeof(hpint64), compare);
  return;
}

void nsideToLevels(long nside, long *lenArr, long *cumLenArr)
{
  long beltLength     = 4 * nside;
  long nbLevelsInBelt = 2 * nside - 1;
  
  cumLenArr[0] = 0;
  int i, j = 1;
  
  for (i=1; i<=nside; i++, j++) {
    lenArr[j-1]  = 4 * i;
    cumLenArr[j] = cumLenArr[j-1] + lenArr[j-1];
  }
  
  for (i=0; i<nbLevelsInBelt; i++, j++) {
    lenArr[j-1]  = beltLength;
    cumLenArr[j] = cumLenArr[j-1] + lenArr[j-1];
  }
  
  for (i=nside; i>0; i--, j++) {
    lenArr[j-1]  = 4 * i;
    cumLenArr[j] = cumLenArr[j-1] + lenArr[j-1];
  }
  return;
}

void patchToCap(long nside, long patch, int doNest, int *cap)
{
  long patRing = patch;
  if (doNest == 1) nest2ring(nside, patch, &patRing); //-- Require ring
  
  long v1 = 2 * nside * (nside + 1);
  long v2 = 2 * nside * (5 * nside - 1);
  cap[0] = -1 + (int)(patRing < v2) + (int)(patRing < v1);
  return;
}

//----------------------------------------------------------------------
//-- Functions related to Cartesian ordering

int log2_int(int a)
{
  if (a < 1)  return -1;
  if (a == 1) return 0;
  return 1 + log2_int(a >> 1);
}

void ijPixToLocalNest(int resol, int i_pix, int j_pix, int *localNest)
{
  localNest[0] = 0;
  
  int halfNbBits = log2_int(resol);
  int i, k;
  
  for (i=0; i<halfNbBits; i++) {
    k = 1 << i;
    localNest[0] += (i_pix & k) * k;     //-- Bitwise and
    localNest[0] += (j_pix & k) * k * 2; //-- Bitwise and
  }
  return;
}

void localNestToIJPix(int resol, int localNest, int *i_pix, int *j_pix)
{
  i_pix[0] = 0;
  j_pix[0] = 0;
  
  int halfNbBits = log2_int(resol);
  int i, k;
  
  for (i=0; i<halfNbBits; i++) {
    k = 1 << i;
    i_pix[0] += (localNest & (k * k)) / k;
    j_pix[0] += (localNest & (2 * k * k)) / k;
  }
  j_pix[0] /= 2;
  return;
}

void CartesianToLocalNest(int resol, int carte, int *localNest)
{
  int i_pix = carte % resol;
  int j_pix = carte / resol;
  ijPixToLocalNest(resol, i_pix, j_pix, localNest);
  return;
}

void localNestToCartesian(int resol, int localNest, int *carte)
{
  int i_pix, j_pix;
  localNestToIJPix(resol, localNest, &i_pix, &j_pix);
  carte[0] = i_pix + j_pix * resol;
  return;
}

void resolToIJPix(int resol, int *i_pix, int *j_pix)
{
  int length = resol * resol;
  int i;
  for (i=0; i<length; i++) localNestToIJPix(resol, i, &i_pix[i], &j_pix[i]);
  return;
}

void resolToLocalNest(int resol, int *localNest)
{
  int length = resol * resol;
  int i;
  for (i=0; i<length; i++) CartesianToLocalNest(resol, i, &localNest[i]);
  return;
}

void CartesianToRingInPatch(long nside, long patch, int resol, long *pix)
{
  int length = resol * resol;
  long first;
  ring2nest(nside, patch, &first);
  first *= length;
  
  long nsidePix = nside * resol;
  int i, localNest;
  
  for (i=0; i<length; i++) {
    CartesianToLocalNest(resol, i, &localNest);
    pix[i] = first + localNest;
    nest2ring(nsidePix, pix[i], &pix[i]);
  }
  return;
}

//----------------------------------------------------------------------
//-- Functions related to neighbors

void decompose(long nside, long patch, int doNest, int *cap, int *level, int *length, int *off, int *j)
{
  long lenArr[4*nside-1];
  long cumLenArr[4*nside];
  
  nsideToLevels(nside, lenArr, cumLenArr);
  patchToCap(nside, patch, doNest, cap);
  
  long patRing = patch;
  if (doNest == 1) nest2ring(nside, patch, &patRing); //-- Require ring
  
  //-- Determine level
  level[0] = 0;
  while (patRing >= cumLenArr[level[0]+1]) level[0]++;
  
  //-- Determine length, off, j
  int rest = patRing - cumLenArr[level[0]];
  length[0] = lenArr[level[0]] / 4;
  off[0]    = (rest / length[0] + 2) % 4 - 2;
  j[0]      = rest % length[0];
  return;
}

void decompose2(long nside, long patch, int doNest, long lenArr[], long cumLenArr[], long v1, long v2, int *cap, int *level, int *length, int *off, int *j)
{
  long patRing = patch;
  if (doNest == 1) nest2ring(nside, patch, &patRing); //-- Require ring
  
  //-- Determine cap
  cap[0] = -1 + (int)(patRing < v2) + (int)(patRing < v1);
  
  //-- Determine level
  level[0] = 0;
  while (patRing >= cumLenArr[level[0]+1]) level[0]++;
  
  //-- Determine length, off, j
  int rest = patRing - cumLenArr[level[0]];
  length[0] = lenArr[level[0]] / 4;
  off[0]    = (rest / length[0] + 2) % 4 - 2;
  j[0]      = rest % length[0];
  return;
}

void recombine(long nside, int level, int off, int j, long *patch)
{
  if (level == -1) {
    patch[0] = -1;
    return;
  }
  
  long lenArr[4*nside-1];
  long cumLenArr[4*nside];
  
  nsideToLevels(nside, lenArr, cumLenArr);
  off      = (off + 4) % 4;
  patch[0] = cumLenArr[level] + off * lenArr[level] / 4 + j;
  return;
}

void findNeighborsForCap(long nside, int level, int off, int j, long neiInfo[3][8])
{
  //-- neiInfor contains level, length, off, j of 8 neighbors
  //-- in the order of SW, W, NW, N, NE, E, SE, S
  
      neiInfo[0][0] = level+1; neiInfo[1][0] = off;   neiInfo[2][0] = j;
      neiInfo[0][1] = level;   neiInfo[1][1] = off;   neiInfo[2][1] = j-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off;   neiInfo[2][2] = j-1;
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j-1;
      neiInfo[0][4] = level-1; neiInfo[1][4] = off;   neiInfo[2][4] = j;
      neiInfo[0][5] = level;   neiInfo[1][5] = off;   neiInfo[2][5] = j+1;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off;   neiInfo[2][6] = j+1;
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j+1;
  
  if (j == 0) {
      neiInfo[0][1] = level+1; neiInfo[1][1] = off-1; neiInfo[2][1] = level+1;
      neiInfo[0][2] = level;   neiInfo[1][2] = off-1; neiInfo[2][2] = level;
      neiInfo[0][3] = level-1; neiInfo[1][3] = off-1; neiInfo[2][3] = level-1;
  }
  
  if (j == level) {
      neiInfo[0][3] = level-1; neiInfo[1][3] = off+1; neiInfo[2][3] = 0;
      neiInfo[0][4] = level;   neiInfo[1][4] = off+1; neiInfo[2][4] = 0;
      neiInfo[0][5] = level+1; neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
  }
  
  if (level == nside - 1) {
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j;
    if (j == 0) {
      neiInfo[0][1] = -1;
    }
    if (j == level) {
      neiInfo[0][5] = -1;
      neiInfo[0][6] = nside;   neiInfo[1][6] = off+1; neiInfo[2][6] = 0;
    }
  }
  
  if (level == 0) {
      neiInfo[0][0] = 1;       neiInfo[1][0] = off;   neiInfo[2][0] = 0;
      neiInfo[0][1] = 1;       neiInfo[1][1] = off-1; neiInfo[2][1] = 1;
      neiInfo[0][2] = 0;       neiInfo[1][2] = off-1; neiInfo[2][2] = 0;
      neiInfo[0][3] = 0;       neiInfo[1][3] = off-2; neiInfo[2][3] = 0;
      neiInfo[0][4] = 0;       neiInfo[1][4] = off+1; neiInfo[2][4] = 0;
      neiInfo[0][5] = 1;       neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
      neiInfo[0][6] = 1;       neiInfo[1][6] = off;   neiInfo[2][6] = 1;
      neiInfo[0][7] = 2;       neiInfo[1][7] = off;   neiInfo[2][7] = 1;
  }
  
  int i;
  for (i=0; i<8; i++) neiInfo[1][i] = (neiInfo[1][i] + 4) % 4;
  return; 
}

void findNeighborsForBelt(long nside, int level, int off, int j, long neiInfo[3][8])
{
  //-- neiInfor contains level, length, off, j of 8 neighbors
  //-- in the order of SW, W, NW, N, NE, E, SE, S
  
  if (level % 2 == 0) {
      neiInfo[0][0] = level+1; neiInfo[1][0] = off;   neiInfo[2][0] = j-1;
      neiInfo[0][1] = level;   neiInfo[1][1] = off;   neiInfo[2][1] = j-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off;   neiInfo[2][2] = j-1;
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j;
      neiInfo[0][4] = level-1; neiInfo[1][4] = off;   neiInfo[2][4] = j;
      neiInfo[0][5] = level;   neiInfo[1][5] = off;   neiInfo[2][5] = j+1;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off;   neiInfo[2][6] = j;
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j;
  
    if (j == 0) {
      neiInfo[0][0] = level+1; neiInfo[1][0] = off-1; neiInfo[2][0] = nside-1;
      neiInfo[0][1] = level;   neiInfo[1][1] = off-1; neiInfo[2][1] = nside-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off-1; neiInfo[2][2] = nside-1;
    }
    
    if (j == nside - 1) {
      neiInfo[0][5] = level;   neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
    }
    
    if (level == nside) {
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j-1;
      if (j == 0) {
      neiInfo[0][3] = -1;
      }
    }
    if (level == 3*nside-2) {
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j-1;
      if (j == 0) {
      neiInfo[0][7] = -1;
      }
    }
  }
  
  else {
      neiInfo[0][0] = level+1; neiInfo[1][0] = off;   neiInfo[2][0] = j;
      neiInfo[0][1] = level;   neiInfo[1][1] = off;   neiInfo[2][1] = j-1;
      neiInfo[0][2] = level-1; neiInfo[1][2] = off;   neiInfo[2][2] = j;
      neiInfo[0][3] = level-2; neiInfo[1][3] = off;   neiInfo[2][3] = j;
      neiInfo[0][4] = level-1; neiInfo[1][4] = off;   neiInfo[2][4] = j+1;
      neiInfo[0][5] = level;   neiInfo[1][5] = off;   neiInfo[2][5] = j+1;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off;   neiInfo[2][6] = j+1;
      neiInfo[0][7] = level+2; neiInfo[1][7] = off;   neiInfo[2][7] = j;
      
    if (j == 0) {
      neiInfo[0][1] = level;   neiInfo[1][1] = off-1; neiInfo[2][1] = nside-1;
    }
      
    if (j == nside - 1) {
      neiInfo[0][4] = level-1; neiInfo[1][4] = off+1; neiInfo[2][4] = 0;
      neiInfo[0][5] = level;   neiInfo[1][5] = off+1; neiInfo[2][5] = 0;
      neiInfo[0][6] = level+1; neiInfo[1][6] = off+1; neiInfo[2][6] = 0;
    }
  }
  
  int i;
  for (i=0; i<8; i++) neiInfo[1][i] = (neiInfo[1][i] + 4) % 4;
  return; 
}

void findNeighbors(long nside, long patch, int doNest, long neighbors[8])
{
  long neiInfo[4][8];
  int cap, level, length, off, j, i;
  
  long patRing = patch;
  if (doNest == 1) nest2ring(nside, patch, &patRing); //-- Require ring
  decompose(nside, patRing, doNest, &cap, &level, &length, &off, &j);
  
  if      (cap > 0)  findNeighborsForCap(nside, level, off, j, neiInfo);
  else if (cap == 0) findNeighborsForBelt(nside, level, off, j, neiInfo);
  else  {
    level = 4*nside-2 - level;
    findNeighborsForCap(nside, level, off, j, neiInfo);
    for (i=0; i<8; i++) neiInfo[0][i] = (neiInfo[0][i] == -1) ? -1 : 4*nside-2 - neiInfo[0][i];
    for (i=0; i<3; i++) {
      neighbors[0]  = neiInfo[i][0];
      neiInfo[i][0] = neiInfo[i][2];
      neiInfo[i][2] = neighbors[0];
      neighbors[0]  = neiInfo[i][3];
      neiInfo[i][3] = neiInfo[i][7];
      neiInfo[i][7] = neighbors[0];
      neighbors[0]  = neiInfo[i][4];
      neiInfo[i][4] = neiInfo[i][6];
      neiInfo[i][6] = neighbors[0];
    }
  }
  
  for (i=0; i<8; i++) recombine(nside, neiInfo[0][i], neiInfo[1][i], neiInfo[2][i], &neighbors[i]);
  return; 
}

//----------------------------------------------------------------------
//-- Functions related to boundary

void getLimits(long nside, long patch, int doNest, double limits[4])
{
  //-- Return z_min, z_max, phi_min, phi_max
  
  int cap, level, length, off, j;
  
  long patRing = patch;
  if (doNest == 1) nest2ring(nside, patch, &patRing); //-- Require ring
  decompose(nside, patRing, doNest, &cap, &level, &length, &off, &j);
  
  if (cap == 1) {
    if (level == nside-1) limits[0] = 2.0/3.0 * (2.0 - (double)(level+2) / (double)nside);
    else                  limits[0] = 1.0 - pow((double)(level+2) / (double)nside, 2) / 3.0;
    limits[1] = 1.0 - pow((double)level / (double)nside, 2) / 3.0;
    limits[2] = HALF_PI * ((double)j / (double)length + off);
    limits[3] = HALF_PI * ((double)(j+1) / (double)length + off);
  }
  
  else if (cap == 0) {
    limits[0] = 2.0/3.0 * (2.0 - (double)(level+2) / (double)nside);
    limits[1] = 2.0/3.0 * (2.0 - (double)level / (double)nside);
    limits[2] = HALF_PI * (((double)j - 0.5 * (double)((level+1)%2)) / (double)nside + off);
    limits[3] = HALF_PI * (((double)j - 0.5 * (double)((level+1)%2) + 1) / (double)nside + off);
  }
  
  else {
    limits[0] = -(1.0 - pow((double)(4*nside-2 - level) / (double)nside, 2) / 3.0);
    if (level == 3*nside-1) limits[1] = 2.0/3.0 * (2.0 - (double)level / (double)nside);
    else                    limits[1] = -(1.0 - pow((double)(4*nside - level) / (double)nside, 2) / 3.0);
    limits[2] = HALF_PI * ((double)j / (double)length + off);
    limits[3] = HALF_PI * ((double)(j+1) / (double)length + off);
  }
  return;
}

void adjustRA(double *RA, int N, int mode, int sign)
{
  double half;
  if (mode == 1)      half = 180.0;
  else if (mode == 0) half = PI;
  else                half = 2.0;
  
  int i;
  if (sign > 0) {
    for (i=0; i<N; i++) RA[i] = fmod(RA[i], 2.0*half);
  }
  else if (sign < 0) {
    for (i=0; i<N; i++) RA[i] = fmod(RA[i], 2.0*half) - 2.0*half;
  }
  else {
    for (i=0; i<N; i++) RA[i] = fmod(RA[i] + half, 2.0*half) - half;
  }
  return;
}

double HEALPixCapLine(double psi, double l)
{
  return 1.0 - pow(l / psi, 2) / 3.0;
}

double HEALPixBeltLine(double psi, double l)
{
  return 4.0/3.0 * (0.5 - l + psi);
}

double HEALPixCapArea(double psi_0, double psi, double l, double z_0)
{
  return (1.0-z_0) * (psi-psi_0) + (l*l)/3.0 * (1.0/psi-1.0/psi_0);
}

double HEALPixCapInverseArea_plus(double A, double cst1, double cst2, double cst3)
{
  return (A + cst1 + sqrt(pow(A + cst1, 2) - cst2)) * cst3;
}

double HEALPixCapInverseArea_minus(double A, double cst1, double cst2, double cst3)
{
  return (A + cst1 - sqrt(pow(A + cst1, 2) - cst2)) * cst3;
}

double HEALPixBeltInverseArea(double A, double cst4, double cst5)
{
  return 0.75 * (-cst4 - sqrt(cst5 + 8.0/3.0*A));
}

void capSampling(gsl_rng *generator, long nside, int N, int level, int length, int j, double z_0, double *pos1, double *pos2)
{
  //-- Determine corners
  double dummy        = -1.0;
  double j2           = (double)j;
  double length2      = (double)length;
  double cornerPsi[4] = {dummy, j2/length2, (j2+1.0)/(length2+1.0), (j2+1.0)/length2};
  if (level > 0)         cornerPsi[0] = j2 / (length2-1.0);
  if (level == nside -1) cornerPsi[2] = 0.5 * (cornerPsi[1] + cornerPsi[3]);
  
  //-- Determine psi
  double psi0Arr[4]   = {cornerPsi[1], 1.0-cornerPsi[3], 1.0-cornerPsi[2], cornerPsi[2]};
  if (j == 0)     psi0Arr[0] = dummy;
  if (j == level) psi0Arr[1] = dummy;
  
  //-- Determine l
  double dl   = 1.0 / (double)nside;
  double lArr[4];
  lArr[0] = j2 * dl;
  lArr[1] = (length2-1.0-j2) * dl;
  lArr[2] = lArr[1] + dl;
  lArr[3] = lArr[0] + dl;
  
  double halfDenom = 3.0 * pow((double)nside, 2);
  double A_0       = 1.0 / (2.0 * halfDenom);
  double AArr[4]   = {A_0, A_0, A_0, A_0};
  
  //-- Determine areas
  if (level == 0 && nside > 1) {
    AArr[2] = -HEALPixCapArea(1.0-cornerPsi[2], 1.0-cornerPsi[1], lArr[2], z_0);
    AArr[3] = -HEALPixCapArea(    cornerPsi[2],     cornerPsi[3], lArr[3], z_0);
  }
  else {
    if (level < nside - 1) {
      AArr[2] = -HEALPixCapArea(1.0-cornerPsi[2], 1.0-cornerPsi[1], lArr[2], z_0);
      AArr[3] = -HEALPixCapArea(    cornerPsi[2],     cornerPsi[3], lArr[3], z_0);
    }
    
    if (j == 0) {
      AArr[0] = 0.0;
      AArr[1] = HEALPixCapArea(1.0-cornerPsi[3], 1.0-cornerPsi[0], lArr[1], z_0);
    }
    else if (j == level) {
      AArr[0] = HEALPixCapArea(    cornerPsi[1],     cornerPsi[0], lArr[0], z_0);
      AArr[1] = 0.0;
    }
    else {
      AArr[0] = HEALPixCapArea(    cornerPsi[1],     cornerPsi[0], lArr[0], z_0);
      AArr[1] = HEALPixCapArea(1.0-cornerPsi[3], 1.0-cornerPsi[0], lArr[1], z_0);
    }
  }
  
  double cst3 = 0.5 / (1.0 - z_0);
  double cst1[4], cst2[4], cst4[4], cst5[4];
  double ACum[5];
  int i;
  
  //-- Original sampling functions:
  //--   cSampP = lambda A, l, psi_0: (A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0) + np.sqrt((A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0))**2 - 4.0*(1-z_0)*l**2/3.0)) * 0.5 / (1.0 - z_0)
  //--   cSampM = lambda A, l, psi_0: (A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0) - np.sqrt((A +(1.0-z_0)*psi_0+l**2/(3.0*psi_0))**2 - 4.0*(1-z_0)*l**2/3.0)) * 0.5 / (1.0 - z_0)
  //-- Belows are simplification.
  ACum[0] = 0.0;
  for (i=0; i<4; i++) {
    ACum[i+1] = ACum[i] + AArr[i];
    cst1[i] = (1.0 - z_0) * psi0Arr[i] + pow(lArr[i], 2) / (3.0 * psi0Arr[i]);
    cst2[i] = 4.0/3.0 * (1 - z_0) * pow(lArr[i], 2);
    cst4[i] = 4.0/3.0 * (0.5 - lArr[i]) - z_0;
    cst5[i] = pow(cst4[i] + 4.0/3.0*psi0Arr[i], 2);
  }
  
  int count = 0;
  double p, q, rest, psi, z;
  
  while (count < N) {
    //-- Sample and determine which zone
    p = gsl_ran_flat(generator, 0.0, ACum[4]);
    q = gsl_ran_flat(generator, 0.0, 1.0);
    for (i=0; i<4; i++) {
      if (ACum[i+1] > p) break;
    }
    rest = p - ACum[i];
    
    //-- Transform into psi & z
    if (i < 2) {
      if (level == 0) {
	psi = rest * halfDenom;
	z   = q / halfDenom;
      }
      else {
	psi = HEALPixCapInverseArea_plus(rest, cst1[i], cst2[i], cst3);
	z   = (HEALPixCapLine(psi, lArr[i]) - z_0) * q;
      }
    }
    else {
      if (level == nside - 1) {
	psi = HEALPixBeltInverseArea(-rest, cst4[i], cst5[i]);
	z   = (HEALPixBeltLine(psi, lArr[i]) - z_0) * q;
      }
      else {
	psi = HEALPixCapInverseArea_minus(-rest, cst1[i], cst2[i], cst3);
	z   = (HEALPixCapLine(psi, lArr[i]) - z_0) * q;
      }
    }
    
    //-- Flip to the original referential
    pos1[count] = (i == 1 || i == 2) ? 1.0 - psi: psi;
    pos2[count] = z;
    count++;
  }
  return;
}

void beltSampling(gsl_rng *generator, long nside, int N, double *pos1, double *pos2)
{
  int count = 0;
  double psi_half = 0.5;
  double z_half   = 2.0 / 3.0;
  double factor   = 1.0 / (double)nside;
  double psi, z;
  
  while (count < N) {
    psi = gsl_ran_flat(generator, -psi_half, psi_half);
    z   = gsl_ran_flat(generator, 0.0, z_half);
    
    if (4.0 * fabs(psi) > 3.0 * z) {
      if (psi > 0.0) psi -= psi_half;
      else           psi += psi_half;
    }
    else z -= z_half;
    
    pos1[count] = psi * factor;
    pos2[count] = z * factor;
    count++;
  }
  return;
}

void patchSampling(gsl_rng *generator, long nside, long patch, int N, int doNest, double *pos1, double *pos2)
{
  int cap, level, length, off, j;
  double ctrTheta, ctrPhi, z0;
  
  //-- Decompose & get z0
  decompose(nside, patch, doNest, &cap, &level, &length, &off, &j);
  if (doNest == 1) pix2ang_nest(nside, patch, &ctrTheta, &ctrPhi);
  else             pix2ang_ring(nside, patch, &ctrTheta, &ctrPhi);
  z0 = cos(ctrTheta);
  
  int i;
  
  //-- Sample
  if (cap == 1) {
    capSampling(generator, nside, N, level, length, j, z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i]  = (pos1[i] + off) * HALF_PI;
      pos2[i] += z0;
    }
  }
  else if (cap == 0) {
    beltSampling(generator, nside, N, pos1, pos2);
    adjustRA(&ctrPhi, 1, 0, 0); //-- N = 1, mode = 0, sign = 0
    
    for (i=0; i<N; i++) {
      pos1[i]  = pos1[i] * HALF_PI + ctrPhi;
      pos2[i] += z0;
    }
  }
  else {
    capSampling(generator, nside, N, 4*nside-2-level, length, length-1-j, -z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i] = (1.0 - pos1[i] + off) * HALF_PI;
      pos2[i] = z0 - pos2[i];
    }
  }
  
  //-- Conversion
  for (i=0; i<N; i++) pos2[i]  = HALF_PI - acos(pos2[i]);
  return;
}

void patchSampling2(gsl_rng *generator, long nside, long patch, int N, int doNest, double *pos1, double *pos2)
{
  double limits[4];
  int count = 0;
  double z, theta, phi;
  long patch2;
  
  getLimits(nside, patch, doNest, limits);
  
  while (count < N) {
    z     = gsl_ran_flat(generator, limits[0], limits[1]);
    phi   = gsl_ran_flat(generator, limits[2], limits[3]);
    theta = acos(z);
    if (doNest == 1) ang2pix_nest(nside, theta, phi, &patch2);
    else             ang2pix_ring(nside, theta, phi, &patch2);
    
    if (patch2 == patch) {
      pos2[count] = HALF_PI - theta;
      count++;
    }
  }
  return;
}

void patchSampling3(gsl_rng *generator, long nside, long patch, int N, int doNest, long lenArr[], long cumLenArr[], long v1, long v2, double *pos1, double *pos2)
{
  int cap, level, length, off, j;
  double ctrTheta, ctrPhi, z0;
  
  //-- Decompose & get z0
  decompose2(nside, patch, doNest, lenArr, cumLenArr, v1, v2, &cap, &level, &length, &off, &j);
  if (doNest == 1) pix2ang_nest(nside, patch, &ctrTheta, &ctrPhi);
  else             pix2ang_ring(nside, patch, &ctrTheta, &ctrPhi);
  z0 = cos(ctrTheta);
  
  int i;
  
  //-- Sample
  if (cap == 1) {
    capSampling(generator, nside, N, level, length, j, z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i]  = (pos1[i] + off) * HALF_PI;
      pos2[i] += z0;
    }
  }
  else if (cap == 0) {
    beltSampling(generator, nside, N, pos1, pos2);
    adjustRA(&ctrPhi, 1, 0, 0); //-- N = 1, mode = 0, sign = 0
    
    for (i=0; i<N; i++) {
      pos1[i]  = pos1[i] * HALF_PI + ctrPhi;
      pos2[i] += z0;
    }
  }
  else {
    capSampling(generator, nside, N, 4*nside-2-level, length, length-1-j, -z0, pos1, pos2);
    
    for (i=0; i<N; i++) {
      pos1[i] = (1.0 - pos1[i] + off) * HALF_PI;
      pos2[i] = z0 - pos2[i];
    }
  }
  
  //-- Conversion
  for (i=0; i<N; i++) pos2[i] = HALF_PI - acos(pos2[i]);
  return;
}

void patchSampling4(gsl_rng *generator, long nside, int cap, int level, int length, int off, int j, double z0, double ctrPhi, double *pos1, double *pos2)
{
  //-- Sample, N = 1
  if (cap == 1) {
    capSampling(generator, nside, 1, level, length, j, z0, pos1, pos2);
    
    pos1[0]  = (pos1[0] + off) * HALF_PI;
    pos2[0] += z0;
  }
  else if (cap == 0) {
    beltSampling(generator, nside, 1, pos1, pos2);
    adjustRA(&ctrPhi, 1, 0, 0); //-- N = 1, mode = 0, sign = 0
    
    pos1[0]  = pos1[0] * HALF_PI + ctrPhi;
    pos2[0] += z0;
  }
  else {
    capSampling(generator, nside, 1, 4*nside-2-level, length, length-1-j, -z0, pos1, pos2);
    
    pos1[0] = (1.0 - pos1[0] + off) * HALF_PI;
    pos2[0] = z0 - pos2[0];
  }
  
  //-- Conversion
  pos2[0] = HALF_PI - acos(pos2[0]);
  return;
}
#endif

//----------------------------------------------------------------------

