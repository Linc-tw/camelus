

  /*******************************************************
   **  main.h						**
   **  Version 2018.03.11				**
   **							**
   **  Copyright (C) 2018 - Chieh-An Lin		**
   **  GNU GPLv3 - https://www.gnu.org/licenses/	**
   *******************************************************/


#include "commonHeader.h"

#ifndef __CAMELUS_MAIN__
#define __CAMELUS_MAIN__

#ifdef __CAMELUS_USE_FITS__
#include "FITSFunctions.h"
#endif

#ifdef __CAMELUS_USE_HEALPIX__
#include "HEALPixFunctions.h"
#endif

#ifdef __CAMELUS_USE_HEALPIX_CXX__
#include "wrapper.h"
#endif

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
#include "ABC.h"

#ifdef __CAMELUS_BETA_MODE__
#include "peakVI.h"
#endif

#ifdef __CAMELUS_USE_LHF__
#include "LHF.h"
#endif


int main(int argc, char *argv[]);
void printInstructions(int task, int doHelp);
void printDetails(int task, int doHeader, int doHelp);
void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err);

#endif

