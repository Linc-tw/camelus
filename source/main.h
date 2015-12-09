

  /*******************************
   **  main.h			**
   **  Chieh-An Lin		**
   **  Version 2015.12.09	**
   *******************************/


#ifndef __main__
#define __main__

#include <mpi.h>

#include "commonHeader.h"
#include "FITSFunctions.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"
#include "constraint.h"
#include "ABC.h"

#define __releaseMenu__

//#include "FSL10.h"
//#include "paperI.h"
//#include "paperII.h"
//#include "paperIII.h"


int main(int argc, char *argv[]);
void printInstructions(int task, int printHeader);
void MPI_finalize(int MPISize, int MPIInd);
void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err);


#endif

