

  /*******************************
   **  main.h			**
   **  Chieh-An Lin		**
   **  Version 2016.03.20	**
   *******************************/


#ifndef __main__
#define __main__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"
#include "constraint.h"
#include "ABC.h"

#ifndef __releaseMenu__
  #include <mpi.h>
  #include "FITSFunctions.h"
  #include "paperIII.h"
#endif


int main(int argc, char *argv[]);
void printInstructions(int task, int printHeader);
void MPI_terminate(int MPISize, int MPIInd);
void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err);


#endif

