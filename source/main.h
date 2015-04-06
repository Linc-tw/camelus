

  /*******************************
   **  main.h			**
   **  Chieh-An Lin		**
   **  Version 2015.04.05	**
   *******************************/


#ifndef __main__
#define __main__

#include "commonHeader.h"
#include "peakParameters.h"
#include "haloSampling.h"
#include "rayTracing.h"
#include "smoothing.h"
#include "peakSelection.h"
#include "ABC.h"

#define __releaseMenu__
//#define __completeMenu__

//#include "FSL10.h"
//#include "paperI.h"
//#include "paperII.h"


int main(int argc, char *argv[]);
void sandbox(cosmo_hm *cmhm, peak_param *peak, error **err);


#endif