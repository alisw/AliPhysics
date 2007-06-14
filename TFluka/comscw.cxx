#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Fdblprc.h"  //(DBLPRC) fluka include
#include "Fiounit.h"  //(IOUNIT) fluka common
#include "Fflkmat.h"  //(FLKMAT) fluka common
#include "Fscohlp.h"  //(SCOHLP) fluka common
#include "Fsouevt.h"  //(SOUEVT) fluka common
#ifndef WIN32
#define comscw comscw_
#define type_of_call
#else
#define comscw COMSCW 
#define type_of_call  _stdcall
#endif
const Double_t  oned  = 1.00000e+00;
const Double_t  factr = 1.60217e-07;
extern "C" double type_of_call comscw(int& ij, double& xa, double& ya, double& za,
				    int& mreg, double& rull, int& llo, int& icall)

{
//
//  Return comscw 
//
   SCOHLP.lsczer = kFALSE;
   if (FLKMAT.rho[FLKMAT.medium[mreg]] <  1.0e-15) 
{   
   return 1.;
}
   if (SCOHLP.iscrng == 1 ) {
//    printf("comscw  %f %f %f %f %f \n", FLKMAT.rho[FLKMAT.medium[mreg]], 4.1e+15*rull*factr/FLKMAT.rho[FLKMAT.medium[mreg]], xa, ya, za);
      return factr/FLKMAT.rho[FLKMAT.medium[mreg]];
}
      return 1.;
}
