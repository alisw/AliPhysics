#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TVirtualMCStack.h"

#ifndef WITH_ROOT
#include "TFluka.h"
#else
#include "TFlukaGeo.h"
#endif

// Fluka include
#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Fdblprc.h"  //(DBLPRC) fluka common
#include "Ftrackr.h"  //(TRACKR) fluka common

#ifndef WIN32
# define mgdraw mgdraw_
#else
# define mgdraw MGDRAW
#endif

extern "C" {
void mgdraw(Int_t& icode, Int_t& mreg)
{
    TFluka* fluka =  (TFluka*) gMC;
//    Int_t verbosityLevel = fluka->GetVerbosityLevel();
//
//  Make sure that stack has currrent track Id
    Int_t trackId = TRACKR.ispusr[mkbmx2-1];
    TVirtualMCStack* cppstack = fluka->GetStack();
    cppstack->SetCurrentTrack(trackId);
//
//    
    fluka->SetMreg(mreg);
    fluka->SetNewreg(mreg);
    fluka->SetIcode(icode);
    fluka->SetCaller(4);
    
//    if (verbosityLevel >= 3) {
//      cout << endl << " !!! I am in mgdraw - calling Stepping()" << endl;
//      cout << endl << " Track Id =" << trackId << endl;
//    }
    if (TRACKR.jtrack == -1)
    printf("Cerenkov photon: region# %6d icode %6d \n", mreg, icode);

    
    (TVirtualMCApplication::Instance())->Stepping();
    fluka->SetTrackIsNew(kFALSE);
} // end of mgdraw
} // end of extern "C"

