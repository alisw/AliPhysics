#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TVirtualMCStack.h"
#include "TFluka.h"
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
    Int_t verbosityLevel = fluka->GetVerbosityLevel();
//
//  Make sure that stack has currrent track Id
    Int_t trackId = TRACKR.ispusr[mkbmx2-1];
    TVirtualMCStack* cppstack = fluka->GetStack();
    cppstack->SetCurrentTrack(trackId);
//
//    
    Int_t oldreg = ((TFluka*) gMC)->GetMreg();
    if (oldreg != mreg) {
//
//  Boundary Crossing
//
	fluka->SetNewreg(mreg);
	if (oldreg == -1) fluka->SetMreg(mreg);
	if (verbosityLevel >= 3)
	    printf("Boundary Crossing %d %d \n", oldreg, mreg);
    } else {
	fluka->SetMreg(mreg);
	fluka->SetNewreg(mreg);
	if (verbosityLevel >= 3)
	    printf("Normal step %d %d \n", oldreg, mreg);
    }
    fluka->SetIcode(icode);
    if (verbosityLevel >= 3) {
	cout << endl << " !!! I am in mgdraw - calling Stepping()" << endl;
	cout << endl << " Track Id =" << trackId << endl;
    }
    
    fluka->FutoTest();

    if (oldreg != mreg) {
//
//  Double step for boundary crossing
//
	fluka->SetTrackIsExiting();
	(TVirtualMCApplication::Instance())->Stepping();
	fluka->SetMreg(mreg);
	fluka->SetTrackIsEntering();
	(TVirtualMCApplication::Instance())->Stepping();
	fluka->SetTrackIsInside();
    } else {
	(TVirtualMCApplication::Instance())->Stepping();
    }
} // end of mgdraw
} // end of extern "C"

