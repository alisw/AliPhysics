#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TVirtualMCStack.h"

#include "TFluka.h"

// Fluka include
#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Fdblprc.h"  //(DBLPRC) fluka common
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fopphst.h"  //(OPPHST) fluka common
#include "Fstack.h"   //(STACK)  fluka common

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
//
    Int_t trackId = -1;
    TVirtualMCStack* cppstack = fluka->GetStack();
    
    if (TRACKR.jtrack == -1) {
	trackId = OPPHST.louopp[OPPHST.lstopp];
	if (trackId == 0) {
	    trackId = STACK.ispark[STACK.lstack][mkbmx2-1];
	}
    } else {
	trackId = TRACKR.ispusr[mkbmx2-1];
    }
    
    cppstack->SetCurrentTrack(trackId);
//
//    
    fluka->SetMreg(mreg);
    fluka->SetNewreg(mreg);
    fluka->SetIcode(icode);
    fluka->SetCaller(4);

    if (!TRACKR.ispusr[mkbmx2 - 2]) {
	//
	// Single step
	if (verbosityLevel >= 3) {
	    cout << endl << " !!! I am in mgdraw - calling Stepping(): " << icode << endl;
	    cout << endl << " Track Id = " << trackId << " region = " << mreg << endl;
	}
	(TVirtualMCApplication::Instance())->Stepping();
	fluka->SetTrackIsNew(kFALSE);
    } else {
	//
	// Tracking is being resumed after secondary tracking
	//
	// Reset flag
	TRACKR.ispusr[mkbmx2 - 2] = 0;
//
	fluka->SetTrackIsNew(kTRUE);

	if (verbosityLevel >= 3) {
	    cout << endl << " !!! I am in mgdraw - resuming Stepping(): " << trackId << endl;
	}

	(TVirtualMCApplication::Instance())->Stepping();
	fluka->SetTrackIsNew(kFALSE);

	if (verbosityLevel >= 3) {
	    cout << endl << " !!! I am in mgdraw - first Stepping() after resume: " << icode << endl;
	    cout << endl << " Track Id = " << trackId << " region = " << mreg << endl;
	}
	(TVirtualMCApplication::Instance())->Stepping();
    }
    
    
} // end of mgdraw
} // end of extern "C"

