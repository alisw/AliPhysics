#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TFluka.h"
#ifndef WIN32
# define mgdraw mgdraw_
#else
# define mgdraw MGDRAW
#endif

extern "C" {
void mgdraw(Int_t& icode, Int_t& mreg)
{
    Int_t oldreg = ((TFluka*) gMC)->GetMreg();
    if (oldreg != mreg) {
//
//  Boundary Crossing
//
	((TFluka*) gMC)->SetNewreg(mreg);
	if (oldreg == -1) 
	    ((TFluka*) gMC)->SetMreg(mreg);
	printf("Boundary Crossing %d %d \n", oldreg, mreg);
    } else {
	((TFluka*) gMC)->SetMreg(mreg);
	((TFluka*) gMC)->SetNewreg(mreg);
	printf("Normal step %d %d \n", oldreg, mreg);
    }
    ((TFluka*) gMC)->SetIcode(icode);
    cout << endl << " !!! I am in mgdraw - calling Stepping()" << endl;
    ((TFluka*) gMC)->FutoTest();

    if (oldreg != mreg) {
//
//  Double step for boundary crossing
//
	((TFluka*) gMC)->SetTrackIsExiting();
	(TVirtualMCApplication::Instance())->Stepping();
	((TFluka*) gMC)->SetMreg(mreg);
	((TFluka*) gMC)->SetTrackIsEntering();
	(TVirtualMCApplication::Instance())->Stepping();
	((TFluka*) gMC)->SetTrackIsInside();
    } else {
	(TVirtualMCApplication::Instance())->Stepping();
    }
} // end of mgdraw
} // end of extern "C"

