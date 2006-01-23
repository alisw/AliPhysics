#include <Riostream.h>

#include "TFluka.h"
#include "TFlukaCodes.h"
#include "TFlukaMCGeometry.h"
#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fltclcm.h"  //(LTCLCM) fluka common
#ifndef WIN32
# define bxdraw bxdraw_
#else
# define bxdraw BXDRAW
#endif
extern "C" {
void bxdraw(Int_t& icode, Int_t& mreg, Int_t& newreg,
            Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
    TFluka* fluka = (TFluka*) gMC;
    Int_t oldlttc = LTCLCM.mlatm1;
    Int_t newlttc = LTCLCM.newlat;
    fluka->SetIcode((FlukaProcessCode_t)icode);
    fluka->SetNewreg(newreg,newlttc);
    fluka->SetXsco(xsco);
    fluka->SetYsco(ysco);
    fluka->SetZsco(zsco);
    Int_t verbosityLevel = fluka->GetVerbosityLevel();
    Bool_t debug = (verbosityLevel>=3)?kTRUE:kFALSE;
//
// Double step for boundary crossing
//
    fluka->SetTrackIsNew(kFALSE); // has to be called BEFORE Stepping()
    if (mreg != fluka->GetDummyRegion()) {
       if (debug) printf("bxdraw (ex) \n");
       fluka->SetTrackIsExiting();
       fluka->SetCaller(kBXExiting);
       fluka->SetMreg(mreg,oldlttc);
       (TVirtualMCApplication::Instance())->Stepping(); 
    }
    if (newreg != fluka->GetDummyRegion()) {
       if (debug) printf("bxdraw (en) \n");
       fluka->SetCaller(kBXEntering);
       fluka->SetTrackIsEntering();
       if (fluka->GetDummyBoundary() == 1) fluka->SetDummyBoundary(2);
       fluka->SetMreg(newreg,newlttc);
       (TVirtualMCApplication::Instance())->Stepping();
    }      

} // end of bxdraw
} // end of extern "C"

