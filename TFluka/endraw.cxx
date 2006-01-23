#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TGeoMaterial.h"
#include "TGeoManager.h"
#include "TFlukaCerenkov.h"

#include "TFluka.h"
#include "TFlukaCodes.h"

#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fltclcm.h"  //(LTCLCM) fluka common
#include "Fpaprop.h"  //(PAPROP) fluka common
#ifndef WIN32
# define endraw endraw_
#else
# define endraw ENDRAW
#endif
extern "C" {
void endraw(Int_t& icode, Int_t& mreg, Double_t& rull, Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  TFluka* fluka = (TFluka*) gMC;
  // nothing to do if particle in dummy region
  if (mreg == fluka->GetDummyRegion()) return;
  Int_t verbosityLevel = fluka->GetVerbosityLevel();
  Bool_t debug = (verbosityLevel >= 3)? kTRUE : kFALSE;
  Int_t mlttc = LTCLCM.mlatm1;
  fluka->SetCaller(kENDRAW);
  fluka->SetRull(rull);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);
  fluka->SetMreg(mreg, mlttc);
  
  Float_t edep = rull;
  
  if (icode == kKASKADinelarecoil) {
    if (debug) cout << " For icode=" << icode << " Stepping is NOT called" << endl;
    return;
  }

  if (TRACKR.jtrack == -1) {
// Handle quantum efficiency the G3 way
      if (debug) printf("endraw: Cerenkov photon depositing energy: %d %e\n", mreg, rull);
      TGeoMaterial* material = (gGeoManager->GetCurrentVolume())->GetMaterial();
      TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
      if (cerenkov) {
	  Double_t eff = (cerenkov->GetQuantumEfficiency(rull));
	  if (gRandom->Rndm() > eff) {
	      edep = 0.;
	  }
      }
  }

  if (icode != kEMFSCOstopping1 && icode != kEMFSCOstopping2) {
      fluka->SetIcode((FlukaProcessCode_t)icode);
      fluka->SetRull(edep);
      (TVirtualMCApplication::Instance())->Stepping();
  } else {
  //
  // For icode 21,22 the particle has fallen below thresshold.
  // This has to be signalled to the StepManager() 
  //
      fluka->SetRull(edep);
      fluka->SetIcode((FlukaProcessCode_t) icode);
      (TVirtualMCApplication::Instance())->Stepping();
      fluka->SetTrackIsNew(kFALSE);
      fluka->SetIcode((FlukaProcessCode_t)icode);
      fluka->SetRull(0.);
      (TVirtualMCApplication::Instance())->Stepping();
  }
} // end of endraw
} // end of extern "C"

