#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TGeoMaterial.h"
#include "TGeoManager.h"
#include "TFlukaCerenkov.h"

#include "TFluka.h"

#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#ifndef WIN32
# define endraw endraw_
#else
# define endraw ENDRAW
#endif
extern "C" {
void endraw(Int_t& icode, Int_t& mreg, Double_t& rull, Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  TFluka* fluka = (TFluka*) gMC;
  Int_t verbosityLevel = fluka->GetVerbosityLevel();
  Bool_t debug = (verbosityLevel >= 3)? kTRUE : kFALSE;
  fluka->SetCaller(3);
  fluka->SetRull(rull);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);
  fluka->SetMreg(mreg);
  
  Float_t edep = rull;
  
  if (icode == 11) {
    if (debug) cout << " For icode=" << icode << " Stepping is NOT called" << endl;
    return;
  }
  if (TRACKR.jtrack == -1) {
// Handle quantum efficiency the G3 way
      if (debug) printf("endraw: Cerenkov photon depositing energy: %d %e\n", mreg, rull);
      TGeoMaterial* material = (gGeoManager->GetCurrentVolume())->GetMaterial();
      // Int_t nmat = material->GetIndex();
      TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
      if (cerenkov) {
	  Double_t eff = (cerenkov->GetQuantumEfficiency(rull));
	  if (gRandom->Rndm() > eff) {
	      edep = 0.;
	  }
      }
  }
  if (debug) printf("endraw: Depositing energy for : %d %e icode: %d \n", TRACKR.ispusr[mkbmx2-1], rull, icode);

  if (icode != 21 && icode != 22) {
      fluka->SetIcode(icode);
      fluka->SetRull(edep);
      (TVirtualMCApplication::Instance())->Stepping();
  } else {
  //
  // for icode 21,22 the particle has fallen below thresshold 
  // This has to be signalled to the StepManager() 
  //
      fluka->SetRull(edep);
      fluka->SetIcode(20);
      (TVirtualMCApplication::Instance())->Stepping();
      fluka->SetTrackIsNew(kFALSE);
      fluka->SetIcode(icode);
      fluka->SetRull(0.);
      (TVirtualMCApplication::Instance())->Stepping();
  }
} // end of endraw
} // end of extern "C"

