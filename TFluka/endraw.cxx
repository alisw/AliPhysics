#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TGeoMaterial.h"
#include "TGeoManager.h"
#include "TFlukaCerenkov.h"

#ifndef WITH_ROOT
#include "TFluka.h"
#else
#include "TFlukaGeo.h"
#endif

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
  fluka->SetCaller(3);
  fluka->SetIcode(icode);
  fluka->SetRull(rull);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);
  fluka->SetMreg(mreg);
  if (icode == 11) {
    cout << " For icode=" << icode << " Stepping is NOT called" << endl;
    return;
  }
  if (TRACKR.jtrack == -1) {
// Handle quantum efficiency the G3 way
      printf("endraw: Cerenkov photon depositing energy: %d %e\n", mreg, rull);
      TGeoMaterial* material = (gGeoManager->GetCurrentVolume())->GetMaterial();
      Int_t nmat = material->GetIndex();
      TFlukaCerenkov*  cerenkov = dynamic_cast<TFlukaCerenkov*> (material->GetCerenkovProperties());
      if (cerenkov) {
	  Double_t eff = (cerenkov->GetQuantumEfficiency(rull));
	  if (gRandom->Rndm() > eff) {
	      rull = 0.;
	      fluka->SetRull(rull);
	  }
      }
  }
  (TVirtualMCApplication::Instance())->Stepping();
  fluka->SetTrackIsNew(kFALSE);
} // end of endraw
} // end of extern "C"

