#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TFluka.h"
#include <TLorentzVector.h>
#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#ifndef WIN32
# define usdraw usdraw_
#else
# define usdraw USDRAW
#endif
extern "C" {
void usdraw(Int_t& icode, Int_t& mreg, 
            Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  TFluka *fluka = (TFluka*)gMC;
  Int_t verbosityLevel = fluka->GetVerbosityLevel();
  Bool_t debug = (verbosityLevel >= 3)? kTRUE : kFALSE;
  fluka->SetCaller(6);
  fluka->SetIcode(icode);

  if (fluka->IsTrackDisappeared()) {
      TRACKR.ispusr[mkbmx2 - 2] = 1;
// Save properties at point where particle disappears in case this is only an interruption
      TLorentzVector p;
      gMC->TrackMomentum(p);
      
      TRACKR.spausr[0] = xsco;               // x
      TRACKR.spausr[1] = ysco;               // y
      TRACKR.spausr[2] = zsco;               // z
      TRACKR.spausr[3] = gMC->TrackTime();   // t
      TRACKR.spausr[4] = p[0];               // px
      TRACKR.spausr[5] = p[1];               // py
      TRACKR.spausr[6] = p[2];               // pz
      TRACKR.spausr[7] = p[3];               // e
      TRACKR.spausr[8] = gMC->TrackLength(); // Length 
  }

  fluka->SetMreg(mreg);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);

  if (debug) printf("USDRAW: Number of track segments:%6d %6d %6d %10.3e\n", TRACKR.ntrack, TRACKR.mtrack, icode, TRACKR.atrack);

  (TVirtualMCApplication::Instance())->Stepping();
  fluka->SetTrackIsNew(kFALSE);
  
} // end of usdraw
} // end of extern "C"

