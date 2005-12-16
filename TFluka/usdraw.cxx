#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TFluka.h"
#include <TLorentzVector.h>
#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Femfstk.h"  //(EMFSTK) fluka common
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

  if (icode/100 == 2) {
      for (Int_t npnw = EMFSTK.npstrt-1; npnw <= EMFSTK.npemf-1; npnw++) {
	  if (EMFSTK.iespak[npnw][mkbmx2-1] ==  TRACKR.ispusr[mkbmx2 - 1] ) {
	      EMFSTK.iespak[npnw][mkbmx2 - 2] = 1;
// Save properties at point where particle disappears in case this is only an interruption
	      TLorentzVector p;
	      gMC->TrackMomentum(p);
	      EMFSTK.espark[npnw][0] = xsco;               // x
	      EMFSTK.espark[npnw][1] = ysco;               // y
	      EMFSTK.espark[npnw][2] = zsco;               // z
	      EMFSTK.espark[npnw][3] = gMC->TrackTime();   // t
	      EMFSTK.espark[npnw][4] = p[0];               // px
	      EMFSTK.espark[npnw][5] = p[1];               // py
	      EMFSTK.espark[npnw][6] = p[2];               // pz
	      EMFSTK.espark[npnw][7] = p[3];               // e
	      EMFSTK.espark[npnw][8] = gMC->TrackLength(); // Length 
	  } // Track found in stack
      } // Loop over emf stack 
  } // Electromagnetic process
  


  fluka->SetMreg(mreg);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);

  if (debug) printf("USDRAW: Number of track segments:%6d %6d %6d %10.3e\n", TRACKR.ntrack, TRACKR.mtrack, icode, TRACKR.atrack);

  (TVirtualMCApplication::Instance())->Stepping();
  fluka->SetTrackIsNew(kFALSE);
  
} // end of usdraw
} // end of extern "C"

