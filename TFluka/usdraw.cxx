#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TFluka.h"
#include "TFlukaCodes.h"
#include <TLorentzVector.h>
#include "Fdimpar.h"  //(DIMPAR) fluka include
#include "Ftrackr.h"  //(TRACKR) fluka common
#include "Fltclcm.h"  //(LTCLCM) fluka common
#include "Femfstk.h"  //(EMFSTK) fluka common
#include "Fflkstk.h"  //(FLKSTK) fluka common
#ifndef WIN32
# define usdraw usdraw_
#else
# define usdraw USDRAW
#endif

#include "TGeoManager.h" // <- delete

extern "C" {
void usdraw(Int_t& icode, Int_t& mreg, 
            Double_t& xsco, Double_t& ysco, Double_t& zsco)
{
  TFluka *fluka = (TFluka*)gMC;
  // nothing to do if particle inside dummy region
  if (mreg == fluka->GetDummyRegion()) return;
  Int_t verbosityLevel = fluka->GetVerbosityLevel();
  Bool_t debug = (verbosityLevel >= 3)? kTRUE : kFALSE;
  fluka->SetCaller(kUSDRAW);
  fluka->SetIcode((FlukaProcessCode_t) icode);

//
// Catch paused tracks
//
  if (icode/100 == kEMFSCO) {
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

  if (icode == kKASKADdray || icode == kKASKADpair || icode == kKASKADbrems) {
	TRACKR.ispusr[mkbmx2-2] = 1;
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
  
  //
  //
  Int_t mlttc = TRACKR.lt1trk; //LTCLCM.mlatm1;
  fluka->SetMreg(mreg, mlttc);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);

    // check region lattice consistency (debug Ernesto)
    // *****************************************************
    Int_t nodeId;
    Int_t volId = fluka->CurrentVolID(nodeId);
    Int_t crtlttc = gGeoManager->GetCurrentNodeId()+1;
    if(verbosityLevel>=3 && mreg != volId  && !gGeoManager->IsOutside() ) {
       cout << "  usdraw:   track=" << TRACKR.ispusr[mkbmx2-1] << " pdg=" << fluka->PDGFromId(TRACKR.jtrack)
            << " icode=" << icode << " gNstep=" << fluka->GetNstep() << endl
            << "               fluka   mreg=" << mreg << " mlttc=" << mlttc << endl
            << "               TGeo   volId=" << volId << " crtlttc=" << crtlttc << endl
            << "     common TRACKR   lt1trk=" << TRACKR.lt1trk << " lt2trk=" << TRACKR.lt2trk << endl
            << "     common LTCLCM   newlat=" << LTCLCM.newlat << " mlatld=" <<  LTCLCM.mlatld << endl
            << "                     mlatm1=" << LTCLCM.mlatm1 << " mltsen=" <<  LTCLCM.mltsen << endl
            << "                     mltsm1=" << LTCLCM.mltsm1 << " mlattc=" << LTCLCM.mlattc << endl;
        if( mlttc == crtlttc ) cout << "   *************************************************************" << endl;
    }
    // *****************************************************

  if (debug) printf("USDRAW: Number of track segments:%6d %6d icode=%d tof=%10.3e track=%d pdg=%d mreg=%d\n",
  TRACKR.ntrack, TRACKR.mtrack, icode, TRACKR.atrack, TRACKR.ispusr[mkbmx2-1], fluka->PDGFromId(TRACKR.jtrack), mreg );

  TVirtualMCStack* cppstack = fluka->GetStack();
  cppstack->SetCurrentTrack( TRACKR.ispusr[mkbmx2-1] );

  (TVirtualMCApplication::Instance())->Stepping();
  fluka->SetTrackIsNew(kFALSE);

 
} // end of usdraw
} // end of extern "C"

