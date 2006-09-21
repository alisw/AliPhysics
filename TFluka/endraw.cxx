#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TGeoMaterial.h"
#include "TGeoManager.h"
#include <TParticle.h>
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
  Int_t mlttc = (icode==kKASKADinelarecoil) ? TRACKR.lt2trk : TRACKR.lt1trk; //LTCLCM.mlatm1;
  fluka->SetCaller(kENDRAW);
  fluka->SetRull(rull);
  fluka->SetXsco(xsco);
  fluka->SetYsco(ysco);
  fluka->SetZsco(zsco);
  fluka->SetMreg(mreg, mlttc);

  Float_t edep = rull;
  
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

  TVirtualMCStack* cppstack = fluka->GetStack();

  if (debug) {
     cout << "ENDRAW For icode=" << icode 
          << " track=" << TRACKR.ispusr[mkbmx2-1] << " pdg=" << fluka->PDGFromId(TRACKR.jtrack)
          << " edep="<< edep << " mreg=" << mreg << endl;
  }

    // check region lattice consistency (debug Ernesto)
    // *****************************************************
    Int_t nodeId;
    Int_t volId = fluka->CurrentVolID(nodeId);
    Int_t crtlttc = gGeoManager->GetCurrentNodeId()+1;
    if( mreg != volId  && !gGeoManager->IsOutside() ) {
       cout << "  endraw:   track=" << TRACKR.ispusr[mkbmx2-1] << " pdg=" << fluka->PDGFromId(TRACKR.jtrack)
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


  if (icode != kEMFSCOstopping1 && icode != kEMFSCOstopping2) {
      fluka->SetIcode((FlukaProcessCode_t)icode);
      fluka->SetRull(edep);
      if (icode == kKASKADelarecoil && TRACKR.ispusr[mkbmx2-5]) {
         //  Elastic recoil and in stuprf npprmr > 0,
         //  the secondary being loaded is actually still the interacting particle
         cppstack->SetCurrentTrack( TRACKR.ispusr[mkbmx2-4] );
         //      cout << "endraw elastic recoil track=" << TRACKR.ispusr[mkbmx2-1] << " parent=" << TRACKR.ispusr[mkbmx2-4]
         //           << endl;
      }
      else
          cppstack->SetCurrentTrack(TRACKR.ispusr[mkbmx2-1] );
      (TVirtualMCApplication::Instance())->Stepping();
  } else {
  //
  // For icode 21,22 the particle has fallen below thresshold.
  // This has to be signalled to the StepManager() 
  //
      cppstack->SetCurrentTrack( TRACKR.ispusr[mkbmx2-1] );
      fluka->SetRull(edep);
      fluka->SetIcode((FlukaProcessCode_t) icode);
      (TVirtualMCApplication::Instance())->Stepping();
      fluka->SetTrackIsNew(kFALSE);
      fluka->SetIcode((FlukaProcessCode_t)icode);
      fluka->SetRull(0.);
      (TVirtualMCApplication::Instance())->Stepping();
//      cppstack->SetCurrentTrack( saveTrackId );

  }
} // end of endraw
} // end of extern "C"

