#include "AliReconstruction.h"
#include "TString.h"
#include "Riostream.h"
#include "AliPHOSGetter.h"
#include "AliEMCALGetter.h"

void reco(TString opt, TString name) 
{
  AliReconstruction rec ; 
  if ( !opt.Contains("T") )
    rec.SetRunTracking(kFALSE) ;
  if ( opt.Contains("R") ) 
    rec.SetRunReconstruction(name.Data()) ; 
  if ( !opt.Contains("E") )
    rec.SetFillESD("") ; 
  else 
    rec.SetFillESD(name.Data()) ; 
  rec.Run() ;

  if ( name.Contains("PHOS") ) {
    cout << ">>>>>>>>>>>> PHOS " << endl ; 
    AliPHOSGetter * gime = AliPHOSGetter::Instance("galice.root") ; 
    Int_t event ; 
    for (event = 0; event < gime->MaxEvent(); event++) {
      cout << "event # " << event << endl ; 
      gime->Event(event, "RP") ; 
      cout << "   EMC RecPoints  # " << gime->EmcRecPoints()->GetEntries() << endl ; 
      cout << "   CPV RecPoints  # " << gime->CpvRecPoints()->GetEntries() << endl ; 
      cout << "   Track Segments # " << gime->TrackSegments()->GetEntries() << endl ; 
      cout << "   Rec Particles  # " << gime->RecParticles()->GetEntries() << endl ; 
    }
  } 
 if ( name.Contains("EMCAL") ) {
    cout << ">>>>>>>>>>>> EMCAL " << endl ; 
    AliEMCALGetter * gime = AliEMCALGetter::Instance("galice.root") ; 
    Int_t event ; 
    for (event = 0; event < gime->MaxEvent(); event++) {
      cout << "event # " << event << endl ; 
      gime->Event(event, "RP") ; 
      cout << "       RecPoints  # " << gime->ECARecPoints()->GetEntries() << endl ; 
      cout << "   Rec Particles  # " << gime->RecParticles()->GetEntries() << endl ; 
    }
 } 
}   
