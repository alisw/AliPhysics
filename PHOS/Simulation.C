#include "AliSimulation.h"
#include "TString.h"
#include "AliPHOSGetter.h"
#include "AliEMCALGetter.h"
#include "Riostream.h"

void simu(TString opt, TString name) 
{
  AliSimulation sim ; 
  // Generation and simulation
  if ( !opt.Contains("G") )
    sim.SetRunGeneration(kFALSE) ;
  // Making SDigits 
  if ( !opt.Contains("S") )
    sim.SetMakeSDigits("") ; 
  else 
    sim.SetMakeSDigits(name.Data()) ;
  // Making Digits 
  if ( !opt.Contains("D") )
    sim.SetMakeDigits("") ; 
  else 
    sim.SetMakeDigits(name.Data()) ;    
  //Merging
  // to implement 
  sim.Run() ;  
  // Checking result
  if ( name.Contains("PHOS") ) {
    cout << ">>>>>>>>>>>> PHOS " << endl ; 
    AliPHOSGetter * gime = AliPHOSGetter::Instance("galice.root") ;
    Int_t event ; 
    for (event = 0; event < gime->MaxEvent(); event++) {
      cout << "event # " << event << endl ; 
      gime->Event(event, "SD") ; 
      cout << "  SDigits # " << gime->SDigits()->GetEntries() << endl ; 
      cout << "   Digits # " << gime->Digits()->GetEntries() << endl ; 
    }
  }
  if ( name.Contains("EMCAL") ) {
    cout << ">>>>>>>>>>>> EMCAL " << endl ; 
    AliEMCALGetter * gime = AliEMCALGetter::Instance("galice.root") ;
    Int_t event ; 
    for (event = 0; event < gime->MaxEvent(); event++) {
      cout << "event # " << event << endl ; 
      gime->Event(event, "SD") ; 
      cout << "  SDigits # " << gime->SDigits()->GetEntries() << endl ; 
      cout << "   Digits # " << gime->Digits()->GetEntries() << endl ; 
    }
  }
}
