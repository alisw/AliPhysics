#include "TFile.h"
#include "AliPHOSGetter.h"
#include "Riostream.h"
#include "AliESD.h"

void Ana() 
{
  AliPHOSGetter * gime = AliPHOSGetter::Instance("galice.root") ; 
  Int_t nEvent = gime->MaxEvent() ;  
  Int_t event ; 
  AliESD * esd ;
  for (event = 0 ; event < nEvent; event++) {
    esd = gime->ESD(event) ; 
    esd->Print(); 
    
  }
}
