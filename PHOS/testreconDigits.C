#include "AliPHOSDigitizer.h"
#include "AliPHOSGetter.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "AliPHOSHit.h"
#include "TFolder.h"
#include "TStopwatch.h"
#include "TObjArray.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSSDigitizer.h"

void testreconDigits(Int_t nevent = 1, const char *config="testconfig.C")
{ const Float_t maxDigits = 0 ;
  const Float_t widDigits = TMath::Sqrt(maxDigits) ;
  AliPHOSDigitizer *d = new AliPHOSDigitizer("galice.root","test suite");
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  d->ExecuteTask("deb"); 
  Float_t nDigits = (Float_t)(gime->Digitizer()->GetDigitsInRun()) / gime->MaxEvent();
  // cout << "# of Digits " << gime->Digits()->GetEntries() << endl ;



}
