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
{ const Float_t maxDigits = 3483.41 ;
  const Float_t widDigits = TMath::Sqrt(maxDigits) ;
  AliPHOSDigitizer *d = new AliPHOSDigitizer("galice.root","test suite");
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  d->ExecuteTask("deb"); 
  Float_t nDigits = (Float_t)(gime->Digitizer()->GetDigitsInRun()) / gime->MaxEvent();
   cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"la valeur de nDigits est "<<nDigits<<endl;
  cerr<<"__________________________________________________________________"<<endl;
  if ( nDigits < maxDigits-widDigits || nDigits > maxDigits+widDigits ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"       MESS ==> Error detected in the Digits process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/galice.root galice.root | mail -s "PHOS INSTALLATION ERROR" schutz@in2p3.fr");
  }

}
