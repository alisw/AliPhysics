#include "AliPHOSDigitizer.h"
#include "AliPHOSGetter.h"
#include "TSystem.h"
#include "AliPHOSDigit.h"


void testreconDigits(Int_t nevent = 1, const char *config="testconfig.C")
{ const Float_t maxDigits = 3489.41 ;
  const Float_t widDigits = TMath::Sqrt(maxDigits) ;
  AliPHOSDigitizer *d = new AliPHOSDigitizer("testPHOS.root","test suite");
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  d->ExecuteTask("deb"); 
  Float_t nDigits = (Float_t)(gime->Digitizer()->GetDigitsInRun()) / gime->MaxEvent();
  
  if ( nDigits < maxDigits-widDigits || nDigits > maxDigits+widDigits ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the Digits process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' schutz@in2p3.fr");
  }
 cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Digits process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;
}

