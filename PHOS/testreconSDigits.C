#include "AliPHOSGetter.h"
#include "TSystem.h"
#include "AliPHOSDigit.h"
#include "AliPHOSSDigitizer.h"

void testreconSDigits(Int_t nevent = 1, const char *config="testconfig.C")
{
  cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr<<"           MESS ==> Beginning of the PHOS reconstruction. "<<endl;
  cerr<<" ___________________________________________________________________ "<<endl;
  const Float_t maxSDigits = 62.89 ;
  const Float_t widSDigits = TMath::Sqrt(maxSDigits) ;
  AliPHOSSDigitizer *sd = new AliPHOSSDigitizer("testPHOS.root","test suite");
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  sd->ExecuteTask("deb"); 
  Float_t nSDigits =  (Float_t) (gime->SDigitizer()->GetSDigitsInRun()) / gime->MaxEvent();
  
 
   if ( nSDigits < maxSDigits-widSDigits || nSDigits > maxSDigits+widSDigits ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the SDigits process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' schutz@in2p3.fr");
 }
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> SDigits process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;
}
