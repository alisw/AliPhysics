#include "AliPHOSPID.h"
#include "AliPHOSGetter.h"
#include "AliPHOSPIDv1.h"
#include "TSystem.h"

void testreconRecParticles(Int_t nevent = 1, const char *config="testconfig.C")

{ 

  const Float_t maxRecParticles = 1 ;
  const Float_t widRecParticles = TMath::Sqrt(maxRecParticles) ;
  TString name = "test suite" ;


  AliPHOSPID * pid = new AliPHOSPIDv1("testPHOS.root",name.Data());
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  pid->ExecuteTask("deb");   
  Float_t nRecParticles =  (Float_t) (gime->PID(name.Data())->GetRecParticlesInRun())/gime->MaxEvent();
 
 
 if ( nRecParticles < maxRecParticles-0.25 || nRecParticles > maxRecParticles+0.25 ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the RecParticles process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' schutz@in2p3.fr");
 }
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> RecParticles process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;

 cerr<<"_____________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> reconstruction ended successfully."<<endl;
  cerr<<"_____________________________________________________________________"<<endl;

}
