#include "AliPHOSPID.h"
#include "AliPHOSGetter.h"

void testreconRecParticles(Int_t nevent = 1, const char *config="testconfig.C")

{

AliPHOSGetter * gime = AliPHOSGetter::GetInstance("galice.root") ;
 AliPHOSPID * pid = new AliPHOSPIDv1("galice.root","test suite");
 pid->ExecuteTask("deb");
 //cout << "# of trackSegments " << gime->trackSegments()<<GetEntries()<< endl;
 cerr<<"_____________________________________________________________________"<<endl;
 cerr<<" "<<endl;
 cerr<<"      MESS ==> reconstruction ended successfully."<<endl;
 cerr<<"_____________________________________________________________________"<<endl;




}
