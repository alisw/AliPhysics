#include "AliPHOSClusterizer.h"
#include "AliPHOSGetter.h"
#include "AliPHOSClusterizerv1.h"
#include "TSystem.h"


void testreconRecPoints(Int_t nevent = 1, const char *config="testconfig.C")
{
  const Float_t maxRecPoints = 222.83 ;
  const Float_t widRecPoints = TMath::Sqrt(maxRecPoints) ;
  TString name = "test suite" ;
  AliPHOSClusterizer * cluster = new  AliPHOSClusterizerv1("testPHOS.root", name.Data());
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance();
  cluster->ExecuteTask("deb");
  Float_t nRecPoints =  (Float_t) (gime->Clusterizer(name.Data())->GetRecPointsInRun()) / gime->MaxEvent();
 
   if ( nRecPoints < maxRecPoints-widRecPoints || nRecPoints > maxRecPoints+widRecPoints ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the Clusterizing process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' schutz@in2p3.fr");
 }
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Cluster process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;

 
}
