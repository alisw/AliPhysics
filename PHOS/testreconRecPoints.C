#include "AliPHOSClusterizer.h"
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

void testreconRecPoints(Int_t nevent = 1, const char *config="testconfig.C")
{
  const Float_t maxRecPoints = 223.26 ;
  const Float_t widRecPoints = TMath::Sqrt(maxRecPoints) ;
  TString name = "test suite" ;
  AliPHOSClusterizer * cluster = new  AliPHOSClusterizerv1("galice.root", name.Data());
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(); 
  cluster->ExecuteTask("deb");
  TString fullName = name + cluster->Version() ;  
  Float_t nRecPoints =  (Float_t) (gime->Clusterizer(fullName.Data())->GetRecPointsInRun()) / gime->MaxEvent();
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"nRecPoints vaut "<<nRecPoints<<endl;
  cerr<<"__________________________________________________________________"<<endl; 
   if ( nRecPoints < maxRecPoints-widRecPoints || nRecPoints > maxRecPoints+widRecPoints ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"       MESS ==> Error detected in the Clusterizing process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/galice.root galice.root | mail -s "PHOS INSTALLATION ERROR" schutz@in2p3.fr");
  }

 
}
