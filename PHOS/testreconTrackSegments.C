#include " AliPHOSTrackSegmentMaker.h"
#include "AliPHOSGetter.h"
#include "TSystem.h"

void testreconTrackSegments(Int_t nevent = 1, const char *config="testconfig.C")

{ 
  const Float_t maxTrackSegments = 1 ;
  const Float_t widTrackSegments = TMath::Sqrt(maxTrackSegments) ;
  TString name = "test suite" ;
 
  AliPHOSTrackSegmentMaker * tracks = new AliPHOSTrackSegmentMakerv1("galice.root",name.Data());
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  tracks->ExecuteTask("deb");
  TString fullName = name + tracks->Version() ;
  Float_t nTrackSegments =  (Float_t) (gime->TrackSegmentMaker(fullName.Data())->GetTrackSegmentsInRun()) / gime->MaxEvent();
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"nTrackSegments vaut "<<nTrackSegments<<endl;
  cerr<<"__________________________________________________________________"<<endl; 
   if ( nTrackSegments < maxTrackSegments-0.15 || nTrackSegments > maxTrackSegments+0.15 ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"       MESS ==> Error detected in the TrackSegments process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/galice.root galice.root | mail -s "PHOS INSTALLATION ERROR" schutz@in2p3.fr");
  }

 
}
