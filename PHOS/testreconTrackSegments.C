#include " AliPHOSTrackSegmentMaker.h"
#include " AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSGetter.h"
#include "TSystem.h"

void testreconTrackSegments(Int_t nevent = 1, const char *config="testconfig.C")

{ 
  const Float_t maxTrackSegments = 1 ;
  const Float_t widTrackSegments = TMath::Sqrt(maxTrackSegments) ;
  TString name = "test suite" ;
 
  AliPHOSTrackSegmentMaker * tracks = new AliPHOSTrackSegmentMakerv1("testPHOS.root",name.Data());
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  tracks->ExecuteTask("deb");
  Float_t nTrackSegments =  (Float_t) (gime->TrackSegmentMaker(name.Data())->GetTrackSegmentsInRun()) / gime->MaxEvent();
 
   if ( nTrackSegments < maxTrackSegments-0.25 || nTrackSegments > maxTrackSegments+0.25 ) {
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the TrackSegments process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
   // gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' schutz@in2p3.fr");
 }
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> TrackSegments process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;
}

