#include " AliPHOSTrackSegmentMaker.h"
#include "AliPHOSGetter.h"

void testreconTrackSegments(Int_t nevent = 1, const char *config="testconfig.C")

{
 AliPHOSGetter * gime = AliPHOSGetter::GetInstance("galice.root") ;
 AliPHOSTrackSegmentMaker * tracks = new AliPHOSTrackSegmentMakerv1("galice.root","test suite");
 tracks->ExecuteTask("deb");
 //cout << "# of trackSegments " << gime->trackSegments()<<GetEntries()<< endl;

}
