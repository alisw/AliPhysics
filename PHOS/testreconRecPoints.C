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
;
 AliPHOSClusterizer * cluster = new  AliPHOSClusterizerv1("galice.root","test suite");
 AliPHOSGetter * gime = AliPHOSGetter::GetInstance(); 
 cluster->ExecuteTask("deb");
 Float_t nRecPoints =  (Float_t) (gime->RecPoints()->GetRecPointsInRun()) / gime->MaxEvent();
}
