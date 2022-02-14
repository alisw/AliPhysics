#include "/Users/stephanstiefelmaier/work/projects/a_helper/usefullFunctions.C"

#if !defined(__CINT__) || defined(__CLING__)
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGGA/GammaConv/macros/AddTask_V0Reader.C>
#include <PWGGA/GammaConv/macros/AddTask_PhotonQA.C>
#include <PWGGA/GammaConv/macros/AddTask_GammaConvV1_PbPb.C>

#endif

void PrintPhotons(AliConversionPhotonCuts::TMapPhotonBool &theMap){
  cout << "printPhotons start\nGetTrackLabelPositive() GetTrackLabelNegative() GetChi2perNDF()\n";

  for (auto it = theMap.begin(); it!= theMap.end(); ++it){
    cout << it->first->GetTrackLabelPositive() << " " << it->first->GetTrackLabelNegative() << " " << it->first->GetChi2perNDF() << endl;
  }
 cout << "end\n";
}

Bool_t TestSharedElectronCut(AliConversionPhotonCuts &thePhotonCuts){

  // define photon sample        x       x    x       x
  vector<Int_t> posLabels     = {0,  2,  5,   2,  8,  7,  5 };
  vector<Int_t> negLabels     = {1,  3,  1,   4,  7,  9, 12 };
  vector<Double_t> chi2s      = {5., 4., 3.,  4., 8., 9., 1.};
  vector<Bool_t> expectations = {1,  0,  1,   1,  0,  1,  0 }; // 1 means should get kicked if working as expected
  Size_t nPhotons = posLabels.size();

  // create photon sample
  AliConversionPhotonCuts::TMapPhotonBool lSample_before;
  AliConversionPhotonCuts::TMapPhotonBool lSample_targetAfter;

  // lets contruct the AliAODConversionPhotons into a TClonesArray. This should (I hope) ensure
  // a fixed order in memory of the photons each time this code runs
  TClonesArray lArray("AliAODConversionPhoton", 100);
  for (auto i : RangeFrom0(nPhotons)){

    new(lArray[i]) AliAODConversionPhoton();
    AliAODConversionPhoton *lPhot = dynamic_cast<AliAODConversionPhoton*>(lArray[i]);
    lPhot->SetTrackLabels(posLabels[i], negLabels[i]);
    lPhot->SetChi2perNDF(chi2s[i]);

    lSample_before.insert({lPhot, kTRUE});
    if (!expectations[i]){ lSample_targetAfter.insert({lPhot, kTRUE}); }
  }

  // check if indeed the order got preserved
  Size_t i = 0;
  for (const auto &iPhotBool : lSample_before){
    if (iPhotBool.first != lArray[i]) { break; }
    ++i;
  }
  if (i!=nPhotons) {
    cout << "INFO: TestSharedElectronCut(): The order did not get preserved. If this test fails, this might be the reason.\n";
  }

  // apply cut
  thePhotonCuts.RemovePhotonsWithSharedTracks(lSample_before);
  lArray.Delete();

  // check if the right photons were removed
  return lSample_before==lSample_targetAfter;
}


void runTests_GammaConvROOT6(){

  AliConversionPhotonCuts lPhotCuts;
  cout << (TestSharedElectronCut(lPhotCuts) ? "TestSharedElectronCut OK." : "TestSharedElectronCut failed!" )  << endl;
  return;
}
