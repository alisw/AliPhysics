/* $Id$ */

//
// Script to make correction maps for dndeta measurements using the
// dNdEtaCorrection class.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"

void makeCorrection2(Char_t* dataDir, Int_t nRuns=20, Int_t offset = 0, Bool_t debug = kFALSE, const Char_t* option = "")
{
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");

  TChain* chain = CreateESDChain(dataDir, nRuns, offset);

  fEsdTrackCuts = new AliESDtrackCuts();
  fEsdTrackCuts->DefineHistograms(1);

  fEsdTrackCuts->SetMinNClustersTPC(50);
  fEsdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  fEsdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  fEsdTrackCuts->SetRequireTPCRefit(kTRUE);

  fEsdTrackCuts->SetMinNsigmaToVertex(3);
  fEsdTrackCuts->SetAcceptKingDaughters(kFALSE);

  chain->GetUserInfo()->Add(fEsdTrackCuts);

  AliLog::SetClassDebugLevel("AlidNdEtaCorrectionSelector", AliLog::kInfo);
  AliLog::SetClassDebugLevel("AliSelectorRL", AliLog::kInfo);

  TString selector("AlidNdEtaCorrectionSelector.cxx+");
  if (debug != kFALSE)
    selector += "g";

  chain->Process(selector, option);
}
