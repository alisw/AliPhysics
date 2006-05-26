/* $Id$ */

//
// Script to make correction maps for dndeta measurements using the
// dNdEtaCorrection class.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"

void makeCorrection2(Char_t* dataDir, Int_t nRuns=20, Int_t offset = 0)
{
  gSystem->Load("libPWG0base");

  TChain* chain = CreateESDChainFromDir(dataDir, nRuns, offset, kFALSE);

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
  chain->Process("AlidNdEtaCorrectionSelector.cxx+");
}
