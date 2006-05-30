/* $Id$ */

//
// Script to test the dN/dEta analysis using the dNdEtaAnalysis and
// dNdEtaCorrection classes. Note that there is a cut on the events,
// so the measurement will be biassed.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"

testAnalysis2(Char_t* dataDir, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE)
{
  gSystem->Load("libPWG0base");

  TChain* chain = CreateESDChainFromDir(dataDir, nRuns, offset);

  // ########################################################
  // selection of esd tracks
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
  esdTrackCuts->DefineHistograms(1);

  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);

  esdTrackCuts->SetMinNsigmaToVertex(3);
  esdTrackCuts->SetAcceptKingDaughters(kFALSE);

  chain->GetUserInfo()->Add(esdTrackCuts);

  if (aMC == kFALSE)
  {
    dNdEtaCorrection* dNdEtaCorrection = new dNdEtaCorrection();
    dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");
    dNdEtaCorrection->RemoveEdges(2, 0, 2);

    chain->GetUserInfo()->Add(dNdEtaCorrection);
  }

  TString selectorName = ((aMC == kFALSE) ? "AlidNdEtaAnalysisESDSelector" : "AlidNdEtaAnalysisMCSelector");
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx++";
  if (aDebug != kFALSE)
    selectorName += "g";

  TStopwatch timer;
  timer.Start();

  chain->Process(selectorName);

  timer.Stop();
  timer.Print();
}
