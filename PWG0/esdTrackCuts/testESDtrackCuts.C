/* $Id$ */

//
// script to run the AliMultiplicityESDSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

TChain* testESDtrackCuts(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, Char_t* proofServer = "lxb6046")
{
  if (aProof)
    connectProof(proofServer);

  TString libraries("libEG;libGeom;libESD;libPWG0base");
  TString packages("PWG0base");

  if (!prepareQuery(libraries, packages, 1))
    return;

  // selection of esd tracks
  AliESDtrackCuts* esdTrackCutsAll = new AliESDtrackCuts("esdTrackCutsAll");

  esdTrackCutsAll->DefineHistograms(1);
  esdTrackCutsAll->SetMinNClustersTPC(50);
  esdTrackCutsAll->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsAll->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsAll->SetRequireTPCRefit(kTRUE);
  esdTrackCutsAll->SetMinNsigmaToVertex(3);
  esdTrackCutsAll->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsAll->SetAcceptKingDaughters(kFALSE);

  TList inputList;
  inputList.Add(esdTrackCutsAll);

  TChain* chain = CreateESDChain(data, nRuns, offset);

  TString selectorName = "AliTestESDtrackCutsSelector";
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx++";

  if (aDebug != kFALSE)
    selectorName += "g";

  Int_t result = executeQuery(chain, &inputList, selectorName);

  if (result != 0)
  {
    printf("ERROR: Executing process failed with %d.\n", result);
    return;
  }
}

void draw(const char* dir, const char* fileName = "trackCuts.root")
{
  /*
   draw("esdTrackCutsAll")
   draw("fEsdTrackCutsPri")
   draw("fEsdTrackCutsSec")
   draw("fEsdTrackCutsPlusZ")
   draw("fEsdTrackCutsMinusZ")
   draw("fEsdTrackCutsPos")
   draw("fEsdTrackCutsNeg")
  */

  gSystem->Load("libPWG0base");

  TFile::Open(fileName);

  AliESDtrackCuts* cuts = new AliESDtrackCuts(dir, dir);
  cuts->LoadHistograms();

  cuts->DrawHistograms();
}
