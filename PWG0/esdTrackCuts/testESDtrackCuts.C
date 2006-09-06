/* $Id$ */

//
// script to run the AliMultiplicityESDSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

TChain* testESDtrackCuts(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE)
{
  if (aProof)
    connectProof("proof01@lxb6046");

  TString libraries("libEG;libGeom;libESD;libPWG0base;libVMC;libMinuit;libSTEER;libPWG0dep;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;lib
EGPythia6;libAliPythia6");
  TString packages("PWG0base;PWG0dep");

  if (!prepareQuery(libraries, packages, kTRUE))
    return;

  // selection of esd tracks
  AliESDtrackCuts* esdTrackCutsAll = new AliESDtrackCuts("esdTrackCutsAll");
  AliESDtrackCuts* esdTrackCutsPri = new AliESDtrackCuts("esdTrackCutsPri");
  AliESDtrackCuts* esdTrackCutsSec = new AliESDtrackCuts("esdTrackCutsSec");

  esdTrackCutsAll->DefineHistograms(1);
  esdTrackCutsAll->SetMinNClustersTPC(50);
  esdTrackCutsAll->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsAll->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsAll->SetRequireTPCRefit(kTRUE);
  esdTrackCutsAll->SetMinNsigmaToVertex(3);
  esdTrackCutsAll->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsAll->SetAcceptKingDaughters(kFALSE);

  esdTrackCutsPri->DefineHistograms(4);
  esdTrackCutsPri->SetMinNClustersTPC(50);
  esdTrackCutsPri->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsPri->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsPri->SetRequireTPCRefit(kTRUE);
  esdTrackCutsPri->SetMinNsigmaToVertex(3);
  esdTrackCutsPri->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsPri->SetAcceptKingDaughters(kFALSE);

  esdTrackCutsSec->DefineHistograms(2);
  esdTrackCutsSec->SetMinNClustersTPC(50);
  esdTrackCutsSec->SetMaxChi2PerClusterTPC(3.5);
  esdTrackCutsSec->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  esdTrackCutsSec->SetRequireTPCRefit(kTRUE);
  esdTrackCutsSec->SetMinNsigmaToVertex(3);
  esdTrackCutsSec->SetRequireSigmaToVertex(kTRUE);
  esdTrackCutsSec->SetAcceptKingDaughters(kFALSE);


  TList inputList;
  inputList.Add(esdTrackCutsAll);
  inputList.Add(esdTrackCutsPri);
  inputList.Add(esdTrackCutsSec);

  TChain* chain = CreateESDChain(data, nRuns, offset);

  TString selectorName = "AliTestESDtrackCutsSelector";
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx+";

  if (aDebug != kFALSE)
    selectorName += "g";

  Int_t result = executeQuery(chain, &inputList, selectorName);

  if (result != 0)
  {
    printf("ERROR: Executing process failed with %d.\n", result);
    return;
  }
}

