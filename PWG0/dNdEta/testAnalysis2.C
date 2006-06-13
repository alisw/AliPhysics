/* $Id$ */

//
// Script to test the dN/dEta analysis using the dNdEtaAnalysis and
// dNdEtaCorrection classes. Note that there is a cut on the events,
// so the measurement will be biassed.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"

void testAnalysis2(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE)
{
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libESD");
  gSystem->Load("libPWG0base");
  if (aMC != kFALSE)
    gSystem->Load("libPWG0dep");

  gROOT->ProcessLine(".L CreatedNdEta.C");
  gROOT->ProcessLine(".L CreateCuts.C");

  TChain* chain = 0;
  TVirtualProof* proof = 0;
  if (aProof == kFALSE)
    chain = CreateESDChainFromDir(data, nRuns, offset);
  else
  {
    chain = CreateESDChainFromList(data, nRuns, offset);
    proof = gROOT->Proof("jgrosseo@lxb6046");

    if (!proof)
    {
      printf("ERROR: PROOF connection not established.\n");
      return;
    }

    if (proof->EnablePackage("ESD"))
    {
      printf("ERROR: ESD package could not be enabled.\n");
      return;
    }

    if (proof->EnablePackage("PWG0base"))
    {
      printf("ERROR: PWG0base package could not be enabled.\n");
      return;
    }

    if (aMC != kFALSE)
    {
      if (proof->EnablePackage("PWG0dep"))
      {
        printf("ERROR: PWG0dep package could not be enabled.\n");
        return;
      }
    }

    //chain->SetProof(proof);
  }

  // ########################################################
  // selection of esd tracks
  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  chain->GetUserInfo()->Add(esdTrackCuts);
  if (proof)
    proof->AddInput(esdTrackCuts);

  if (aMC == kFALSE)
  {
    AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection();
    dNdEtaCorrection->LoadHistograms("correction_map.root","dndeta_correction");
    //dNdEtaCorrection->RemoveEdges(2, 0, 2);

    chain->GetUserInfo()->Add(dNdEtaCorrection);
    if (proof)
      proof->AddInput(dNdEtaCorrection);
  }

  TString selectorName = ((aMC == kFALSE) ? "AlidNdEtaAnalysisESDSelector" : "AlidNdEtaAnalysisMCSelector");
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  // workaround for a bug in PROOF that only allows header files for .C files
  // please create symlink from <selector>.cxx to <selector>.C
  if (proof != kFALSE)
    selectorName += ".C+";
  else
    selectorName += ".cxx+";

  if (aDebug != kFALSE)
    selectorName += "g";

  TStopwatch timer;
  timer.Start();

  Long64_t result = -1;

  if (proof != kFALSE)
    result = chain->MakeTDSet()->Process(selectorName);
  else
    result = chain->Process(selectorName);

  if (result != 0)
  {
    printf("ERROR: Executing process failed with %d.\n", result);
    return;
  }

  timer.Stop();
  timer.Print();

  CreatedNdEta(aMC ? kFALSE : kTRUE, aMC ? "analysis_mc.root" : "analysis_esd.root");
}

