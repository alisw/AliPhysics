/* $Id$ */

//
// script to run the AliMultiplicityESDSelector
//

#include "../CreateESDChain.C"

void runMultiplicitySelector(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE)
{
  TStopwatch timer;
  timer.Start();

  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libESD");
  gSystem->Load("libPWG0base");
  if (aMC != kFALSE)
    gSystem->Load("libPWG0dep");

  gROOT->ProcessLine(".L CreatedNdEta.C");
  gROOT->ProcessLine(".L CreateCuts.C");
  gROOT->ProcessLine(".L drawPlots.C");

  TChain* chain = CreateESDChain(data, nRuns, offset);
  TVirtualProof* proof = 0;

  if (aProof != kFALSE)
  {
    proof = TProof::Open("jgrosseo@lxb6046");

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

  TString selectorName = ((aMC == kFALSE) ? "AliMultiplicityESDSelector" : "AliMultiplicityMCSelector");
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx++";

  if (aDebug != kFALSE)
    selectorName += "g";

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
}

