/* $Id$ */

//
// Script to run the creation of input for systematics
//

#include "../CreateESDChain.C"

void makeSystematics(Char_t* dataDir, Int_t nRuns=20, Int_t offset = 0, Bool_t debug = kFALSE, const Char_t* option = "")
{
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");

  gROOT->ProcessLine(".L CreateCuts.C");

  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  TChain* chain = CreateESDChain(dataDir, nRuns, offset);
  chain->GetUserInfo()->Add(esdTrackCuts);

  TString selector("AlidNdEtaSystematicsSelector.cxx+");
  if (debug != kFALSE)
    selector += "g";

  chain->Process(selector, option);
}
