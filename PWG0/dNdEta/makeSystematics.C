/* $Id$ */

//
// Script to run the creation of input for systematics
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void makeSystematics(Char_t* dataDir, Int_t nRuns=20, Int_t offset = 0, Bool_t debug = kFALSE, Bool_t aProof = kFALSE, const Char_t* option = "")
{
  if (aProof)
    connectProof("proof01@lxb6046");

  TString libraries("libEG;libGeom;libESD;libVMC;libMinuit;libSTEER;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;libEGPythia6;libAliPythia6;libPWG0base;libPWG0dep");
  TString packages("PWG0base;PWG0dep");

  if (!prepareQuery(libraries, packages, kTRUE))
    return;

  gROOT->ProcessLine(".L CreateCuts.C");

  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  TList inputList;
  inputList.Add(esdTrackCuts);

  TChain* chain = CreateESDChain(dataDir, nRuns, offset);

  TString selector("AlidNdEtaSystematicsSelector.cxx++");
  if (debug != kFALSE)
    selector += "g";

  Int_t result = executeQuery(chain, &inputList, selector, option);
}
