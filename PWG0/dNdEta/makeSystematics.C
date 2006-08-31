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

  TString libraries("libEG;libGeom;libESD;libPWG0base;libVMC;libMinuit;libSTEER;libPWG0dep;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;libEGPythia6;libAliPythia6");
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

void runAnalysisWithDifferentMaps(Char_t* dataDir, Int_t nRuns=20, Int_t offset = 0, Bool_t debug = kFALSE, Bool_t proof = kFALSE)
{
  // This functions runs the dN/dEta analysis with different correction maps to gather systematics
  // It runs with the "normal map", and with 4 other different cases where particle species are enhanced
  // or reduced.
  // The normal map is expected in correction_map.root, created by AlidNdEtaCorrectionSelector
  // The others in new_compositions.root in the folders (K|p)(Boosted|Reduced), created
  //   by AlidNdEtaSystematicsSelector and Composition() out of drawSystematics.C

  gROOT->ProcessLine(".L testAnalysis2.C");

  gSystem->Exec("rm analysis_esd.root");

  testAnalysis2(dataDir, nRuns, offset, kFALSE, debug, proof, "correction_map.root", "dndeta_correction");
  if (gSystem->Exec("mv analysis_esd.root systematics_dndeta_reference.root") != 0)
  {
    printf("systematics_dndeta_reference.root failed\n");
    return;
  }

  testAnalysis2(dataDir, nRuns, offset, kFALSE, debug, proof, "new_compositions.root", "KBoosted");
  if (gSystem->Exec("mv analysis_esd.root systematics_dndeta_KBoosted.root") != 0)
  {
    printf("systematics_dndeta_KBoosted.root failed\n");
    return;
  }

  testAnalysis2(dataDir, nRuns, offset, kFALSE, debug, proof, "new_compositions.root", "KReduced");
  if (gSystem->Exec("mv analysis_esd.root systematics_dndeta_KReduced.root") != 0)
  {
    printf("systematics_dndeta_KReduced.root failed\n");
    return;
  }

  testAnalysis2(dataDir, nRuns, offset, kFALSE, debug, proof, "new_compositions.root", "pBoosted");
  if (gSystem->Exec("mv analysis_esd.root systematics_dndeta_pBoosted.root") != 0)
  {
    printf("systematics_dndeta_pBoosted.root failed\n");
    return;
  }

  testAnalysis2(dataDir, nRuns, offset, kFALSE, debug, proof, "new_compositions.root", "pReduced");
  if (gSystem->Exec("mv analysis_esd.root systematics_dndeta_pReduced.root") != 0)
  {
    printf("systematics_dndeta_pReduced.root failed\n");
    return;
  }
}
