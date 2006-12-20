/* $Id$ */

//
// script to run the AliVertexSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void runVertexSelector(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aProof = kFALSE)
{
  if (aProof)
    connectProof("jgrosseo@lxb6046");

  TString libraries("libEG;libGeom;libESD;libPWG0base");
  TString packages("PWG0base");

  //libraries += ";libVMC;libMinuit;libSTEER;libPWG0dep;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;libEGPythia6;libAliPythia6";

  if (!prepareQuery(libraries, packages, 1))
    return;

  TChain* chain = CreateESDChain(data, nRuns, offset);

  TList inputList;

  gROOT->ProcessLine(".L CreateCuts.C");

  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  inputList.Add(esdTrackCuts);

  executeQuery(chain, &inputList, "AliVertexSelector.cxx+");
}

