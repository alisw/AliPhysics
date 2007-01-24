/* $Id$ */

//
// script to run the AliMultiplicityESDSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void runMultiplicitySelector(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE)
{
  if (aProof)
    connectProof("proof01@lxb6046");

  TString libraries("libEG;libGeom;libESD;libPWG0base");
  TString packages("PWG0base");

  if (aMC != kFALSE)
  {
    libraries += ";libVMC;libMinuit;libSTEER;libPWG0dep;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;libEGPythia6;libAliPythia6";
    packages += ";PWG0dep";
  }

  if (!prepareQuery(libraries, packages, kTRUE))
    return;

  gROOT->ProcessLine(".L CreateCuts.C");
  gROOT->ProcessLine(".L drawPlots.C");

  // selection of esd tracks
  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  TList inputList;
  inputList.Add(esdTrackCuts);

  TChain* chain = CreateESDChain(data, nRuns, offset);

  TString selectorName = ((aMC == kFALSE) ? "AliMultiplicityESDSelector" : "AliMultiplicityMCSelector");
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

void draw(const char* fileName = "multiplicityMC.root")
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  mult->DrawHistograms();
}


void* fit(const char* fileName = "multiplicityMC.root")
{
  gSystem->Load("libPWG0base");

  AliMultiplicityCorrection* mult = new AliMultiplicityCorrection("Multiplicity", "Multiplicity");

  TFile::Open(fileName);
  mult->LoadHistograms("Multiplicity");

  mult->ApplyMinuitFit(3, kFALSE);

  return mult;
}
