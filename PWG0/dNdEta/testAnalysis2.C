/* $Id$ */

//
// Script to test the dN/dEta analysis using the dNdEtaAnalysis and
// dNdEtaCorrection classes. Note that there is a cut on the events,
// so the measurement will be biassed.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void testAnalysis2(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction", const char* option = "", const char* proofServer = "jgrosseo@lxb6046")
{
  if (aProof)
    connectProof(proofServer);

  TString libraries("libEG;libGeom;libESD;libPWG0base");
  TString packages("PWG0base");
  if (aMC != kFALSE)
  {
    libraries += ";libVMC;libMinuit;libSTEER;libEVGEN;libFASTSIM;libmicrocern;libpdf;libpythia6;libEGPythia6;libAliPythia6;libPWG0dep";
    packages += ";PWG0dep";
  }

  if (!prepareQuery(libraries, packages, kTRUE))
    return;

  //TODO somehow prevent loading several times
  gROOT->ProcessLine(".L CreateCuts.C");
  gROOT->ProcessLine(".L drawPlots.C");

  TChain* chain = CreateESDChain(data, nRuns, offset);


  TList inputList;

  if (aMC == kFALSE)
  {
    // selection of esd tracks
    AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
    if (!esdTrackCuts)
    {
      printf("ERROR: esdTrackCuts could not be created\n");
      return;
    }

    inputList.Add(esdTrackCuts);

    AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
    dNdEtaCorrection->LoadHistograms(correctionMapFile, correctionMapFolder);
    dNdEtaCorrection->ReduceInformation();
    dNdEtaCorrection->SetName("dndeta_correction");
    dNdEtaCorrection->SetTitle("dndeta_correction");

    inputList.Add(dNdEtaCorrection);
  }

  TString selectorName = ((aMC == kFALSE) ? "AlidNdEtaAnalysisESDSelector" : "AlidNdEtaAnalysisMCSelector");
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx+";

  if (aDebug != kFALSE)
    selectorName += "g";

  Int_t result = executeQuery(chain, &inputList, selectorName, option);

  if (result >= 0)
  {
    dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");

    TFile* file = TFile::Open(aMC ? "analysis_mc.root" : "analysis_esd.root");
    if (!file)
    {
      cout << "Error. File not found" << endl;
      return;
    }
    fdNdEtaAnalysis->LoadHistograms();
    fdNdEtaAnalysis->DrawHistograms();

    dNdEta(kTRUE);
  }
}
