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

void rundNdEtaAnalysis(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aMC = kFALSE, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction", const char* option = "", const char* proofServer = "lxb6046")
{
  if (aProof)
    connectProof(proofServer);

  TString libraries("libEG;libGeom;libESD;libPWG0base");
  TString packages("PWG0base");

  if (!prepareQuery(libraries, packages, kTRUE))
    return;

  gROOT->ProcessLine(".L CreateCuts.C");
  gROOT->ProcessLine(".L drawPlots.C");

  TChain* chain = CreateESDChain(data, nRuns, offset);

  TList inputList;

  // selection of esd tracks
  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  inputList.Add(esdTrackCuts);

  TString selectorName = ((aMC == kFALSE) ? "AlidNdEtaAnalysisESDSelector" : "AlidNdEtaAnalysisMCSelector");
  AliLog::SetClassDebugLevel(selectorName, AliLog::kInfo);

  selectorName += ".cxx+";

  if (aDebug != kFALSE)
    selectorName += "+g";

  Int_t result = executeQuery(chain, &inputList, selectorName, option);

  if (result >= 0)
  {
    if (aMC)
    {
      dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");

      TFile* file = TFile::Open("analysis_mc.root");

      if (!file)
      {
        cout << "Error. File not found" << endl;
        return;
      }
      fdNdEtaAnalysis->LoadHistograms();
      fdNdEtaAnalysis->DrawHistograms(kTRUE);
    }
    else
      FinishAnalysisAll("analysis_esd_raw.root", "analysis_esd.root", correctionMapFile, correctionMapFolder);
  }
}

