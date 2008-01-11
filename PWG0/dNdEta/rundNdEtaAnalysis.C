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

void loadlibs()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libPWG0base");
}

void FinishAnalysisAll(const char* dataInput = "analysis_esd_raw.root", const char* dataOutput = "analysis_esd.root", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction")
{
  loadlibs();

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  TFile::Open(correctionMapFile);
  dNdEtaCorrection->LoadHistograms();

  TFile* file = TFile::Open(dataInput);

  if (!file)
  {
    cout << "Error. File not found" << endl;
    return;
  }

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis("dndeta", "dndeta");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kINEL);
  fdNdEtaAnalysis->DrawHistograms(kTRUE);
  TFile* file2 = TFile::Open(dataOutput, "RECREATE");
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTr", "dndetaTr");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kVertexReco);
  fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTrVtx", "dndetaTrVtx");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3, AlidNdEtaCorrection::kTrack2Particle);
  fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();

  file->cd();
  fdNdEtaAnalysis = new dNdEtaAnalysis("dndetaTracks", "dndetaTracks");
  fdNdEtaAnalysis->LoadHistograms("fdNdEtaAnalysisESD");
  fdNdEtaAnalysis->Finish(0, 0.3, AlidNdEtaCorrection::kNone);
  fdNdEtaAnalysis->DrawHistograms(kTRUE);
  file2->cd();
  fdNdEtaAnalysis->SaveHistograms();
}

void* FinishAnalysis(const char* analysisFile = "analysis_esd.root", const char* analysisDir = "dndeta", const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction", Bool_t useUncorrected = kFALSE, Bool_t simple = kFALSE)
{
  loadlibs();

  TFile* file = TFile::Open(analysisFile);

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(analysisDir, analysisDir);
  fdNdEtaAnalysis->LoadHistograms();

  if (correctionMapFile)
  {
    AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
    TFile::Open(correctionMapFile);
    dNdEtaCorrection->LoadHistograms();

    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.2, AlidNdEtaCorrection::kINEL);
    fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kINEL);
    //fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0, AlidNdEtaCorrection::kTrack2Particle);
  }
  else
    fdNdEtaAnalysis->Finish(0, 0, AlidNdEtaCorrection::kNone);

  fdNdEtaAnalysis->DrawHistograms(simple);

  TH1* hist = fdNdEtaAnalysis->GetdNdEtaHistogram(1);
  Int_t binLeft = hist->GetXaxis()->FindBin(-0.5);
  Int_t binRight = hist->GetXaxis()->FindBin(0.5);
  Float_t value1 = hist->Integral(binLeft, binRight);

  hist = fdNdEtaAnalysis->GetdNdEtaHistogram(2);
  Float_t value2 = hist->Integral(binLeft, binRight);

  if (value2 > 0)
    printf("Ratio is %f, values are %f %f\n", value1 / value2, value1, value2);

  return fdNdEtaAnalysis;
}
