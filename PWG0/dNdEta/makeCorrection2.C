/* $Id$ */

//
// Script to make correction maps for dndeta measurements using the
// dNdEtaCorrection class.
//
// implementation with TSelector
//

#include "../CreateESDChain.C"
#include "../PWG0Helper.C"

void makeCorrection2(Char_t* dataDir, Int_t nRuns=20, Int_t offset = 0, Bool_t debug = kFALSE, Bool_t aProof = kFALSE, const Char_t* option = "", const Char_t* proofConnect = "proof01@lxb6046")
{
  if (aProof)
    connectProof(proofConnect);

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

  TString selector("AlidNdEtaCorrectionSelector.cxx++");
  if (debug != kFALSE)
    selector += "g";

  Int_t result = executeQuery(chain, &inputList, selector, option);
}

void VerifyCorrection(Int_t draw = 0, const char* correctionMapFile = "correction_map.root", const char* correctionMapFolder = "dndeta_correction", const char* analysisDir = "dndetaESD", const char* checkDir = "dndetaMC")
{
  gSystem->Load("libPWG0base");

  AlidNdEtaCorrection* dNdEtaCorrection = new AlidNdEtaCorrection(correctionMapFolder, correctionMapFolder);
  dNdEtaCorrection->LoadHistograms(correctionMapFile, correctionMapFolder);
  dNdEtaCorrection->SetName("dndeta_correction");
  dNdEtaCorrection->SetTitle("dndeta_correction");

  TFile::Open(correctionMapFile);

  dNdEtaAnalysis* fdNdEtaAnalysis = new dNdEtaAnalysis(analysisDir, analysisDir);
  fdNdEtaAnalysis->LoadHistograms();

  // correct with track2particle
  TH3F* hist = fdNdEtaAnalysis->GetHistogram();
  for (Int_t x=1; x<=hist->GetNbinsX(); ++x)
    for (Int_t y=1; y<=hist->GetNbinsY(); ++y)
      for (Int_t z=1; z<=hist->GetNbinsZ(); ++z)
      {
        Float_t correction = dNdEtaCorrection->GetTrack2ParticleCorrection(hist->GetXaxis()->GetBinCenter(x), hist->GetYaxis()->GetBinCenter(y), hist->GetZaxis()->GetBinCenter(z));
        Float_t value = hist->GetBinContent(x, y, z);

        hist->SetBinContent(x, y, z, correction * value);
      }

  fdNdEtaAnalysis->Finish(dNdEtaCorrection, 0.3);

  dNdEtaAnalysis* fdNdEtaAnalysisMC = new dNdEtaAnalysis(checkDir, checkDir);
  fdNdEtaAnalysisMC->LoadHistograms();

  fdNdEtaAnalysisMC->Finish(0, 0.3);

  if (draw == 0)
    fdNdEtaAnalysis->DrawHistograms();
  else if (draw == 1)
    fdNdEtaAnalysisMC->DrawHistograms();
  else
  {
    for (Int_t i=0; i<dNdEtaAnalysis::kVertexBinning; ++i)
      fdNdEtaAnalysis->GetdNdEtaHistogram(i)->Divide(fdNdEtaAnalysisMC->GetdNdEtaHistogram(i));

    fdNdEtaAnalysis->DrawHistograms();
  }
}
