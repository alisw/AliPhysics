
#include "ROOT/RDFHistoModels.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TFile.h"

#include "AliAnalysisTaskAO2Dconverter.h"

#include "AliAnalysisTaskAO2Dconverter.h"

void read()
{
  const Char_t* name = AliAnalysisTaskAO2Dconverter::TreeName[AliAnalysisTaskAO2Dconverter::kTracks];
  ROOT::RDataFrame d(name, "AO2D.root");
  ROOT::RDF::TH2DModel betamodel("beta", ";#it{p} (GeV/#it{c});TOF #beta;", 1000, 0.1, 5, 1000, 0, 2);
  auto h = d.Define("beta", "fTOFsignal > 0.1 ? fLength/(fTOFsignal * (TMath::C() * 1.e2 / 1.e12)) : -10").Histo2D(betamodel, "fTPCinnerP", "beta");
  new TCanvas();
  gPad->SetLogz();
  gPad->SetLogx();
  auto hdrawn = h->DrawCopy("COLZ");
  hdrawn->SetDirectory(0);
}
