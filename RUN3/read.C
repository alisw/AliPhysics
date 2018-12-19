
#include "ROOT/TDataFrame.hxx"
#include "TCanvas.h"
#include "TFile.h"

#include "AliAnalysisTaskAO2Dconverter.h"

using namespace ROOT::Experimental; // TDataFrame's namespace

void read()
{
  const Char_t* name = AliAnalysisTaskAO2Dconverter::TreeName[AliAnalysisTaskAO2Dconverter::kTracks];
  TDataFrame d(name, "AO2D.root");
  auto h = d.Define("beta", "fTOFsignal > 0.1 ? fLength/(fTOFsignal * (TMath::C() * 1.e2 / 1.e12)) : -10").Histo2D(::TH2F{ "beta", ";#it{p} (GeV/#it{c});TOF #beta;", 1000, 0.1, 5, 1000, 0, 2 }, "fTPCinnerP", "beta");
  new TCanvas();
  gPad->SetLogz();
  gPad->SetLogx();
  auto hdrawn = h->DrawCopy("COLZ");
  hdrawn->SetDirectory(0);
}
