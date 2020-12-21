
#include "ROOT/RDFHistoModels.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TFile.h"

#include "AliAnalysisTaskAO2Dconverter.h"

#include "AliExternalTrackParam.h"

void read(Bool_t isMC = kTRUE, const Char_t* fname = "AO2D.root")
{
  ROOT::RDataFrame dtrk(AliAnalysisTaskAO2Dconverter::TreeName[AliAnalysisTaskAO2Dconverter::kTracks].Data(), fname);
  ROOT::RDF::TH2DModel betamodel("beta", ";#it{p} (GeV/#it{c});TOF #beta;", 1000, 0.1, 5, 1000, 0, 2);
  auto h = dtrk.Define("param", "std::array<double, 5> p{fY, fZ, fSnp, fTgl, fSigned1Pt}; return p;")
               .Define("covar", "std::array<double, 15> p{fCYY, fCZY, fCZZ, fCSnpY, fCSnpZ, fCSnpSnp, fCTglY, fCTglZ, fCTglSnp, fCTglTgl, fC1PtY, fC1PtZ, fC1PtSnp, fC1PtTgl, fC1Pt21Pt2}; return p;")
               .Define("P", "AliExternalTrackParam((double)fX, (double)fAlpha, param.data(), covar.data()).GetP()")
               .Define("beta", "fTOFsignal > 0.1 ? fLength/(fTOFsignal * (TMath::C() * 1.e2 / 1.e12)) : -10")
               .Histo2D(betamodel, "P", "beta");
  new TCanvas();
  gPad->SetLogz();
  gPad->SetLogx();
  auto hdrawn = h->DrawCopy("COLZ");
  hdrawn->SetDirectory(0);
  if (isMC) {
    ROOT::RDataFrame dev(AliAnalysisTaskAO2Dconverter::TreeName[AliAnalysisTaskAO2Dconverter::kEvents].Data(), fname);
    auto hij = [](Short_t id) { Printf("Is %sHijiing", TESTBIT(id, AliAnalysisTaskAO2Dconverter::kAliGenCocktailEventHeader) ? "" : "not "); };
    dev.Foreach(hij, { "fGeneratorID" });
  }
}
