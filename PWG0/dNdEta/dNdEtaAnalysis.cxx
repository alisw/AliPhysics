#include "dNdEtaAnalysis.h"

#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>

//____________________________________________________________________
ClassImp(dNdEtaAnalysis)

//____________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(Char_t* name) {

  fName = TString(name);

  hEtaVsVtx  = new TH2F("eta_vs_vtx","",80,-20,20,120,-6,6);
  hEtaVsVtx->SetXTitle("vtx z [cm]");
  hEtaVsVtx->SetYTitle("#eta");

  hEtaVsVtxUncorrected = dynamic_cast<TH2F*> (hEtaVsVtx->Clone("eta_vs_vtx_uncorrected"));
  hVtx       = hEtaVsVtx->ProjectionX("vtx");
  hdNdEta    = hEtaVsVtx->ProjectionY("dNdEta");

  hdNdEta->SetYTitle("dN/d#eta");

  hEtaVsVtx->Sumw2();
  hVtx->Sumw2();
}

//____________________________________________________________________
void
dNdEtaAnalysis::FillTrack(Float_t vtx, Float_t eta, Float_t weight) {
  hEtaVsVtx->Fill(vtx, eta, weight);
  hEtaVsVtxUncorrected->Fill(vtx,eta);
}

//____________________________________________________________________
void
dNdEtaAnalysis::FillEvent(Float_t vtx) {
  hVtx->Fill(vtx);
}

//____________________________________________________________________
void
dNdEtaAnalysis::Finish() {  

  // first normalize with n events (per vtx)
  for (Int_t i_vtx=0; i_vtx<=hVtx->GetNbinsX(); i_vtx++) {
    Float_t nEvents      = hVtx->GetBinContent(i_vtx);
    Float_t nEventsError = hVtx->GetBinError(i_vtx);
    
    if (nEvents==0) continue;
    
    for (Int_t i_eta=0; i_eta<=hEtaVsVtx->GetNbinsY(); i_eta++) {
      Float_t value = hEtaVsVtx->GetBinContent(i_vtx, i_eta)/nEvents;
      if (value==0) continue;
      Float_t error = hEtaVsVtx->GetBinError(i_vtx, i_eta)/nEvents;
      error = TMath::Sqrt(TMath::Power(hEtaVsVtx->GetBinError(i_vtx, i_eta)/
				       hEtaVsVtx->GetBinContent(i_vtx, i_eta),2) +
      			  TMath::Power(nEventsError/nEvents,2));
      hEtaVsVtx->SetBinContent(i_vtx, i_eta, value);
      hEtaVsVtx->SetBinError(i_vtx,   i_eta, error);
    }
  }

  // then take the wieghted average for each eta
  // is this the right way to do it???
  for (Int_t i_eta=0; i_eta<=hEtaVsVtx->GetNbinsY(); i_eta++) {
    Float_t sum           = 0;
    Float_t sumError2     = 0;
    Int_t   nMeasurements = 0;    

    Float_t sum_xw = 0;
    Float_t sum_w  = 0;
    
    for (Int_t i_vtx=0; i_vtx<=hVtx->GetNbinsX(); i_vtx++) {
      if (hVtx->GetBinContent(i_vtx)==0)             continue;
      if (hEtaVsVtx->GetBinContent(i_vtx, i_eta)==0) continue;

      Float_t w = 1/TMath::Power(hEtaVsVtx->GetBinError(i_vtx, i_eta),2);
      sum_xw = sum_xw + hEtaVsVtx->GetBinContent(i_vtx, i_eta)*w;
      sum_w  = sum_w + w;

      sum = sum + hEtaVsVtx->GetBinContent(i_vtx, i_eta);
      sumError2 = sumError2 + TMath::Power(hEtaVsVtx->GetBinError(i_vtx, i_eta),2);      
      nMeasurements++;
    }
    Float_t dndeta = 0;
    Float_t error  = 0;

    if (nMeasurements!=0) {
      dndeta = sum/Float_t(nMeasurements);
      error  = TMath::Sqrt(sumError2)/Float_t(nMeasurements);
      
      dndeta = sum_xw/sum_w;
      error  = 1/TMath::Sqrt(sum_w);

      dndeta = dndeta/hdNdEta->GetBinWidth(i_eta);
      error  = error/hdNdEta->GetBinWidth(i_eta);

      hdNdEta->SetBinContent(i_eta, dndeta);
      hdNdEta->SetBinError(i_eta, error);
    }

  }
}



//____________________________________________________________________
void
dNdEtaAnalysis::SaveHistograms() {

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());
  
  hEtaVsVtx  ->Write();
  hVtx       ->Write();
  hdNdEta    ->Write();

  gDirectory->cd("../");
}

//____________________________________________________________________
void dNdEtaAnalysis::DrawHistograms()
{
  TCanvas* canvas = new TCanvas("dNdEtaAnalysis", "dNdEtaAnalysis", 800, 800);
  canvas->Divide(2, 2);

  canvas->cd(1);
  if (hEtaVsVtx)
    hEtaVsVtx->Draw("COLZ");

  canvas->cd(2);
  if (hEtaVsVtxUncorrected)
    hEtaVsVtxUncorrected->Draw("COLZ");

  canvas->cd(3);
  if (hVtx)
    hVtx->Draw();

  canvas->cd(4);
  if (hdNdEta)
    hdNdEta->Draw();
}
