/* $Id$ */

#include "dNdEtaAnalysis.h"

#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCollection.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>

#include "dNdEtaCorrection.h"

//____________________________________________________________________
ClassImp(dNdEtaAnalysis)

//____________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(Char_t* name, Char_t* title) :
TNamed(name, title)
{
  hEtaVsVtx  = new TH2F(Form("%s_eta_vs_vtx", name),"",80,-20,20,120,-6,6);
  hEtaVsVtx->SetXTitle("vtx z [cm]");
  hEtaVsVtx->SetYTitle("#eta");

  hEtaVsVtxUncorrected = dynamic_cast<TH2F*> (hEtaVsVtx->Clone(Form("%s_eta_vs_vtx_uncorrected", name)));
  hEtaVsVtxCheck = dynamic_cast<TH2F*> (hEtaVsVtx->Clone(Form("%s_eta_vs_vtx_check", name)));
  hVtx       = hEtaVsVtx->ProjectionX(Form("%s_vtx", name));
  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    hdNdEta[i]    = hEtaVsVtx->ProjectionY(Form("%s_dNdEta_%d", name, i));
    hdNdEta[i]->SetYTitle("dN/d#eta");
  }

  hEtaVsVtx->Sumw2();
  hVtx->Sumw2();
}

//____________________________________________________________________
void
dNdEtaAnalysis::FillTrack(Float_t vtx, Float_t eta, Float_t c) {
  hEtaVsVtxUncorrected->Fill(vtx,eta);
  hEtaVsVtxCheck->Fill(vtx, eta, c);
}

//____________________________________________________________________
void
dNdEtaAnalysis::FillEvent(Float_t vtx) {
  hVtx->Fill(vtx);
}

//____________________________________________________________________
void dNdEtaAnalysis::Finish(dNdEtaCorrection* correction)
{
  // correct with correction values if available
  // TODO what do we do with the error?
  if (!correction)
    printf("INFO: No correction applied\n");

  // this can be replaced by TH2F::Divide if we agree that the binning will be always the same
  for (Int_t iVtx=0; iVtx<=hEtaVsVtxUncorrected->GetNbinsX(); iVtx++)
  {
    for (Int_t iEta=0; iEta<=hEtaVsVtxUncorrected->GetNbinsY(); iEta++)
    {
      Float_t correctionValue = 1;
      if (correction)
        correctionValue = correction->GetCorrection(hEtaVsVtxUncorrected->GetXaxis()->GetBinCenter(iVtx), hEtaVsVtxUncorrected->GetYaxis()->GetBinCenter(iEta));

      Float_t value = hEtaVsVtxUncorrected->GetBinContent(iVtx, iEta);
      Float_t error = hEtaVsVtxUncorrected->GetBinError(iVtx, iEta);

      Float_t correctedValue = value * correctionValue;
      Float_t correctedError = error * correctionValue;

      if (correctedValue != 0)
      {
        hEtaVsVtx->SetBinContent(iVtx, iEta, correctedValue);
        hEtaVsVtx->SetBinError(iVtx, iEta, correctedError);
      }
    }
  }

  // normalize with n events (per vtx)
  for (Int_t iVtx=0; iVtx<=hVtx->GetNbinsX(); iVtx++) {
    Float_t nEvents      = hVtx->GetBinContent(iVtx);
    Float_t nEventsError = hVtx->GetBinError(iVtx);

    if (nEvents==0) continue;

    for (Int_t iEta=0; iEta<=hEtaVsVtx->GetNbinsY(); iEta++) {
      Float_t value = hEtaVsVtx->GetBinContent(iVtx, iEta) / nEvents;
      if (value==0) continue;
      Float_t error = hEtaVsVtx->GetBinError(iVtx, iEta)/nEvents;
      error = TMath::Sqrt(TMath::Power(hEtaVsVtx->GetBinError(iVtx, iEta)/
				       hEtaVsVtx->GetBinContent(iVtx, iEta),2) +
      			  TMath::Power(nEventsError/nEvents,2));
      hEtaVsVtx->SetBinContent(iVtx, iEta, value);
      hEtaVsVtx->SetBinError(iVtx,   iEta, error);
    }

    //debug
    for (Int_t iEta=0; iEta<=hEtaVsVtxCheck->GetNbinsY(); iEta++) {
      Float_t value = hEtaVsVtxCheck->GetBinContent(iVtx, iEta) / nEvents;
      if (value==0) continue;
      Float_t error = hEtaVsVtxCheck->GetBinError(iVtx, iEta)/nEvents;
      error = TMath::Sqrt(TMath::Power(hEtaVsVtxCheck->GetBinError(iVtx, iEta)/
				       hEtaVsVtxCheck->GetBinContent(iVtx, iEta),2) +
      			  TMath::Power(nEventsError/nEvents,2));
      hEtaVsVtxCheck->SetBinContent(iVtx, iEta, value);
      hEtaVsVtxCheck->SetBinError(iVtx,   iEta, error);
    }
  }

  // then take the wieghted average for each eta
  // is this the right way to do it???
  for (Int_t iEta=0; iEta<=hEtaVsVtx->GetNbinsY(); iEta++) {
    Float_t sum           = 0;
    Float_t sumError2     = 0;
    Int_t   nMeasurements = 0;    

    Float_t sumXw = 0;
    Float_t sumW  = 0;
    
    // do we have several histograms for different vertex positions?
    Int_t vertexBinWidth = hVtx->GetNbinsX() / (kVertexBinning-1);
    for (Int_t vertexPos=0; vertexPos<kVertexBinning; ++vertexPos)
    {
      Int_t vertexBinBegin = 0;
      Int_t vertexBinEnd = hVtx->GetNbinsX();

      // the first histogram is always for the whole vertex range
      if (vertexPos > 0)
      {
        vertexBinBegin = vertexBinWidth * (vertexPos-1);
        vertexBinEnd = vertexBinBegin + vertexBinWidth;
      }

      for (Int_t iVtx=vertexBinBegin; iVtx<=vertexBinEnd; iVtx++) {
        if (hVtx->GetBinContent(iVtx)==0)             continue;
        if (hEtaVsVtxCheck->GetBinContent(iVtx, iEta)==0) continue;

        Float_t w = 1/TMath::Power(hEtaVsVtx->GetBinError(iVtx, iEta),2);
        sumXw = sumXw + hEtaVsVtxCheck->GetBinContent(iVtx, iEta)*w;
        sumW  = sumW + w;

        sum = sum + hEtaVsVtxCheck->GetBinContent(iVtx, iEta);
        sumError2 = sumError2 + TMath::Power(hEtaVsVtxCheck->GetBinError(iVtx, iEta),2);      
        nMeasurements++;
      }
      Float_t dndeta = 0;
      Float_t error  = 0;

      if (nMeasurements!=0) {
        dndeta = sum/Float_t(nMeasurements);
        error  = TMath::Sqrt(sumError2)/Float_t(nMeasurements);

        dndeta = sumXw/sumW;
        error  = 1/TMath::Sqrt(sumW);

        dndeta = dndeta/hdNdEta[vertexPos]->GetBinWidth(iEta);
        error  = error/hdNdEta[vertexPos]->GetBinWidth(iEta);

        hdNdEta[vertexPos]->SetBinContent(iEta, dndeta);
        hdNdEta[vertexPos]->SetBinError(iEta, error);
      }
    }
  }
}

//____________________________________________________________________
void
dNdEtaAnalysis::SaveHistograms() {

  gDirectory->mkdir(GetName());
  gDirectory->cd(GetName());
  
  hEtaVsVtx  ->Write();
  hEtaVsVtxUncorrected->Write();
  hVtx       ->Write();
  for (Int_t i=0; i<kVertexBinning; ++i)
    hdNdEta[i]    ->Write();

  gDirectory->cd("../");
}

//____________________________________________________________________
void dNdEtaAnalysis::DrawHistograms()
{
  TCanvas* canvas = new TCanvas("dNdEtaAnalysis", "dNdEtaAnalysis", 1200, 800);
  canvas->Divide(3, 2);

  canvas->cd(1);
  if (hEtaVsVtx)
    hEtaVsVtx->Draw("COLZ");

  canvas->cd(2);
  if (hEtaVsVtxCheck)
    hEtaVsVtxCheck->Draw("COLZ");

  canvas->cd(3);
  if (hEtaVsVtxUncorrected)
    hEtaVsVtxUncorrected->Draw("COLZ");

  canvas->cd(4);
  TH2F* clone = (TH2F*) hEtaVsVtxCheck->Clone("clone");
  clone->Divide(hEtaVsVtx);
  clone->Draw("COLZ");

  canvas->cd(5);
  if (hVtx)
    hVtx->Draw();

  canvas->cd(6);
  if (hdNdEta[0])
    hdNdEta[0]->Draw();

    // histograms for different vertices?
  if (kVertexBinning > 0)
  {
      // doesnt work, but i dont get it, giving up...
    /*TCanvas* canvas2 =*/ new TCanvas("dNdEtaAnalysisVtx", "dNdEtaAnalysisVtx", 450, 450);
    //Int_t yPads = (Int_t) TMath::Ceil(((Double_t) kVertexBinning - 1) / 2);
    //printf("%d\n", yPads);
    //canvas2->Divide(2, yPads);

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    for (Int_t i=1; i<kVertexBinning; ++i)
    {
      //canvas2->cd(i-1);
      //printf("%d\n", i);
      if (hdNdEta[i])
      {
        hdNdEta[i]->SetLineColor(i);
        hdNdEta[i]->Draw((i == 1) ? "" : "SAME");
        legend->AddEntry(hdNdEta[i], Form("Vtx Bin %d", i-1));
      }
    }

    legend->Draw();
  }
}

Long64_t dNdEtaAnalysis::Merge(TCollection* list)
{
  // Merges a list of dNdEtaAnalysis objects with this one.
  // This is needed for PROOF.
  // Returns the number of merged objects (including this)

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // sub collections
  const Int_t nCollections = kVertexBinning + 3;
  TList* collections[nCollections];
  for (Int_t i=0; i<nCollections; ++i)
    collections[i] = new TList;

  Int_t count = 0;
  while ((obj = iter->Next()))
  {
    dNdEtaAnalysis* entry = dynamic_cast<dNdEtaAnalysis*> (obj);
    if (entry == 0)
      continue;

    collections[0]->Add(entry->hEtaVsVtx);
    collections[1]->Add(entry->hEtaVsVtxUncorrected);
    collections[2]->Add(entry->hVtx);

    for (Int_t i=0; i<kVertexBinning; ++i)
      collections[3+i]->Add(entry->hdNdEta[i]);

    ++count;
  }

  hEtaVsVtx->Merge(collections[0]);
  hEtaVsVtxUncorrected->Merge(collections[1]);
  hVtx->Merge(collections[2]);
  for (Int_t i=0; i<kVertexBinning; ++i)
    hdNdEta[i]->Merge(collections[3+i]);

  for (Int_t i=0; i<nCollections; ++i)
    delete collections[i];

  return count+1;
}
