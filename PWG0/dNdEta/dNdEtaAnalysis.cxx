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
dNdEtaAnalysis::FillTrack(Float_t vtx, Float_t eta) {
  hEtaVsVtxUncorrected->Fill(vtx,eta);
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

  for (Int_t iEta=0; iEta<=hEtaVsVtx->GetNbinsY(); iEta++)
  {
    // do we have several histograms for different vertex positions?
    Int_t vertexBinWidth = hVtx->GetNbinsX() / (kVertexBinning-1);
    for (Int_t vertexPos=0; vertexPos<kVertexBinning; ++vertexPos)
    {
      Int_t vertexBinBegin = 1;
      Int_t vertexBinEnd = hVtx->GetNbinsX() + 1;

      // the first histogram is always for the whole vertex range
      if (vertexPos > 0)
      {
        vertexBinBegin = 1 + vertexBinWidth * (vertexPos-1);
        vertexBinEnd = vertexBinBegin + vertexBinWidth;
      }

      Float_t totalEvents = hVtx->Integral(vertexBinBegin, vertexBinEnd - 1);
      if (totalEvents == 0)
      {
        printf("WARNING: No events for hist %d %d %d\n", vertexPos, vertexBinBegin, vertexBinEnd);
        continue;
      }

      Float_t sum = 0;
      Float_t sumError2 = 0;
      for (Int_t iVtx = vertexBinBegin; iVtx < vertexBinEnd; iVtx++)
      {
        if (hEtaVsVtx->GetBinContent(iVtx, iEta) != 0)
        {
          sum = sum + hEtaVsVtx->GetBinContent(iVtx, iEta);
          sumError2 = sumError2 + TMath::Power(hEtaVsVtx->GetBinError(iVtx, iEta),2);
        }
      }

      Float_t dndeta = sum / totalEvents;
      Float_t error  = TMath::Sqrt(sumError2) / totalEvents;

      dndeta = dndeta/hdNdEta[vertexPos]->GetBinWidth(iEta);
      error  = error/hdNdEta[vertexPos]->GetBinWidth(iEta);

      hdNdEta[vertexPos]->SetBinContent(iEta, dndeta);
      hdNdEta[vertexPos]->SetBinError(iEta, error);
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

void dNdEtaAnalysis::LoadHistograms()
{
  gDirectory->cd(GetName());

  hEtaVsVtx = dynamic_cast<TH2F*> (gDirectory->Get(hEtaVsVtx->GetName()));
  hEtaVsVtxUncorrected = dynamic_cast<TH2F*> (gDirectory->Get(hEtaVsVtxUncorrected->GetName()));

  hVtx = dynamic_cast<TH1D*> (gDirectory->Get(hVtx->GetName()));

  for (Int_t i=0; i<kVertexBinning; ++i)
    hdNdEta[i] = dynamic_cast<TH1D*> (gDirectory->Get(hdNdEta[i]->GetName()));

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

    for (Int_t i=0; i<kVertexBinning; ++i)
    {
      //canvas2->cd(i-1);
      //printf("%d\n", i);
      if (hdNdEta[i])
      {
        hdNdEta[i]->SetLineColor(i+1);
        hdNdEta[i]->Draw((i == 0) ? "" : "SAME");
        legend->AddEntry(hdNdEta[i], (i == 0) ? "Vtx All" : Form("Vtx Bin %d", i-1));
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
