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
  // constructor

  fEtaVsVtx  = new TH2F(Form("%s_eta_vs_vtx", name),"",80,-20,20,120,-6,6);
  fEtaVsVtx->SetXTitle("vtx z [cm]");
  fEtaVsVtx->SetYTitle("#eta");

  fEtaVsVtxUncorrected = dynamic_cast<TH2F*> (fEtaVsVtx->Clone(Form("%s_eta_vs_vtx_uncorrected", name)));
  fVtx       = fEtaVsVtx->ProjectionX(Form("%s_vtx", name));
  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i]    = fEtaVsVtx->ProjectionY(Form("%s_dNdEta_%d", name, i));
    fdNdEta[i]->SetYTitle("dN/d#eta");
  }

  fEtaVsVtx->Sumw2();
  fVtx->Sumw2();
}

//____________________________________________________________________
void dNdEtaAnalysis::FillTrack(Float_t vtx, Float_t eta)
{
  // fills a track into the histograms

  fEtaVsVtxUncorrected->Fill(vtx,eta);
}

//____________________________________________________________________
void dNdEtaAnalysis::FillEvent(Float_t vtx)
{
  // fills an event into the histograms

  fVtx->Fill(vtx);
}

//____________________________________________________________________
void dNdEtaAnalysis::Finish(dNdEtaCorrection* correction)
{
  // correct with correction values if available

  // TODO what do we do with the error?
  if (!correction)
    printf("INFO: No correction applied\n");

  // this can be replaced by TH2F::Divide if we agree that the binning will be always the same
  for (Int_t iVtx=0; iVtx<=fEtaVsVtxUncorrected->GetNbinsX(); iVtx++)
  {
    for (Int_t iEta=0; iEta<=fEtaVsVtxUncorrected->GetNbinsY(); iEta++)
    {
      Float_t correctionValue = 1;
      if (correction)
        correctionValue = correction->GetCorrection(fEtaVsVtxUncorrected->GetXaxis()->GetBinCenter(iVtx), fEtaVsVtxUncorrected->GetYaxis()->GetBinCenter(iEta));

      Float_t value = fEtaVsVtxUncorrected->GetBinContent(iVtx, iEta);
      Float_t error = fEtaVsVtxUncorrected->GetBinError(iVtx, iEta);

      Float_t correctedValue = value * correctionValue;
      Float_t correctedError = error * correctionValue;

      if (correctedValue != 0)
      {
        fEtaVsVtx->SetBinContent(iVtx, iEta, correctedValue);
        fEtaVsVtx->SetBinError(iVtx, iEta, correctedError);
      }
    }
  }

  for (Int_t iEta=0; iEta<=fEtaVsVtx->GetNbinsY(); iEta++)
  {
    // do we have several histograms for different vertex positions?
    Int_t vertexBinWidth = fVtx->GetNbinsX() / (kVertexBinning-1);
    for (Int_t vertexPos=0; vertexPos<kVertexBinning; ++vertexPos)
    {
      Int_t vertexBinBegin = 1;
      Int_t vertexBinEnd = fVtx->GetNbinsX() + 1;

      // the first histogram is always for the whole vertex range
      if (vertexPos > 0)
      {
        vertexBinBegin = 1 + vertexBinWidth * (vertexPos-1);
        vertexBinEnd = vertexBinBegin + vertexBinWidth;
      }

      Float_t totalEvents = fVtx->Integral(vertexBinBegin, vertexBinEnd - 1);
      if (totalEvents == 0)
      {
        printf("WARNING: No events for hist %d %d %d\n", vertexPos, vertexBinBegin, vertexBinEnd);
        continue;
      }

      Float_t sum = 0;
      Float_t sumError2 = 0;
      for (Int_t iVtx = vertexBinBegin; iVtx < vertexBinEnd; iVtx++)
      {
        if (fEtaVsVtx->GetBinContent(iVtx, iEta) != 0)
        {
          sum = sum + fEtaVsVtx->GetBinContent(iVtx, iEta);
          sumError2 = sumError2 + TMath::Power(fEtaVsVtx->GetBinError(iVtx, iEta),2);
        }
      }

      Float_t dndeta = sum / totalEvents;
      Float_t error  = TMath::Sqrt(sumError2) / totalEvents;

      dndeta = dndeta/fdNdEta[vertexPos]->GetBinWidth(iEta);
      error  = error/fdNdEta[vertexPos]->GetBinWidth(iEta);

      fdNdEta[vertexPos]->SetBinContent(iEta, dndeta);
      fdNdEta[vertexPos]->SetBinError(iEta, error);
    }
  }
}

//____________________________________________________________________
void dNdEtaAnalysis::SaveHistograms()
{
  // save the histograms to a directory with the name of this class (retrieved from TNamed)

  gDirectory->mkdir(GetName());
  gDirectory->cd(GetName());

  fEtaVsVtx  ->Write();
  fEtaVsVtxUncorrected->Write();
  fVtx       ->Write();
  for (Int_t i=0; i<kVertexBinning; ++i)
    fdNdEta[i]    ->Write();

  gDirectory->cd("../");
}

void dNdEtaAnalysis::LoadHistograms()
{
  // loads the histograms from a directory with the name of this class (retrieved from TNamed)

  gDirectory->cd(GetName());

  fEtaVsVtx = dynamic_cast<TH2F*> (gDirectory->Get(fEtaVsVtx->GetName()));
  fEtaVsVtxUncorrected = dynamic_cast<TH2F*> (gDirectory->Get(fEtaVsVtxUncorrected->GetName()));

  fVtx = dynamic_cast<TH1D*> (gDirectory->Get(fVtx->GetName()));

  for (Int_t i=0; i<kVertexBinning; ++i)
    fdNdEta[i] = dynamic_cast<TH1D*> (gDirectory->Get(fdNdEta[i]->GetName()));

  gDirectory->cd("../");
}

//____________________________________________________________________
void dNdEtaAnalysis::DrawHistograms()
{
  // draws the histograms

  TCanvas* canvas = new TCanvas("dNdEtaAnalysis", "dNdEtaAnalysis", 800, 800);
  canvas->Divide(2, 2);

  canvas->cd(1);
  if (fEtaVsVtx)
    fEtaVsVtx->Draw("COLZ");

  canvas->cd(2);
  if (fEtaVsVtxUncorrected)
    fEtaVsVtxUncorrected->Draw("COLZ");

  canvas->cd(3);
  if (fVtx)
    fVtx->Draw();

  canvas->cd(4);
  if (fdNdEta[0])
    fdNdEta[0]->Draw();

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
      if (fdNdEta[i])
      {
        fdNdEta[i]->SetLineColor(i+1);
        fdNdEta[i]->Draw((i == 0) ? "" : "SAME");
        legend->AddEntry(fdNdEta[i], (i == 0) ? "Vtx All" : Form("Vtx Bin %d", i-1));
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

    collections[0]->Add(entry->fEtaVsVtx);
    collections[1]->Add(entry->fEtaVsVtxUncorrected);
    collections[2]->Add(entry->fVtx);

    for (Int_t i=0; i<kVertexBinning; ++i)
      collections[3+i]->Add(entry->fdNdEta[i]);

    ++count;
  }

  fEtaVsVtx->Merge(collections[0]);
  fEtaVsVtxUncorrected->Merge(collections[1]);
  fVtx->Merge(collections[2]);
  for (Int_t i=0; i<kVertexBinning; ++i)
    fdNdEta[i]->Merge(collections[3+i]);

  for (Int_t i=0; i<nCollections; ++i)
    delete collections[i];

  return count+1;
}
