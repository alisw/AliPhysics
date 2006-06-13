/* $Id$ */

#include "dNdEtaAnalysis.h"

#include <TFile.h>
#include <TH3F.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCollection.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>

#include "AlidNdEtaCorrection.h"
#include "AliPWG0Helper.h"

//____________________________________________________________________
ClassImp(dNdEtaAnalysis)

//____________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(Char_t* name, Char_t* title) :
  TNamed(name, title),
  fData(0),
  fDataUncorrected(0),
  fNEvents(0),
  fVtx(0)
{
  // constructor

  fData  = new TH3F(Form("%s_analysis", name),"dNdEtaAnalysis",80,-20,20,120,-6,6,100, 0, 10);
  fData->SetXTitle("vtx z [cm]");
  fData->SetYTitle("#eta");
  fData->SetZTitle("p_{T}");

  fDataUncorrected = dynamic_cast<TH3F*> (fData->Clone(Form("%s_analysis_uncorrected", name)));
  fVtx       = dynamic_cast<TH1D*> (fData->Project3D("x"));

  fdNdEta[0] = dynamic_cast<TH1D*> (fData->Project3D("y"));
  for (Int_t i=1; i<kVertexBinning; ++i)
  {
    fdNdEta[i]    = dynamic_cast<TH1D*> (fdNdEta[0]->Clone(Form("%s_%d", fdNdEta[0]->GetName(), i)));
    fdNdEta[i]->SetYTitle("dN/d#eta");
  }

  fData->Sumw2();
  fVtx->Sumw2();
}

//____________________________________________________________________
dNdEtaAnalysis::~dNdEtaAnalysis()
{
  // destructor

  delete fData;
  fData = 0;

  delete fDataUncorrected;
  fDataUncorrected = 0;

  delete fVtx;
  fVtx = 0;

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    delete fdNdEta[i];
    fdNdEta[i] = 0;
  }
}

//_____________________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(const dNdEtaAnalysis &c) :
  TNamed(c),
  fData(0),
  fDataUncorrected(0),
  fNEvents(0),
  fVtx(0)
{
  //
  // dNdEtaAnalysis copy constructor
  //

  ((dNdEtaAnalysis &) c).Copy(*this);
}

//_____________________________________________________________________________
dNdEtaAnalysis &dNdEtaAnalysis::operator=(const dNdEtaAnalysis &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((dNdEtaAnalysis &) c).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void dNdEtaAnalysis::Copy(TObject &c) const
{
  //
  // Copy function
  //

  dNdEtaAnalysis& target = (dNdEtaAnalysis &) c;

  target.fData = dynamic_cast<TH3F*> (fData->Clone());
  target.fDataUncorrected = dynamic_cast<TH3F*> (fDataUncorrected->Clone());
  target.fVtx = dynamic_cast<TH1D*> (fVtx->Clone());

  for (Int_t i=0; i<kVertexBinning; ++i)
    target.fdNdEta[i] = dynamic_cast<TH1D*> (fdNdEta[i]->Clone());

  target.fNEvents = fNEvents;

  TNamed::Copy((TNamed &) c);
}

//____________________________________________________________________
void dNdEtaAnalysis::FillTrack(Float_t vtx, Float_t eta, Float_t pt, Float_t weight)
{
  // fills a track into the histograms

  fDataUncorrected->Fill(vtx, eta, pt);
  fData->Fill(vtx, eta, pt, weight);
}

//____________________________________________________________________
void dNdEtaAnalysis::FillEvent(Float_t vtx, Float_t weight)
{
  // fills an event into the histograms

  fVtx->Fill(vtx, weight);

  fNEvents += weight;
}

//____________________________________________________________________
void dNdEtaAnalysis::Finish(AlidNdEtaCorrection* correction, Float_t ptCut)
{
  // correct with correction values if available

  // TODO what do we do with the error?
  if (!correction)
    printf("INFO: No correction applied\n");

  // In fData we have the track2particle and vertex reconstruction efficiency correction already applied

  // integrate over pt (with pt cut)
  fData->GetZaxis()->SetRange(fData->GetZaxis()->FindBin(ptCut), fData->GetZaxis()->GetNbins());
  TH2D* vtxVsEta = dynamic_cast<TH2D*> (fData->Project3D("yx2"));
  if (vtxVsEta == 0)
  {
    printf("ERROR: pt integration failed\n");
    return;
  }

  new TCanvas;
  vtxVsEta->Draw("COLZ");

  for (Int_t iEta=0; iEta<=vtxVsEta->GetNbinsY(); iEta++)
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
        if (vtxVsEta->GetBinContent(iVtx, iEta) != 0)
        {
          sum = sum + vtxVsEta->GetBinContent(iVtx, iEta);
          sumError2 = sumError2 + TMath::Power(vtxVsEta->GetBinError(iVtx, iEta),2);
        }
      }

      Float_t ptCutOffCorrection = correction->GetMeasuredFraction(ptCut, vtxVsEta->GetYaxis()->GetBinCenter(iEta));
      //ptCutOffCorrection = 1;
      if (ptCutOffCorrection <= 0)
      {
        printf("UNEXPECTED: ptCutOffCorrection is %f for hist %d %d %d\n", ptCutOffCorrection, vertexPos, vertexBinBegin, vertexBinEnd);
        continue;
      }

      Float_t dndeta = sum / totalEvents / ptCutOffCorrection;
      Float_t error  = TMath::Sqrt(sumError2) / totalEvents / ptCutOffCorrection;

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

  fData  ->Write();
  AliPWG0Helper::CreateProjections(fData);
  fDataUncorrected->Write();
  AliPWG0Helper::CreateProjections(fDataUncorrected);

  fVtx       ->Write();
  for (Int_t i=0; i<kVertexBinning; ++i)
    fdNdEta[i]    ->Write();

  gDirectory->cd("../");
}

void dNdEtaAnalysis::LoadHistograms()
{
  // loads the histograms from a directory with the name of this class (retrieved from TNamed)

  gDirectory->cd(GetName());

  fData = dynamic_cast<TH3F*> (gDirectory->Get(fData->GetName()));
  fDataUncorrected = dynamic_cast<TH3F*> (gDirectory->Get(fDataUncorrected->GetName()));

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
  if (fData)
    fData->Draw("COLZ");

  canvas->cd(2);
  if (fDataUncorrected)
    fDataUncorrected->Draw("COLZ");

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

    collections[0]->Add(entry->fData);
    collections[1]->Add(entry->fDataUncorrected);
    collections[2]->Add(entry->fVtx);

    for (Int_t i=0; i<kVertexBinning; ++i)
      collections[3+i]->Add(entry->fdNdEta[i]);

    ++count;
  }

  fData->Merge(collections[0]);
  fDataUncorrected->Merge(collections[1]);
  fVtx->Merge(collections[2]);
  for (Int_t i=0; i<kVertexBinning; ++i)
    fdNdEta[i]->Merge(collections[3+i]);

  for (Int_t i=0; i<nCollections; ++i)
    delete collections[i];

  return count+1;
}
