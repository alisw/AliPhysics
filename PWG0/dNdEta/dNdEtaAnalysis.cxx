/* $Id$ */

#include "dNdEtaAnalysis.h"

#include <TFile.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCollection.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TLine.h>

#include "AlidNdEtaCorrection.h"
#include <AliCorrection.h>
#include <AliPWG0Helper.h>
#include <AliCorrectionMatrix2D.h>
#include <AliCorrectionMatrix3D.h>

//____________________________________________________________________
ClassImp(dNdEtaAnalysis)

//____________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis() :
  TNamed(),
  fData(0),
  fPtDist(0)
{
  // default constructor

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i] = 0;
    fdNdEtaPtCutOffCorrected[i] = 0;
  }
}

//____________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(Char_t* name, Char_t* title) :
  TNamed(name, title),
  fData(0),
  fPtDist(0)
{
  // constructor

  Float_t binLimitsPt[] = {0.0, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 5.0, 10.0, 100.0};

  fData = new AliCorrection("Analysis", Form("%s Analysis", title));

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fdNdEta[0] = new TH1F("dNdEta", "dN_{ch}/d#eta;#eta;dN_{ch}/d#eta", 20, -2, 2);

  fdNdEtaPtCutOffCorrected[0] = dynamic_cast<TH1F*> (fdNdEta[0]->Clone(Form("%s_corrected", fdNdEta[0]->GetName())));

  for (Int_t i=1; i<kVertexBinning; ++i)
  {
    fdNdEta[i]    = dynamic_cast<TH1F*> (fdNdEta[0]->Clone(Form("%s_%d", fdNdEta[0]->GetName(), i)));
    fdNdEtaPtCutOffCorrected[i]    = dynamic_cast<TH1F*> (fdNdEtaPtCutOffCorrected[0]->Clone(Form("%s_%d", fdNdEtaPtCutOffCorrected[0]->GetName(), i)));
  }

  fPtDist = new TH1F("Pt", "p_{T} distribution;p_{T} [GeV/c];#frac{dN}{d#eta dp_{T}} [c/GeV]", 28, binLimitsPt);

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
dNdEtaAnalysis::~dNdEtaAnalysis()
{
  // destructor

  if (fData)
  {
    delete fData;
    fData = 0;
  }

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    if (fdNdEta[i])
    {
      delete fdNdEta[i];
      fdNdEta[i] = 0;
    }
    if (fdNdEtaPtCutOffCorrected[i])
    {
      delete fdNdEtaPtCutOffCorrected[i];
      fdNdEtaPtCutOffCorrected[i] = 0;
    }
  }

  if (fPtDist)
  {
    delete fPtDist;
    fPtDist = 0;
  }
}

//_____________________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(const dNdEtaAnalysis &c) :
  TNamed(c),
  fData(0),
  fPtDist(0)
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

  target.fData = dynamic_cast<AliCorrection*> (fData->Clone());

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    target.fdNdEta[i] = dynamic_cast<TH1F*> (fdNdEta[i]->Clone());
    target.fdNdEtaPtCutOffCorrected[i] = dynamic_cast<TH1F*> (fdNdEtaPtCutOffCorrected[i]->Clone());
  }

  target.fPtDist = dynamic_cast<TH1F*> (fPtDist->Clone());

  TNamed::Copy((TNamed &) c);
}

//____________________________________________________________________
void dNdEtaAnalysis::FillTrack(Float_t vtx, Float_t eta, Float_t pt)
{
  // fills a track into the histograms

  fData->GetTrackCorrection()->FillMeas(vtx, eta, pt);
}

//____________________________________________________________________
void dNdEtaAnalysis::FillEvent(Float_t vtx, Float_t n)
{
  // fills an event into the histograms

  fData->GetEventCorrection()->FillMeas(vtx, n);
}

//____________________________________________________________________
void dNdEtaAnalysis::Finish(AlidNdEtaCorrection* correction, Float_t ptCut, AlidNdEtaCorrection::CorrectionType correctionType)
{
  //
  // correct with the given correction values and calculate dNdEta and pT distribution
  // the corrections that are applied can be steered by the flag correctionType
  //

  // TODO put tag somewhere which corrections have been applied

  // set corrections to 1
  fData->SetCorrectionToUnity();

  if (correction && correctionType != AlidNdEtaCorrection::kNone)
  {
    TH3F* trackCorr = fData->GetTrackCorrection()->GetCorrectionHistogram();
    TH2F* eventCorr = fData->GetEventCorrection()->GetCorrectionHistogram();

    if (correctionType >= AlidNdEtaCorrection::kTrack2Particle)
      trackCorr->Multiply(correction->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetCorrectionHistogram());

    if (correctionType >= AlidNdEtaCorrection::kVertexReco)
    {
      trackCorr->Multiply(correction->GetVertexRecoCorrection()->GetTrackCorrection()->GetCorrectionHistogram());
      eventCorr->Multiply(correction->GetVertexRecoCorrection()->GetEventCorrection()->GetCorrectionHistogram());
    }

    switch (correctionType)
    {
      case AlidNdEtaCorrection::kINEL :
      {
        trackCorr->Multiply(correction->GetTriggerBiasCorrectionINEL()->GetTrackCorrection()->GetCorrectionHistogram());
        eventCorr->Multiply(correction->GetTriggerBiasCorrectionINEL()->GetEventCorrection()->GetCorrectionHistogram());
        break;
      }
      case AlidNdEtaCorrection::kNSD :
      {
        trackCorr->Multiply(correction->GetTriggerBiasCorrectionNSD()->GetTrackCorrection()->GetCorrectionHistogram());
        eventCorr->Multiply(correction->GetTriggerBiasCorrectionNSD()->GetEventCorrection()->GetCorrectionHistogram());
        break;
      }
      case AlidNdEtaCorrection::kND :
      {
        trackCorr->Multiply(correction->GetTriggerBiasCorrectionND()->GetTrackCorrection()->GetCorrectionHistogram());
        eventCorr->Multiply(correction->GetTriggerBiasCorrectionND()->GetEventCorrection()->GetCorrectionHistogram());
        break;
      }
      default : break;
    }
  }
  else
    printf("INFO: No correction applied\n");

  fData->Multiply();

  TH3F* dataHist = fData->GetTrackCorrection()->GetGeneratedHistogram();

  // integrate multiplicity axis out (include under/overflow bins!!!)
  TH2F* tmp = fData->GetEventCorrection()->GetGeneratedHistogram();
  TH1D* vertexHist = tmp->ProjectionX("_px", 0, tmp->GetNbinsY() + 1, "e");

  // create pt hist
  {
    // reset all ranges
    dataHist->GetXaxis()->SetRange(0, 0);
    dataHist->GetYaxis()->SetRange(0, 0);
    dataHist->GetZaxis()->SetRange(0, 0);

    // vtx cut
    Int_t vertexBinBegin = dataHist->GetXaxis()->FindBin(-5);
    Int_t vertexBinEnd = dataHist->GetXaxis()->FindBin(5);
    dataHist->GetXaxis()->SetRange(vertexBinBegin, vertexBinEnd);
    Float_t nEvents = vertexHist->Integral(vertexBinBegin, vertexBinEnd);

    if (nEvents > 0)
    {
      // eta cut
      dataHist->GetYaxis()->SetRange(dataHist->GetYaxis()->FindBin(-0.8), dataHist->GetYaxis()->FindBin(0.8));
      Float_t etaWidth = 1.6;

      TH1D* ptHist = dynamic_cast<TH1D*> (dataHist->Project3D("ze"));

      for (Int_t i=1; i<=fPtDist->GetNbinsX(); ++i)
      {
        Float_t binSize = fPtDist->GetBinWidth(i);
        fPtDist->SetBinContent(i, ptHist->GetBinContent(i) / binSize / nEvents / etaWidth);
        fPtDist->SetBinError(i, ptHist->GetBinError(i) / binSize / nEvents / etaWidth);
      }

      delete ptHist;
    }
    else
      printf("ERROR: nEvents is 0!\n");
  }

  // reset all ranges
  dataHist->GetXaxis()->SetRange(0, 0);
  dataHist->GetYaxis()->SetRange(0, 0);
  dataHist->GetZaxis()->SetRange(0, 0);

  // integrate over pt (with pt cut)
  Int_t ptLowBin = 1;
  if (ptCut > 0)
    ptLowBin = dataHist->GetZaxis()->FindBin(ptCut);

  dataHist->GetZaxis()->SetRange(ptLowBin, dataHist->GetZaxis()->GetNbins());
  printf("range %d %d\n", ptLowBin, dataHist->GetZaxis()->GetNbins());
  TH2D* vtxVsEta = dynamic_cast<TH2D*> (dataHist->Project3D("yx2e"));

  dataHist->GetZaxis()->SetRange(0, 0);
  vtxVsEta->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle());
  vtxVsEta->GetYaxis()->SetTitle(dataHist->GetYaxis()->GetTitle());

  if (vtxVsEta == 0)
  {
    printf("ERROR: pt integration failed\n");
    return;
  }

  const Float_t vertexRange = 4.99;

  for (Int_t iEta=1; iEta<=vtxVsEta->GetNbinsY(); iEta++)
  {
    // do we have several histograms for different vertex positions?
    Int_t vertexBinGlobalBegin = vertexHist->GetXaxis()->FindBin(-vertexRange);
    Int_t vertexBinWidth = (vertexHist->GetXaxis()->FindBin(vertexRange) - vertexBinGlobalBegin + 1) / (kVertexBinning-1);
    //printf("vertexBinGlobalBegin = %d, vertexBinWidth = %d\n", vertexBinGlobalBegin, vertexBinWidth);
    for (Int_t vertexPos=0; vertexPos<kVertexBinning; ++vertexPos)
    {
      Int_t vertexBinBegin = vertexBinGlobalBegin;
      Int_t vertexBinEnd = vertexBinGlobalBegin + vertexBinWidth * (kVertexBinning-1);

      // the first histogram is always for the whole vertex range
      if (vertexPos > 0)
      {
        vertexBinBegin = vertexBinGlobalBegin + vertexBinWidth * (vertexPos-1);
        vertexBinEnd = vertexBinBegin + vertexBinWidth;
      }

      //printf("vertexBinBegin = %d, vertexBinEnd = %d\n", vertexBinBegin, vertexBinEnd);

      Float_t totalEvents = vertexHist->Integral(vertexBinBegin, vertexBinEnd - 1);
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

      Float_t ptCutOffCorrection = 1;
      if (correction && ptCut > 0)
        ptCutOffCorrection = correction->GetMeasuredFraction(correctionType, ptCut, vtxVsEta->GetYaxis()->GetBinCenter(iEta));

      if (ptCutOffCorrection <= 0)
      {
        printf("UNEXPECTED: ptCutOffCorrection is %f for hist %d %d %d\n", ptCutOffCorrection, vertexPos, vertexBinBegin, vertexBinEnd);
        continue;
      }

      printf("Eta: %d Vertex Range: %d %d, Event Count %f, Track Sum: %f, Track Sum corrected: %f\n", iEta, vertexBinBegin, vertexBinEnd, totalEvents, sum, sum / ptCutOffCorrection);

      Float_t dndeta = sum / totalEvents;
      Float_t error  = TMath::Sqrt(sumError2) / totalEvents;

      dndeta = dndeta/fdNdEta[vertexPos]->GetBinWidth(iEta);
      error  = error/fdNdEta[vertexPos]->GetBinWidth(iEta);

      fdNdEta[vertexPos]->SetBinContent(iEta, dndeta);
      fdNdEta[vertexPos]->SetBinError(iEta, error);

      dndeta /= ptCutOffCorrection;
      error  /= ptCutOffCorrection;

      fdNdEtaPtCutOffCorrected[vertexPos]->SetBinContent(iEta, dndeta);
      fdNdEtaPtCutOffCorrected[vertexPos]->SetBinError(iEta, error);
    }
  }
}

//____________________________________________________________________
void dNdEtaAnalysis::SaveHistograms()
{
  // save the histograms to a directory with the name of this class (retrieved from TNamed)

  gDirectory->mkdir(GetName());
  gDirectory->cd(GetName());

  if (fData)
  {
    fData->SaveHistograms();
  }
  else
    printf("dNdEtaAnalysis::SaveHistograms: UNEXPECTED: fData is 0\n");

  if (fPtDist)
    fPtDist       ->Write();
  else
    printf("dNdEtaAnalysis::SaveHistograms: UNEXPECTED: fPtDist is 0\n");

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    if (fdNdEta[i])
      fdNdEta[i]->Write();
    else
      printf("dNdEtaAnalysis::SaveHistograms: UNEXPECTED: fdNdEta[%d] is 0\n", i);

    if (fdNdEtaPtCutOffCorrected[i])
      fdNdEtaPtCutOffCorrected[i]->Write();
    else
      printf("dNdEtaAnalysis::SaveHistograms: UNEXPECTED: fdNdEtaPtCutOffCorrected[%d] is 0\n", i);
  }

  gDirectory->cd("../");
}

void dNdEtaAnalysis::LoadHistograms(const Char_t* dir)
{
  // loads the histograms from a directory with the name of this class (retrieved from TNamed)

  if (!dir)
    dir = GetName();

  gDirectory->cd(dir);

  fData->LoadHistograms();

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i] = dynamic_cast<TH1F*> (gDirectory->Get(fdNdEta[i]->GetName()));
    fdNdEtaPtCutOffCorrected[i] = dynamic_cast<TH1F*> (gDirectory->Get(fdNdEtaPtCutOffCorrected[i]->GetName()));
  }

  fPtDist = dynamic_cast<TH1F*> (gDirectory->Get(fPtDist->GetName()));

  gDirectory->cd("../");
}

//____________________________________________________________________
void dNdEtaAnalysis::DrawHistograms(Bool_t simple)
{
  // draws the histograms
  
  if (!simple)
  {
    if (fData)
      fData->DrawHistograms(GetName());

    TCanvas* canvas = new TCanvas(Form("%s_dNdEtaAnalysis", GetName()), Form("%s_dNdEtaAnalysis", GetName()), 800, 400);
    canvas->Divide(2, 1);

    canvas->cd(1);
    if (fdNdEtaPtCutOffCorrected[0])
      fdNdEtaPtCutOffCorrected[0]->Draw();

    if (fdNdEta[0])
    {
      fdNdEta[0]->SetLineColor(kRed);
      fdNdEta[0]->Draw("SAME");
    }

    canvas->cd(2);
    if (fPtDist)
      fPtDist->Draw();
  }

    // histograms for different vertices?
  if (kVertexBinning > 0)
  {
      // doesnt work, but i dont get it, giving up...
    TCanvas* canvas2 = new TCanvas(Form("%s_dNdEtaAnalysisVtx", GetName()), Form("%s_dNdEtaAnalysisVtx", GetName()), 450, 450);
    TCanvas* canvas3 = 0;
    if (!simple)
      canvas3 = new TCanvas(Form("%s_dNdEtaAnalysisVtx_noptcutoff", GetName()), Form("%s_dNdEtaAnalysisVtx_noptcutoff", GetName()), 450, 450);

    //Int_t yPads = (Int_t) TMath::Ceil(((Double_t) kVertexBinning - 1) / 2);
    //printf("%d\n", yPads);
    //canvas2->Divide(2, yPads);

    TLegend* legend = new TLegend(0.4, 0.7, 0.6, 0.9);

    for (Int_t i=0; i<kVertexBinning; ++i)
    {
      if (fdNdEtaPtCutOffCorrected[i])
      {
        canvas2->cd();

        fdNdEtaPtCutOffCorrected[i]->SetLineColor(i+1);
        fdNdEtaPtCutOffCorrected[i]->Draw((i == 0) ? "" : "SAME");
        legend->AddEntry(fdNdEtaPtCutOffCorrected[i], (i == 0) ? "Vtx All" : Form("Vtx Bin %d", i-1));
      }
      if (canvas3 && fdNdEta[i])
      {
        canvas3->cd();

        fdNdEta[i]->SetLineColor(i+1);
        fdNdEta[i]->Draw((i == 0) ? "" : "SAME");
      }
    }

    canvas2->cd();
    legend->Draw();
    canvas2->SaveAs(Form("%s_%s.gif", canvas2->GetName(), GetName()));

    if (canvas3)
    {
      canvas3->cd();
      legend->Draw();
    }
  }
  
  if (kVertexBinning == 3)
  {
     TH1* clone = dynamic_cast<TH1*> (fdNdEtaPtCutOffCorrected[1]->Clone("clone"));
     TH1* clone2 = dynamic_cast<TH1*> (fdNdEtaPtCutOffCorrected[2]->Clone("clone2"));
     
     if (clone && clone2)
     {
        TCanvas* canvas4 = new TCanvas(Form("%s_dNdEtaAnalysisVtxRatios", GetName()), Form("%s_dNdEtaAnalysisVtxRatios", GetName()), 450, 450);
    
        clone->Divide(fdNdEtaPtCutOffCorrected[0]);
        clone->GetYaxis()->SetRangeUser(0.95, 1.05);
        clone->Draw();
        
        clone2->Divide(fdNdEtaPtCutOffCorrected[0]);
        clone2->Draw("SAME");

        TLine* line = new TLine(-1, 1, 1, 1);
        line->Draw();

        canvas4->SaveAs(Form("%s_%s.gif", canvas4->GetName(), GetName()));
     }
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
  const Int_t nCollections = 2 * kVertexBinning + 2; // 2 standalone hists, two arrays of size kVertexBinning
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
    collections[1]->Add(entry->fPtDist);

    for (Int_t i=0; i<kVertexBinning; ++i)
    {
      collections[2+i]->Add(entry->fdNdEta[i]);
      collections[2+kVertexBinning+i]->Add(entry->fdNdEtaPtCutOffCorrected[i]);
    }

    ++count;
  }

  fData->Merge(collections[0]);
  fPtDist->Merge(collections[1]);
  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i]->Merge(collections[2+i]);
    fdNdEta[i]->Merge(collections[2+kVertexBinning+i]);
  }

  for (Int_t i=0; i<nCollections; ++i)
    delete collections[i];

  return count+1;
}
