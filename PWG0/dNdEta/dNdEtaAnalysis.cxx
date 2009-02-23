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
#include <TParameter.h>

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
  fMult(0),
  fPtDist(0),
  fAnalysisMode(AliPWG0Helper::kInvalid),
  fTag()
{
  // default constructor

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i] = 0;
    fdNdEtaPtCutOffCorrected[i] = 0;
  }
}

//____________________________________________________________________
dNdEtaAnalysis::dNdEtaAnalysis(const Char_t* name, const Char_t* title, AliPWG0Helper::AnalysisMode analysisMode) :
  TNamed(name, title),
  fData(0),
  fMult(0),
  fPtDist(0),
  fAnalysisMode(analysisMode),
  fTag()
{
  // constructor

  fData = new AliCorrection("Analysis", Form("%s Analysis", title), analysisMode);

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fMult = new TH1F("TriggeredMultiplicity", "Triggered Events;raw multiplicity;entries", 1000, -0.5, 999.5);

  TH1* histForBinning = fData->GetTrackCorrection()->GetGeneratedHistogram();
  fdNdEta[0] = new TH1F("dNdEta", "dN_{ch}/d#eta;#eta;dN_{ch}/d#eta", histForBinning->GetNbinsY(), histForBinning->GetYaxis()->GetXbins()->GetArray());

  fdNdEtaPtCutOffCorrected[0] = dynamic_cast<TH1F*> (fdNdEta[0]->Clone(Form("%s_corrected", fdNdEta[0]->GetName())));

  for (Int_t i=1; i<kVertexBinning; ++i)
  {
    fdNdEta[i]    = dynamic_cast<TH1F*> (fdNdEta[0]->Clone(Form("%s_%d", fdNdEta[0]->GetName(), i)));
    fdNdEtaPtCutOffCorrected[i]    = dynamic_cast<TH1F*> (fdNdEtaPtCutOffCorrected[0]->Clone(Form("%s_%d", fdNdEtaPtCutOffCorrected[0]->GetName(), i)));
  }

  fPtDist = new TH1F("Pt", "p_{T} distribution;p_{T} [GeV/c];#frac{dN}{d#eta dp_{T}} [c/GeV]", histForBinning->GetNbinsZ(), histForBinning->GetZaxis()->GetXbins()->GetArray());

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

  if (fMult)
  {
    delete fMult;
    fMult = 0;
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
  fMult(0),
  fPtDist(0),
  fAnalysisMode(AliPWG0Helper::kInvalid),
  fTag()
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
  target.fMult = dynamic_cast<TH1F*> (fMult->Clone());

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    target.fdNdEta[i] = dynamic_cast<TH1F*> (fdNdEta[i]->Clone());
    target.fdNdEtaPtCutOffCorrected[i] = dynamic_cast<TH1F*> (fdNdEtaPtCutOffCorrected[i]->Clone());
  }

  target.fPtDist = dynamic_cast<TH1F*> (fPtDist->Clone());

  target.fAnalysisMode = fAnalysisMode;
  target.fTag = fTag;

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
void dNdEtaAnalysis::FillTriggeredEvent(Float_t n)
{
  // fills a triggered event into the histograms

  fMult->Fill(n);
}

//____________________________________________________________________
void dNdEtaAnalysis::Finish(AlidNdEtaCorrection* correction, Float_t ptCut, AlidNdEtaCorrection::CorrectionType correctionType, const char* tag)
{
  //
  // correct with the given correction values and calculate dNdEta and pT distribution
  // the corrections that are applied can be steered by the flag correctionType
  // the measured result is not used up to a multiplicity of multCut (the bin at multCut is the first that is used!)
  //

  fTag.Form("Correcting dN/deta spectrum (data: %d) >>> %s <<<. Correction type: %d, pt cut: %.2f.", (Int_t) fAnalysisMode, tag, (Int_t) correctionType, ptCut);
  Printf("\n\n%s", fTag.Data());

  // set corrections to 1
  fData->SetCorrectionToUnity();

  if (correction && correctionType != AlidNdEtaCorrection::kNone)
  {
    TH3* trackCorr = fData->GetTrackCorrection()->GetCorrectionHistogram();
    TH2* eventCorr = fData->GetEventCorrection()->GetCorrectionHistogram();

    if (correctionType >= AlidNdEtaCorrection::kTrack2Particle)
      trackCorr->Multiply(correction->GetTrack2ParticleCorrection()->GetTrackCorrection()->GetCorrectionHistogram());

    if (correctionType >= AlidNdEtaCorrection::kVertexReco)
    {
      trackCorr->Multiply(correction->GetVertexRecoCorrection()->GetTrackCorrection()->GetCorrectionHistogram());
      eventCorr->Multiply(correction->GetVertexRecoCorrection()->GetEventCorrection()->GetCorrectionHistogram());

      // set bin with multiplicity 0 to unity (correction has no meaning in this bin)
      for (Int_t i=0; i<=eventCorr->GetNbinsX()+1; i++)
        eventCorr->SetBinContent(i, 1, 1);
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

  TH2F* rawMeasured = (TH2F*) fData->GetEventCorrection()->GetMeasuredHistogram()->Clone("rawMeasured");

  fData->Multiply();

  if (correctionType >= AlidNdEtaCorrection::kVertexReco)
  {
    // There are no events with vertex that have 0 multiplicity, therefore
    //   populate bin with 0 multiplicity with the following idea:
    //     alpha = triggered events with vertex at a given vertex position / all triggered events with vertex
    //     triggered events without vertex and 0 multiplicity at a given vertex position = alpha * all triggered events with 0 multiplicity
    //   afterwards we still correct for the trigger efficiency

    //TH2* measuredEvents = fData->GetEventCorrection()->GetMeasuredHistogram();
    TH2* correctedEvents = fData->GetEventCorrection()->GetGeneratedHistogram();

    TH2* eTrig =    correction->GetVertexRecoCorrection()->GetEventCorrection()->GetGeneratedHistogram();
    TH2* eTrigVtx = correction->GetVertexRecoCorrection()->GetEventCorrection()->GetMeasuredHistogram();
    //TH1* eTrigVtx_projx = eTrigVtx->ProjectionX("eTrigVtx_projx", 2, rawMeasured->GetNbinsY()+1);
    TH1* eTrigVtx_projx = eTrigVtx->ProjectionX("eTrigVtx_projx", 2, rawMeasured->GetNbinsY()+1);

    //new TCanvas; correctedEvents->DrawCopy("TEXT");

    // start above 0 mult. bin with integration
    TH1* vertexDist = rawMeasured->ProjectionX("vertexdist_measured", 2, rawMeasured->GetNbinsY()+1);
    //new TCanvas; vertexDist->DrawCopy();

    Int_t allEventsWithVertex = (Int_t) vertexDist->Integral(0, vertexDist->GetNbinsX()+1); // include under/overflow!
    Int_t triggeredEventsWith0Mult = (Int_t) fMult->GetBinContent(1);

    Printf("%d triggered events with 0 mult. -- %d triggered events with vertex", triggeredEventsWith0Mult, allEventsWithVertex);

    TH1* kineBias = (TH1*) vertexDist->Clone("kineBias");
    kineBias->Reset();

    for (Int_t i = 1; i <= rawMeasured->GetNbinsX(); i++)
    {
      Double_t alpha = (Double_t) vertexDist->GetBinContent(i) / allEventsWithVertex;
      Double_t events = alpha * triggeredEventsWith0Mult;

      if (eTrigVtx_projx->GetBinContent(i) == 0)
        continue;

      Double_t fZ = eTrigVtx_projx->Integral(0, eTrigVtx_projx->GetNbinsX()+1) / eTrigVtx_projx->GetBinContent(i) *
        eTrig->GetBinContent(i, 1) / eTrig->Integral(0, eTrig->GetNbinsX()+1, 1, 1);
      kineBias->SetBinContent(i, fZ);

      events *= fZ;

      // multiply with trigger correction if set above
      events *= fData->GetEventCorrection()->GetCorrectionHistogram()->GetBinContent(i, 1);

      Printf("Bin %d, alpha is %.2f, fZ is %.3f, number of events with 0 mult.: %.2f", i, alpha * 100., fZ, events);

      correctedEvents->SetBinContent(i, 1, events);
    }

    //new TCanvas; correctedEvents->DrawCopy("TEXT");
    //new TCanvas; kineBias->DrawCopy();
  }

  fData->PrintInfo(ptCut);

  TH3* dataHist = fData->GetTrackCorrection()->GetGeneratedHistogram();

  // integrate multiplicity axis out (include under/overflow bins!!!)
  TH2* tmp = fData->GetEventCorrection()->GetGeneratedHistogram();

  TH1D* vertexHist = (TH1D*) tmp->ProjectionX("_px", 0, tmp->GetNbinsY() + 1, "e");

  // create pt hist
  if (fAnalysisMode == AliPWG0Helper::kTPC || fAnalysisMode == AliPWG0Helper::kTPCITS)
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

  // integrate over pt (with pt cut) (TPC, TPCITS) or multiplicity (SPD)
  Int_t ptLowBin = 1;
  if (ptCut > 0 && fAnalysisMode != AliPWG0Helper::kSPD)
    ptLowBin = dataHist->GetZaxis()->FindBin(ptCut);
    
  //new TCanvas; dataHist->DrawCopy();

  //dataHist->Sumw2();
  dataHist->GetZaxis()->SetRange(ptLowBin, dataHist->GetZaxis()->GetNbins()+1);
  printf("pt/multiplicity range %d %d\n", ptLowBin, dataHist->GetZaxis()->GetNbins()+1);
  TH2D* vtxVsEta = dynamic_cast<TH2D*> (dataHist->Project3D("yx2e"));
  
  //new TCanvas; vtxVsEta->Draw("COLZ");

  dataHist->GetZaxis()->SetRange(0, 0);
  vtxVsEta->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle());
  vtxVsEta->GetYaxis()->SetTitle(dataHist->GetYaxis()->GetTitle());

  if (vtxVsEta == 0)
  {
    printf("ERROR: pt/multiplicity integration failed\n");
    return;
  }

  //new TCanvas(tag, tag, 800, 600);  vtxVsEta->DrawCopy("COLZ");

  // clear result histograms
  for (Int_t vertexPos=0; vertexPos<kVertexBinning; ++vertexPos)
  {
    fdNdEta[vertexPos]->Reset();
    fdNdEtaPtCutOffCorrected[vertexPos]->Reset();
  }

  const Float_t vertexRangeBegin[kVertexBinning] = { -9.99,  -9.99,  0.01 };
  const Float_t vertexRangeEnd[kVertexBinning]   = {  9.99,  -0.01,  9.99 };

  for (Int_t iEta=1; iEta<=vtxVsEta->GetNbinsY(); iEta++)
  {
    // loop over vertex ranges
    for (Int_t vertexPos=0; vertexPos<kVertexBinning; ++vertexPos)
    {
      Int_t vertexBinBegin = vertexHist->GetXaxis()->FindBin(vertexRangeBegin[vertexPos]);
      Int_t vertexBinEnd   = vertexHist->GetXaxis()->FindBin(vertexRangeEnd[vertexPos]);

      const Int_t *binBegin = 0;
      const Int_t maxBins = 60;

      // adjust acceptance range
      // produce with drawPlots.C: DetermineAcceptance(...)
      if (fAnalysisMode == AliPWG0Helper::kSPD)
      {
        //const Int_t binBeginSPD[30] = { 18, 16, 15, 14, 13, 13, 12, 11, 9, 7, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1 }; // by eye
        //const Int_t binBeginSPD[30] = { -1, -1, -1, -1, 16, 14, 12, 10, 9, 7, 6, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, -1, -1, -1, -1}; // limit in correction map is 5

        //const Int_t binBegin[30] = { -1, -1, -1, 17, 15, 14, 12, 10, 8, 7, 6, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, -1, -1, -1};  // limit in correction map is 10

        //const Int_t binBeginSPD[30] = { -1, -1, -1, -1, 16, 15, 13, 11, 9, 8, 7, 6, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, -1, -1, -1, -1}; // limit 2
        const Int_t binBeginSPD[maxBins] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 15, 14, 13, 12, 11, 10, 9, 9, 8, 7, 7, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}; //limit 5

        binBegin = binBeginSPD;
      }
      else if (fAnalysisMode == AliPWG0Helper::kTPC)
      {
        //const Int_t binBeginTPC[30] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1}; // limit 5, pt cut off 0.2 mev/c
        const Int_t binBeginTPC[maxBins] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 9, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}; // limit 5

        binBegin = binBeginTPC;
      }
      else if (fAnalysisMode == AliPWG0Helper::kTPCITS)
      {
        // TODO create map
      }

      Int_t vtxBegin = 1;
      Int_t vtxEnd   = maxBins;

      if (binBegin)
      {
        vtxBegin = binBegin[iEta - 1];
        vtxEnd = 18 + 1 - binBegin[maxBins - iEta];
      }
      else
        Printf("WARNING: No acceptance applied!");
      
      // eta range not accessible
      if (vtxBegin == -1)
        continue;
      
      //Printf("%d %d | %d %d", vtxBegin, vertexHist->GetXaxis()->FindBin(GetVtxMin(eta)), vtxEnd, vertexHist->GetXaxis()->FindBin(-GetVtxMin(-eta)));
      //vtxBegin = vertexHist->GetXaxis()->FindBin(GetVtxMin(eta));
      //vtxEnd = vertexHist->GetXaxis()->FindBin(-GetVtxMin(-eta));

      //Float_t eta = vtxVsEta->GetYaxis()->GetBinCenter(iEta);
      //printf("Eta bin: %d (%f) Vertex range: %d; before: %d %d (range) %d %d (acceptance)", iEta, eta, vertexPos, vertexBinBegin, vertexBinEnd, vtxBegin, vtxEnd);
      vertexBinBegin = TMath::Max(vertexBinBegin, vtxBegin);
      vertexBinEnd =   TMath::Min(vertexBinEnd, vtxEnd);
      //Printf(" after:  %d %d", vertexBinBegin, vertexBinEnd);

      // no data for this bin
      if (vertexBinBegin > vertexBinEnd)
      {
        //Printf("Bin empty. Skipped");
        continue;
      }

      Float_t totalEvents = 0;
      Float_t sum = 0;
      Float_t sumError2 = 0;
      Float_t unusedTracks = 0;
      Float_t unusedEvents = 0;
      for (Int_t iVtx = 1; iVtx <= vtxVsEta->GetNbinsX(); iVtx++)
      {
        if (iVtx >= vertexBinBegin && iVtx <= vertexBinEnd)
        {
          if (vtxVsEta->GetBinContent(iVtx, iEta) != 0)
          {
            sum += vtxVsEta->GetBinContent(iVtx, iEta);

            if (sumError2 > 10e30)
              Printf("WARNING: sum of error2 is dangerously large - be prepared for crash... ");

            sumError2 = sumError2 + TMath::Power(vtxVsEta->GetBinError(iVtx, iEta),2);
          }
          totalEvents += vertexHist->GetBinContent(iVtx);
        }
        else
        {
          unusedTracks += vtxVsEta->GetBinContent(iVtx, iEta);
          unusedEvents += vertexHist->GetBinContent(iVtx);
        }
      }

      if (totalEvents == 0)
      {
        printf("WARNING: No events for hist %d %d %d\n", vertexPos, vertexBinBegin, vertexBinEnd);
        continue;
      }

      Float_t ptCutOffCorrection = 1;

      // find pt cut off correction factor
      if (fAnalysisMode != AliPWG0Helper::kSPD)
      {
        if (correction && ptCut > 0)
            ptCutOffCorrection = correction->GetMeasuredFraction(correctionType, ptCut, vtxVsEta->GetYaxis()->GetBinCenter(iEta), vertexBinBegin, vertexBinEnd);

        if (ptCutOffCorrection <= 0)
        {
            printf("UNEXPECTED: ptCutOffCorrection is %f for hist %d %d %d\n", ptCutOffCorrection, vertexPos, vertexBinBegin, vertexBinEnd);
            continue;
        }
      }

      //printf("Eta: %d (%f) Vertex Range: %d %d, Event Count %f, Track Sum: %f, Track Sum corrected: %f \n", iEta, vtxVsEta->GetYaxis()->GetBinCenter(iEta), vertexBinBegin, vertexBinEnd, totalEvents, sum, sum / ptCutOffCorrection);

      Int_t bin = fdNdEta[vertexPos]->FindBin(vtxVsEta->GetYaxis()->GetBinCenter(iEta));
      if (bin > 0 && bin <= fdNdEta[vertexPos]->GetNbinsX())
      {
        Float_t dndeta = sum / totalEvents;
        Float_t error  = TMath::Sqrt(sumError2) / totalEvents;

        dndeta = dndeta / fdNdEta[vertexPos]->GetBinWidth(bin);
        error  = error / fdNdEta[vertexPos]->GetBinWidth(bin);

        fdNdEta[vertexPos]->SetBinContent(bin, dndeta);
        fdNdEta[vertexPos]->SetBinError(bin, error);

        dndeta /= ptCutOffCorrection;
        error  /= ptCutOffCorrection;

        fdNdEtaPtCutOffCorrected[vertexPos]->SetBinContent(bin, dndeta);
        fdNdEtaPtCutOffCorrected[vertexPos]->SetBinError(bin, error);

        //Printf("Bin %d has dN/deta = %f +- %f; %.2f tracks %.2f events (outside acceptance: %.2f tracks, %.2f events)", bin, dndeta, error, sum, totalEvents, unusedTracks, unusedEvents);
      }
    }
  }
}

//____________________________________________________________________
Float_t dNdEtaAnalysis::GetVtxMin(Float_t eta)
{
  // helper function for the SPD acceptance
  // the function returns the beginning of the acceptance window in z vertex position as function of eta
  // to get the maximum: -GetVtxMin(-eta)

  Float_t a[2] = { -15, 0 };
  Float_t b[2] = { 0, -1.4 };
  Float_t c[2] = { 15, -2.2 };

  Float_t meanAB[2];
  meanAB[0] = (b[0] + a[0]) / 2;
  meanAB[1] = (b[1] + a[1]) / 2;

  Float_t meanBC[2];
  meanBC[0] = (c[0] + b[0]) / 2;
  meanBC[1] = (c[1] + b[1]) / 2;

  Float_t mAB = (b[1] - a[1]) / (b[0] - a[0]);
  Float_t mBC = (c[1] - b[1]) / (c[0] - b[0]);

  Float_t bAB = meanAB[1] - mAB * meanAB[0];
  Float_t bBC = meanBC[1] - mBC * meanBC[0];

  if (eta > b[1])
    return 1.0 / mAB * eta - bAB / mAB;

  return 1.0 / mBC * eta - bBC / mBC;
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

  if (fMult)
  {
    fMult->Write();
  }
  else
    printf("dNdEtaAnalysis::SaveHistograms: UNEXPECTED: fMult is 0\n");

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

  TNamed named("fTag", fTag.Data());
  named.Write();

  TParameter<Int_t> param("fAnalysisMode", fAnalysisMode);
  param.Write();

  gDirectory->cd("../");
}

void dNdEtaAnalysis::LoadHistograms(const Char_t* dir)
{
  // loads the histograms from a directory with the name of this class (retrieved from TNamed)

  if (!dir)
    dir = GetName();

  gDirectory->cd(dir);

  fData->LoadHistograms();
  fMult = dynamic_cast<TH1F*> (gDirectory->Get(fMult->GetName()));

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i] = dynamic_cast<TH1F*> (gDirectory->Get(fdNdEta[i]->GetName()));
    fdNdEtaPtCutOffCorrected[i] = dynamic_cast<TH1F*> (gDirectory->Get(fdNdEtaPtCutOffCorrected[i]->GetName()));
  }

  fPtDist = dynamic_cast<TH1F*> (gDirectory->Get(fPtDist->GetName()));

  if (dynamic_cast<TNamed*> (gDirectory->Get("fTag")))
    fTag = (dynamic_cast<TNamed*> (gDirectory->Get("fTag")))->GetTitle();

  if (dynamic_cast<TParameter<Int_t>*> (gDirectory->Get("fAnalysisMode")))
    fAnalysisMode = (AliPWG0Helper::AnalysisMode) (dynamic_cast<TParameter<Int_t>*> (gDirectory->Get("fAnalysisMode")))->GetVal();

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
      fdNdEtaPtCutOffCorrected[0]->DrawCopy();

    if (fdNdEta[0])
    {
      fdNdEta[0]->SetLineColor(kRed);
      fdNdEta[0]->DrawCopy("SAME");
    }

    canvas->cd(2);
    if (fPtDist)
      fPtDist->DrawCopy();
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

    TLegend* legend = new TLegend(0.4, 0.2, 0.6, 0.4);

    for (Int_t i=0; i<kVertexBinning; ++i)
    {
      if (fdNdEtaPtCutOffCorrected[i])
      {
        canvas2->cd();

        fdNdEtaPtCutOffCorrected[i]->SetLineColor(i+1);
        fdNdEtaPtCutOffCorrected[i]->DrawCopy((i == 0) ? "" : "SAME");
        legend->AddEntry(fdNdEtaPtCutOffCorrected[i], (i == 0) ? "Vtx All" : Form("Vtx Bin %d", i-1));
      }
      if (canvas3 && fdNdEta[i])
      {
        canvas3->cd();

        fdNdEta[i]->SetLineColor(i+1);
        fdNdEta[i]->DrawCopy((i == 0) ? "" : "SAME");
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
        clone->DrawCopy();

        clone2->Divide(fdNdEtaPtCutOffCorrected[0]);
        clone2->DrawCopy("SAME");

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
  const Int_t nCollections = 2 * kVertexBinning + 3; // 3 standalone hists, 3 arrays of size kVertexBinning
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
    collections[1]->Add(entry->fMult);
    collections[2]->Add(entry->fPtDist);

    for (Int_t i=0; i<kVertexBinning; ++i)
    {
      collections[3+i]->Add(entry->fdNdEta[i]);
      collections[3+kVertexBinning+i]->Add(entry->fdNdEtaPtCutOffCorrected[i]);
    }

    ++count;
  }

  fData->Merge(collections[0]);
  fMult->Merge(collections[1]);
  fPtDist->Merge(collections[2]);

  for (Int_t i=0; i<kVertexBinning; ++i)
  {
    fdNdEta[i]->Merge(collections[3+i]);
    fdNdEtaPtCutOffCorrected[i]->Merge(collections[3+kVertexBinning+i]);
  }

  for (Int_t i=0; i<nCollections; ++i)
    delete collections[i];

  return count+1;
}


