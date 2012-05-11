#include "AliHighPtDeDxSpectra.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif

using namespace std;

ClassImp(AliHighPtDeDxSpectra);

//
// AliHighPtDeDxSpectra class
//
// This class contains the AliHighPtDeDxSpectra information 
//

//_________________________________________________________
AliHighPtDeDxSpectra::AliHighPtDeDxSpectra():
  TNamed(),
  fDebugLevel(0),
  fUseMcNoVtxCorrection(kTRUE),
  fUseFittedEfficiency(kTRUE),
  fUseBinCorrection(kTRUE),
  fNevents(0),
  hPt(0x0),
  hNevents(0x0),
  hTriggerEfficiency(0x0),
  hEfficiency(0x0),
  hMcNoVtxCorrection(0x0),
  hBinCorrection(0x0),
  hMeanPt(0x0),
  fEfficiency(0x0),
  fBinFit(0x0)
{
  // default constructor - do not use
}

//_________________________________________________________
AliHighPtDeDxSpectra::AliHighPtDeDxSpectra(const char* name, const char* title):
  TNamed(name, title),
  fDebugLevel(0),
  fUseMcNoVtxCorrection(kTRUE),
  fUseFittedEfficiency(kTRUE),
  fUseBinCorrection(kTRUE),
  fNevents(0),
  hPt(0x0),
  hNevents(0x0),
  hTriggerEfficiency(0x0),
  hEfficiency(0x0),
  hMcNoVtxCorrection(0x0),
  hBinCorrection(0x0),
  hMeanPt(0x0),
  fEfficiency(0x0),
  fBinFit(0x0)
{
  // named constructor
}

//_________________________________________________________
AliHighPtDeDxSpectra::~AliHighPtDeDxSpectra()
{
  // histograms
  delete hPt;                // pt spectrum for 
  delete hNevents;           // events
  delete hTriggerEfficiency; // trigger efficency
  delete hEfficiency;        // track correction
  delete hMcNoVtxCorrection;
  delete hBinCorrection;     // bin correction
  delete hMeanPt;            // mean pt of data
  delete fEfficiency;
  delete fBinFit;
}

//_________________________________________________________
TH1D* AliHighPtDeDxSpectra::ConstructBinCorrection(TH1D* histPt, TProfile* histMeanPt)
{
  TH1D* histCorr = dynamic_cast<TH1D*>(histPt->Clone("hBinCorrection"));
  histCorr->SetDirectory(0);
  histCorr->Reset();

  fBinFit = new TF1("fBinFit", "[0]*(1.0+x/[1])**(-[2])", 0.0, 50);
  fBinFit->SetParameters(6, 0.5, 4);
  histPt->Fit(fBinFit, "0N", "", 3.0, 50.0);
  
  const Int_t nPtBins = histPt->GetNbinsX();
  for(Int_t bin = 1; bin < nPtBins; bin++) {
    
    Float_t meanPt = histMeanPt->GetBinContent(bin);
    Float_t binPt  = histPt->GetXaxis()->GetBinCenter(bin);
    Float_t correction = 0;
    if(meanPt>0)
      correction = fBinFit->Eval(meanPt)/fBinFit->Eval(binPt); 
    if(fDebugLevel>5)
      cout << "<pT>: " << meanPt << ", bin pT: " << binPt << ", correction: " << correction << endl;
    histCorr->SetBinContent(bin, correction);
    histCorr->SetBinError(bin, 0.0);
  }

  return histCorr;
}

//_________________________________________________________
TH1D* AliHighPtDeDxSpectra::ConstructTriggerEfficiency(AliHighPtDeDxMc* mc)
{
  TH1D* histEff = dynamic_cast<TH1D*>(mc->GetHistNeventsMcTrig()->Clone("hTriggerEfficiency"));
  histEff->SetDirectory(0);
  histEff->Divide(histEff,  mc->GetHistNeventsMc(), 1, 1, "B");

  if(fDebugLevel>1) {
    
    cout << "Trigger efficiency: " 
	 << Form("%.1f %% (No vtx), %.1f %% (vtx)", 
		 Float_t(histEff->GetBinContent(1)*100.0),
		 Float_t(histEff->GetBinContent(2)*100.0))
	 << endl;
  }
    
  return histEff;
}

//_________________________________________________________
TH1D* AliHighPtDeDxSpectra::GetEventCount(AliHighPtDeDxData* data, AliHighPtDeDxMc* mc)
{
  // Correct
  TH1D* histEvents = dynamic_cast<TH1D*>(data->GetHistNevents()->Clone("hNevents"));
  
  Double_t noVtx = histEvents->GetBinContent(1);
  //  Double_t vtx   = histEvents->GetBinContent(2);
  Double_t vtx   = histEvents->GetBinContent(2);

  if(fDebugLevel>1) {
    cout << "BEFORE correction: " << noVtx << " (no vtx), " << vtx << " (vtx)" << endl;
  }

  if(!fUseMcNoVtxCorrection) {

    TH1D* histVtxStatus = data->GetHistVtxStatus();
    
    Double_t outsideCut = histVtxStatus->GetBinContent(2);
    Double_t insideCut  = histVtxStatus->GetBinContent(3);
    Double_t ratio      = insideCut/(insideCut+outsideCut);
    
    if(fDebugLevel>1) {
      cout << "Inside/Outside/Ratio(I/(I+O)):" << insideCut << "/" << outsideCut << "/" << ratio*100 << " %%" << endl;
    }
    
    noVtx *= ratio;
    histEvents->SetBinContent(1, noVtx);
  
    if(fDebugLevel>1) {
      cout << ", AFTER correction: " << noVtx << endl;
    }
  } else {

    hMcNoVtxCorrection = dynamic_cast<TH1D*>(mc->GetHistNeventsMcTrig()->Clone("hMcEventCorrection"));
    hMcNoVtxCorrection->Divide(mc->GetHistNevents());

    Double_t scaleNoVtx = hMcNoVtxCorrection->GetBinContent(1);
    Double_t scaleVtx   = hMcNoVtxCorrection->GetBinContent(2);

    if(fDebugLevel>1) {
      cout << "Mc event correction: " 
	   << Form("%.1f %% (No vtx), %.1f %% (vtx - not corrected)", 
		   Float_t(scaleNoVtx*100.0),
		   Float_t(scaleVtx*100.0))
	   << endl;
    }

    noVtx *= scaleNoVtx;
    //    vtx  *= scaleVtx;
    histEvents->SetBinContent(1, noVtx);    
    //    histEvents->SetBinContent(2, vtx);    
  }
  
  if(fDebugLevel>1) {
    cout << "AFTER correction: " << noVtx << " (no vtx), " << vtx << " (vtx)" << endl;
  }

  return histEvents;
}


//_________________________________________________________
TH1D* AliHighPtDeDxSpectra::ConstructTrackCorrection(AliHighPtDeDxMc* mc)
{
  TH1D* histEff = dynamic_cast<TH1D*>(mc->GetHistPt()->Clone("hEfficiency"));
  histEff->SetDirectory(0);
  histEff->Divide(histEff,  mc->GetHistPtMc(0, 0), 1, 1, "B");
  //  histEff->Divide(histEff,  mc->GetHistPtMc());
  
  return histEff;
}

//_________________________________________________________
TF1* AliHighPtDeDxSpectra::FitEfficiency(TH1D* histEff)
{
  TF1* func = new TF1("fEfficiency", "[0]*(1-[1]/x)", 0, 50);
  
  histEff->Fit(func, "0N", "", 3, 50); 

  return func;
}

//_________________________________________________________
void AliHighPtDeDxSpectra::GetCorrectedSpectra(AliHighPtDeDxData* data,
					       AliHighPtDeDxMc* mc)
{
  //
  // Noralize spectra
  //

  //
  // Extract the corrected number of events
  //
  Double_t fNeventsVtx = data->GetHistNevents()->GetEntries() - data->GetHistNevents()->GetBinContent(1);
  // hNevents = GetEventCount(data, mc);
  // hTriggerEfficiency = ConstructTriggerEfficiency(mc);
  // hNevents->Divide(hTriggerEfficiency);

  // if(fDebugLevel>1) {
  //   cout << "Events after efficiency: " << hNevents->GetBinContent(1)
  // 	 << " (no vtx) " << hNevents->GetBinContent(2) << " (vtx)" << endl;
  // }

  //  fNevents = hNevents->GetBinContent(1) + hNevents->GetBinContent(2);
  fNevents = fNeventsVtx/0.85;

  //
  // Extract the track correction
  //
  hEfficiency = ConstructTrackCorrection(mc);

  hPt = dynamic_cast<TH1D*>(data->GetHistPt()->Clone("hPt"));
  hPt->SetDirectory(0);

  Double_t etaRange = data->GetEtaHigh() - data->GetEtaLow();

  if(fDebugLevel>1) {
    cout << "N events (corrected): " << fNevents << endl;
    cout << "Eta range: " << etaRange << endl;

    cout << endl << "fNeventsVtx(" << fNeventsVtx << ")/fNevents(" << fNevents << ") = " 
	 << fNeventsVtx/fNevents*100 << " %" << endl << endl;
  }

  hPt->Scale(1.0/etaRange);
  NormalizePt(hPt);
  hPt->Scale(1.0/fNevents);
  
  if (fUseFittedEfficiency) {

    fEfficiency = FitEfficiency(hEfficiency);
    hPt->Divide(fEfficiency);
  } else
    hPt->Divide(hEfficiency);

  if (fUseBinCorrection) {

    hBinCorrection = ConstructBinCorrection(hPt, data->GetHistMeanPt());
    hPt->Divide(hBinCorrection);
  }
}

//_________________________________________________________
void AliHighPtDeDxSpectra::NormalizePt(TH1* hist)
{
  const Int_t n = hist->GetNbinsX();

  for(Int_t bin = 1; bin <= n; bin++) {
    
    const Float_t width = hist->GetXaxis()->GetBinWidth(bin);
    hist->SetBinError(bin, hist->GetBinError(bin)/width);
    hist->SetBinContent(bin, hist->GetBinContent(bin)/width);
  }
}


// //_________________________________________________________
// void AliHighPtDeDxSpectra::FillTrackInfo(Float_t weight) 
// {

//   AliHighPtDeDxBase::FillTrackInfo(weight);
  
//   hPt->Fill(fTrackPt, weight);
// }
