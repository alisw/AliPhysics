#include "AliHighPtDeDxMc.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif

using namespace std;

ClassImp(AliHighPtDeDxMc);

//
// AliHighPtDeDxMc class
//
// This class contains the AliHighPtDeDxMc information 
//

//_________________________________________________________
AliHighPtDeDxMc::AliHighPtDeDxMc():
  AliHighPtDeDxBase(),
  fTrackChargeMc(0),
  fTrackEtaMc(-999),
  fTrackPtMc(-1),
  hVtxStatusMc(0x0),
  hNeventsMc(0x0),
  hNeventsMcTrig(0x0),
  hPtMc(0x0),
  hPtMcNeg(0x0),
  hPtMcPos(0x0),
  hPtPiMc(0x0),
  hPtPiMcNeg(0x0),
  hPtPiMcPos(0x0),
  hPtKMc(0x0),
  hPtKMcNeg(0x0),
  hPtKMcPos(0x0),
  hPtPMc(0x0),
  hPtPMcNeg(0x0),
  hPtPMcPos(0x0),
  hMeanPtMc(0x0)
{
  // default constructor - do not use
}

//_________________________________________________________
AliHighPtDeDxMc::AliHighPtDeDxMc(const char* name, const char* title):
  AliHighPtDeDxBase(name, title),
  fTrackChargeMc(0),
  fTrackEtaMc(-999),
  fTrackPtMc(-1),
  hVtxStatusMc(0x0),
  hNeventsMc(0x0),
  hNeventsMcTrig(0x0),
  hPtMc(0x0),
  hPtMcNeg(0x0),
  hPtMcPos(0x0),
  hPtPiMc(0x0),
  hPtPiMcNeg(0x0),
  hPtPiMcPos(0x0),
  hPtKMc(0x0),
  hPtKMcNeg(0x0),
  hPtKMcPos(0x0),
  hPtPMc(0x0),
  hPtPMcNeg(0x0),
  hPtPMcPos(0x0),
  hMeanPtMc(0x0)
{
  // named constructor
}

//_________________________________________________________
AliHighPtDeDxMc::~AliHighPtDeDxMc()
{
  delete hVtxStatusMc;
  delete hNeventsMc;
  delete hNeventsMcTrig;
  delete hPtMc;
  delete hPtMcNeg;
  delete hPtMcPos;
  delete hPtPiMc;
  delete hPtPiMcNeg;
  delete hPtPiMcPos;
  delete hPtKMc;
  delete hPtKMcNeg;
  delete hPtKMcPos;
  delete hPtPMc;
  delete hPtPMcNeg;
  delete hPtPMcPos;
  delete hMeanPtMc;
}

//_________________________________________________________
void AliHighPtDeDxMc::Init(Int_t nPtBins, Double_t* ptBins)
{
  //
  // Create histograms
  //

  AliHighPtDeDxBase::Init(nPtBins, ptBins);
  
  hVtxStatusMc = new TH1D("hVtxStatus", "Vtx status (MC) - No Vtx = -1, Vtx outside cut = 0, Vtx inside = 1",
			3, -1.5, 1.5);
  hVtxStatusMc->Sumw2();
  hVtxStatusMc->SetDirectory(0);

  hNeventsMc = new TH1D("hNeventsMc", "N events (Mc vtx) - No Vtx = 0, Vtx OK = 1",
			2, 0, 2);
  hNeventsMc->Sumw2();
  hNeventsMc->SetDirectory(0);

  hNeventsMcTrig = new TH1D("hNeventsMcTrig", "N events (MC vtx+trigger) - No Vtx = 0, Vtx OK = 1",
			   2, 0, 2);
  hNeventsMcTrig->Sumw2();
  hNeventsMcTrig->SetDirectory(0);

  hPtMc = new TH1D("hPtMc", "p_{T} input spectrum (MC); p_{T} [GeV/c]; Counts",
		   nPtBins, ptBins);
  hPtMc->Sumw2();
  hPtMc->SetDirectory(0);

  hPtMcNeg = new TH1D("hPtMcNeg", "p_{T} input spectrum (MC) (q<0); p_{T} [GeV/c]; Counts",
		   nPtBins, ptBins);
  hPtMcNeg->Sumw2();
  hPtMcNeg->SetDirectory(0);

  hPtMcPos = new TH1D("hPtMcPos", "p_{T} input spectrum (MC) (q>0); p_{T} [GeV/c]; Counts",
		   nPtBins, ptBins);
  hPtMcPos->Sumw2();
  hPtMcPos->SetDirectory(0);

  hPtPiMc = new TH1D("hPtPiMc", "p_{T} input spectrum for pi (MC); p_{T} [GeV/c]; Counts",
		     nPtBins, ptBins);
  hPtPiMc->Sumw2();
  hPtPiMc->SetDirectory(0);

  hPtPiMcNeg = new TH1D("hPtPiMcNeg", "p_{T} input spectrum for pi (MC) (q<0); p_{T} [GeV/c]; Counts",
			nPtBins, ptBins);
  hPtPiMcNeg->Sumw2();
  hPtPiMcNeg->SetDirectory(0);
  
  hPtPiMcPos = new TH1D("hPtPiMcPos", "p_{T} input spectrum for pi (MC) (q>0); p_{T} [GeV/c]; Counts",
			nPtBins, ptBins);
  hPtPiMcPos->Sumw2();
  hPtPiMcPos->SetDirectory(0);

  hPtKMc = new TH1D("hPtKMc", "p_{T} input spectrum for k (MC); p_{T} [GeV/c]; Counts",
		     nPtBins, ptBins);
  hPtKMc->Sumw2();
  hPtKMc->SetDirectory(0);

  hPtKMcNeg = new TH1D("hPtKMcNeg", "p_{T} input spectrum for k (MC) (q<0); p_{T} [GeV/c]; Counts",
			nPtBins, ptBins);
  hPtKMcNeg->Sumw2();
  hPtKMcNeg->SetDirectory(0);
  
  hPtKMcPos = new TH1D("hPtKMcPos", "p_{T} input spectrum for k (MC) (q>0); p_{T} [GeV/c]; Counts",
			nPtBins, ptBins);
  hPtKMcPos->Sumw2();
  hPtKMcPos->SetDirectory(0);
  
  hPtPMc = new TH1D("hPtPMc", "p_{T} input spectrum for p (MC); p_{T} [GeV/c]; Counts",
		     nPtBins, ptBins);
  hPtPMc->Sumw2();
  hPtPMc->SetDirectory(0);

  hPtPMcNeg = new TH1D("hPtPMcNeg", "p_{T} input spectrum for p (MC) (q<0); p_{T} [GeV/c]; Counts",
			nPtBins, ptBins);
  hPtPMcNeg->Sumw2();
  hPtPMcNeg->SetDirectory(0);
  
  hPtPMcPos = new TH1D("hPtPMcPos", "p_{T} input spectrum for p (MC) (q>0); p_{T} [GeV/c]; Counts",
			nPtBins, ptBins);
  hPtPMcPos->Sumw2();
  hPtPMcPos->SetDirectory(0);

  hMeanPtMc = new TProfile("hMeanPtMc", "mean p_{T}; p_{T} [GeV/c]; mean p_{T}",
			   nPtBins, ptBins);
  hMeanPtMc->SetDirectory(0);
}

//_________________________________________________________
Bool_t AliHighPtDeDxMc::TrackAcceptedMc()
{
  if(fUseEtaCut && (fTrackEtaMc<fEtaLow || fTrackEtaMc>fEtaHigh))
    return kFALSE;
  
  // only accept hadrons = pi(1), K(2), p(3), other (Sigma+ etc. = 6)
  if(fTrackPidMc < 1 || fTrackPidMc == 4 || fTrackPidMc == 5 || fTrackPidMc > 6)
    return kFALSE;
  
  return kTRUE;
}

//_________________________________________________________
void AliHighPtDeDxMc::FillEventInfo() 
{
  //
  // We require that the MC vtx is withion our vertex range
  // (fEventVtxStatusMc==1). In this way we have to accept both
  // fEventVtxStatusMc==0 and fEventVtxStatusMc==1 to take into account
  // migration effects. 
  //
  // The histogram associated with this is hNeventsMc (together with hNevents
  // in the base class)
  //
  // To correct down the "no vtx" events in the data there are two possibilities
  //
  // 1) We can use the data fraction between accepted vtx and all vtx. Note
  // that since the vtx rec efficiency is not constant vs the vtx position we
  // make a small mistake here, but on the other hand we might be less
  // sensitive to the MC.
  //
  // 2) We can use the MC prediction for estimating the fraction. In this way
  // we might be have no problems with the slightly varying vtx efficiency,
  // but we might be more sensitive to the MC.
  //
  // The histogram associated with this method is hNeventsMcTrig (together with hNevents
  // in the base class)
  //
  // TO DO: We need to test which is better, e.g., by using PHOJET on PYTHIA
  // and vice versa.
  //
  AliHighPtDeDxBase::FillEventInfo();
 
  hVtxStatusMc->Fill(fEventVtxStatus);

  if(fEventVtxStatusMc==1) {
    if(fEventVtxStatus==-1) { // no vtx class
      hNeventsMc->Fill(0.5);      
      if(fEventTrigger)
	hNeventsMcTrig->Fill(0.5);
    } else {                  // vtx class
      hNeventsMc->Fill(1.5);
      if(fEventTrigger)
	hNeventsMcTrig->Fill(1.5);
    }
  }
}


//_________________________________________________________
TH1D* AliHighPtDeDxMc::GetPtSpectrum() 
{
  TH1D* histSpectrum = dynamic_cast<TH1D*>(hPtMc->Clone("hPtSpectrum"));

  const Double_t nEvents = hNeventsMc->GetBinContent(1) 
    + hNeventsMc->GetBinContent(2);
  const Double_t etaRange = fEtaHigh - fEtaLow;
  
  histSpectrum->Scale(1.0/etaRange);
  NormalizePt(histSpectrum);
  histSpectrum->Scale(1.0/nEvents);

  return histSpectrum;
}

//_________________________________________________________
void AliHighPtDeDxMc::NormalizePt(TH1* hist)
{
  const Int_t n = hist->GetNbinsX();

  for(Int_t bin = 1; bin <= n; bin++) {
    
    const Float_t width = hist->GetXaxis()->GetBinWidth(bin);
    hist->SetBinError(bin, hist->GetBinError(bin)/width);
    hist->SetBinContent(bin, hist->GetBinContent(bin)/width);
  }
}


// //_________________________________________________________
// void AliHighPtDeDxMc::FillTrackInfo(Float_t weight) 
// {
  
//   AliHighPtDeDxBase::FillTrackInfo(weight);
  
//   hPt->Fill(fTrackPt, weight);
// }

//_________________________________________________________
void AliHighPtDeDxMc::FillTrackInfoMc(Float_t weight) 
{
  
  hPtMc->Fill(fTrackPtMc, weight);
  if(fTrackChargeMc<0)
    hPtMcNeg->Fill(fTrackPtMc, weight);
  else
    hPtMcPos->Fill(fTrackPtMc, weight);

  hMeanPtMc->Fill(fTrackPtMc, fTrackPtMc);

  switch (fTrackPidMc) {
    
  case 1: // pion
    hPtPiMc->Fill(fTrackPtMc, weight);
    if(fTrackChargeMc<0)
      hPtPiMcNeg->Fill(fTrackPtMc, weight);
    else
      hPtPiMcPos->Fill(fTrackPtMc, weight);
    break;
  case 2: // kaon
    hPtKMc->Fill(fTrackPtMc, weight);
    if(fTrackChargeMc<0)
      hPtKMcNeg->Fill(fTrackPtMc, weight);
    else
      hPtKMcPos->Fill(fTrackPtMc, weight);
    break;
  case 3: // proton
    hPtPMc->Fill(fTrackPtMc, weight);
    if(fTrackChargeMc<0)
      hPtPMcNeg->Fill(fTrackPtMc, weight);
    else
      hPtPMcPos->Fill(fTrackPtMc, weight);
    break;
  default:
    break;
  }
}

TH1D* AliHighPtDeDxMc::GetHistPtMc(Int_t pid, Int_t charge)
{
  switch (pid) {
	
  case 0:
    if(charge==0)
      return hPtMc;
    else if(charge<0)
      return hPtMcNeg;
    else
      return hPtMcPos;
    break;
  case 1:
    if(charge==0)
      return hPtPiMc;
    else if(charge<0)
      return hPtPiMcNeg;
    else
      return hPtPiMcPos;
    break;
  case 2:
    if(charge==0)
      return hPtKMc;
    else if(charge<0)
      return hPtKMcNeg;
    else
      return hPtKMcPos;
    break;
  case 3:
    if(charge==0)
      return hPtPMc;
    else if(charge<0)
      return hPtPMcNeg;
    else
      return hPtPMcPos;
    break;
  default:
    cout << "PID: " << pid << " not found" << endl;
    break;
  }
  return 0;
}
