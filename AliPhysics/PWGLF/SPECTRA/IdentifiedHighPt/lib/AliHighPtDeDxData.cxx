#include "AliHighPtDeDxData.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif

#include <TRandom.h>

using namespace std;

ClassImp(AliHighPtDeDxData);

//
// AliHighPtDeDxData class
//
// This class contains the AliHighPtDeDxData information 
//

//_________________________________________________________
AliHighPtDeDxData::AliHighPtDeDxData():
  AliHighPtDeDxBase(),
  fDeDxPi(0x0),
  fDeDxK(0x0),
  fDeDxP(0x0),
  fDeDxE(0x0),
  fSigmaDeDx(0x0),
  hDeltaPiVsPt(0x0),
  hDeltaPiVsPtNeg(0x0),
  hDeltaPiVsPtPos(0x0),
  hDeltaPiVsPtPiGen(0x0),
  hDeltaPiVsPtPiGenNeg(0x0),
  hDeltaPiVsPtPiGenPos(0x0),
  hDeltaPiVsPtKGen(0x0),
  hDeltaPiVsPtKGenNeg(0x0),
  hDeltaPiVsPtKGenPos(0x0),
  hDeltaPiVsPtPGen(0x0),
  hDeltaPiVsPtPGenNeg(0x0),
  hDeltaPiVsPtPGenPos(0x0),
  hDeltaPiVsPtEGen(0x0),
  hDeltaPiVsPtEGenNeg(0x0),
  hDeltaPiVsPtEGenPos(0x0),
  hDeltaPiVsPtPiMc(0x0),
  hDeltaPiVsPtPiMcNeg(0x0),
  hDeltaPiVsPtPiMcPos(0x0),
  hDeltaPiVsPtKMc(0x0),
  hDeltaPiVsPtKMcNeg(0x0),
  hDeltaPiVsPtKMcPos(0x0),
  hDeltaPiVsPtPMc(0x0),
  hDeltaPiVsPtPMcNeg(0x0),
  hDeltaPiVsPtPMcPos(0x0),
  hPtPi(0x0),
  hPtK(0x0),
  hPtP(0x0),
  hPrimaryVsPidVsPt(0x0)
{
  // default constructor - do not use
}

//_________________________________________________________
AliHighPtDeDxData::AliHighPtDeDxData(const char* name, const char* title):
  AliHighPtDeDxBase(name, title),
  fDeDxPi(0x0),
  fDeDxK(0x0),
  fDeDxP(0x0),
  fDeDxE(0x0),
  fSigmaDeDx(0x0),
  hDeltaPiVsPt(0x0),
  hDeltaPiVsPtNeg(0x0),
  hDeltaPiVsPtPos(0x0),
  hDeltaPiVsPtPiGen(0x0),
  hDeltaPiVsPtPiGenNeg(0x0),
  hDeltaPiVsPtPiGenPos(0x0),
  hDeltaPiVsPtKGen(0x0),
  hDeltaPiVsPtKGenNeg(0x0),
  hDeltaPiVsPtKGenPos(0x0),
  hDeltaPiVsPtPGen(0x0),
  hDeltaPiVsPtPGenNeg(0x0),
  hDeltaPiVsPtPGenPos(0x0),
  hDeltaPiVsPtEGen(0x0),
  hDeltaPiVsPtEGenNeg(0x0),
  hDeltaPiVsPtEGenPos(0x0),
  hDeltaPiVsPtPiMc(0x0),
  hDeltaPiVsPtPiMcNeg(0x0),
  hDeltaPiVsPtPiMcPos(0x0),
  hDeltaPiVsPtKMc(0x0),
  hDeltaPiVsPtKMcNeg(0x0),
  hDeltaPiVsPtKMcPos(0x0),
  hDeltaPiVsPtPMc(0x0),
  hDeltaPiVsPtPMcNeg(0x0),
  hDeltaPiVsPtPMcPos(0x0),
  hPtPi(0x0),
  hPtK(0x0),
  hPtP(0x0),
  hPrimaryVsPidVsPt(0x0)
{
  // named constructor
}

//_________________________________________________________
AliHighPtDeDxData::~AliHighPtDeDxData()
{
  // delete fDeDxPi;
  // delete fDeDxK;
  // delete fDeDxP;
  // delete fSigmaDeDx;
  delete hDeltaPiVsPt;
  delete hDeltaPiVsPtNeg;
  delete hDeltaPiVsPtPos;
  delete hDeltaPiVsPtPiGen;
  delete hDeltaPiVsPtPiGenNeg;
  delete hDeltaPiVsPtPiGenPos;
  delete hDeltaPiVsPtKGen;
  delete hDeltaPiVsPtKGenNeg;
  delete hDeltaPiVsPtKGenPos;
  delete hDeltaPiVsPtPGen;
  delete hDeltaPiVsPtPGenNeg;
  delete hDeltaPiVsPtPGenPos;
  delete hDeltaPiVsPtEGen;
  delete hDeltaPiVsPtEGenNeg;
  delete hDeltaPiVsPtEGenPos;
  delete hDeltaPiVsPtPiMc;
  delete hDeltaPiVsPtPiMcNeg;
  delete hDeltaPiVsPtPiMcPos;
  delete hDeltaPiVsPtKMc;
  delete hDeltaPiVsPtKMcNeg;
  delete hDeltaPiVsPtKMcPos;
  delete hDeltaPiVsPtPMc;
  delete hDeltaPiVsPtPMcNeg;
  delete hDeltaPiVsPtPMcPos;
  delete hPtPi;
  delete hPtK;
  delete hPtP;
  delete hPrimaryVsPidVsPt;
}

//_________________________________________________________
void AliHighPtDeDxData::Init(Int_t nPtBins, Double_t* ptBins)
{
  //
  // Create histograms and functions
  //

  //
  // init base class
  //
  AliHighPtDeDxBase::Init(nPtBins, ptBins);

  const Int_t nDeltaPiBins   = 60;
  const Double_t deltaPiLow = -30;
  const Double_t deltaPiHigh = 30;

  hDeltaPiVsPt = new TH2D("hDeltaPiVsPt", "dE/dx-<dE/dx>_{#pi} vs p_{T}; p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPt->Sumw2();
  hDeltaPiVsPt->SetDirectory(0);
  
  hDeltaPiVsPtNeg = new TH2D("hDeltaPiVsPtNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q < 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtNeg->Sumw2();
  hDeltaPiVsPtNeg->SetDirectory(0);

  hDeltaPiVsPtPos = new TH2D("hDeltaPiVsPtPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q > 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPos->Sumw2();
  hDeltaPiVsPtPos->SetDirectory(0);


  //
  // Generated pions
  //

  hDeltaPiVsPtPiGen = new TH2D("hDeltaPiVsPtPiGen", "dE/dx-<dE/dx>_{#pi} vs p_{T} (#pi gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPiGen->Sumw2();
  hDeltaPiVsPtPiGen->SetDirectory(0);

  hDeltaPiVsPtPiGenNeg = new TH2D("hDeltaPiVsPtPiGenNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q < 0) (#pi gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPiGenNeg->Sumw2();
  hDeltaPiVsPtPiGenNeg->SetDirectory(0);

  hDeltaPiVsPtPiGenPos = new TH2D("hDeltaPiVsPtPiGenPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q > 0) (#pi gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPiGenPos->Sumw2();
  hDeltaPiVsPtPiGenPos->SetDirectory(0);

  //
  // Generated kaons
  //

  hDeltaPiVsPtKGen = new TH2D("hDeltaPiVsPtKGen", "dE/dx-<dE/dx>_{#pi} vs p_{T} (K gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtKGen->Sumw2();
  hDeltaPiVsPtKGen->SetDirectory(0);

  hDeltaPiVsPtKGenNeg = new TH2D("hDeltaPiVsPtKGenNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q < 0) (K gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtKGenNeg->Sumw2();
  hDeltaPiVsPtKGenNeg->SetDirectory(0);

  hDeltaPiVsPtKGenPos = new TH2D("hDeltaPiVsPtKGenPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q > 0) (K gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtKGenPos->Sumw2();
  hDeltaPiVsPtKGenPos->SetDirectory(0);

  //
  // Generated protons
  //

  hDeltaPiVsPtPGen = new TH2D("hDeltaPiVsPtPGen", "dE/dx-<dE/dx>_{#pi} vs p_{T} (p gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPGen->Sumw2();
  hDeltaPiVsPtPGen->SetDirectory(0);

  hDeltaPiVsPtPGenNeg = new TH2D("hDeltaPiVsPtPGenNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q < 0) (p gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPGenNeg->Sumw2();
  hDeltaPiVsPtPGenNeg->SetDirectory(0);

  hDeltaPiVsPtPGenPos = new TH2D("hDeltaPiVsPtPGenPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q > 0) (p gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtPGenPos->Sumw2();
  hDeltaPiVsPtPGenPos->SetDirectory(0);

  //
  // Generated electrons
  //

  hDeltaPiVsPtEGen = new TH2D("hDeltaPiVsPtEGen", "dE/dx-<dE/dx>_{#pi} vs p_{T} (e gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtEGen->Sumw2();
  hDeltaPiVsPtEGen->SetDirectory(0);

  hDeltaPiVsPtEGenNeg = new TH2D("hDeltaPiVsPtEGenNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q < 0) (e gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtEGenNeg->Sumw2();
  hDeltaPiVsPtEGenNeg->SetDirectory(0);

  hDeltaPiVsPtEGenPos = new TH2D("hDeltaPiVsPtEGenPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (q > 0) (e gen); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			  nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
  hDeltaPiVsPtEGenPos->Sumw2();
  hDeltaPiVsPtEGenPos->SetDirectory(0);

  //
  // MC
  //
  if(fIsMc) {
    hDeltaPiVsPtPiMc = new TH2D("hDeltaPiVsPtPiMc", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtPiMc->Sumw2();
    hDeltaPiVsPtPiMc->SetDirectory(0);

    hDeltaPiVsPtPiMcNeg = new TH2D("hDeltaPiVsPtPiMcNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC) (q < 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtPiMcNeg->Sumw2();
    hDeltaPiVsPtPiMcNeg->SetDirectory(0);

    hDeltaPiVsPtPiMcPos = new TH2D("hDeltaPiVsPtPiMcPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC) (q > 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtPiMcPos->Sumw2();
    hDeltaPiVsPtPiMcPos->SetDirectory(0);
    
    hDeltaPiVsPtKMc = new TH2D("hDeltaPiVsPtKMc", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			       nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtKMc->Sumw2();
    hDeltaPiVsPtKMc->SetDirectory(0);

    hDeltaPiVsPtKMcNeg = new TH2D("hDeltaPiVsPtKMcNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC) (q < 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtKMcNeg->Sumw2();
    hDeltaPiVsPtKMcNeg->SetDirectory(0);

    hDeltaPiVsPtKMcPos = new TH2D("hDeltaPiVsPtKMcPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC) (q > 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtKMcPos->Sumw2();
    hDeltaPiVsPtKMcPos->SetDirectory(0);
    
    hDeltaPiVsPtPMc = new TH2D("hDeltaPiVsPtPMc", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
			       nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtPMc->Sumw2();
    hDeltaPiVsPtPMc->SetDirectory(0);

    hDeltaPiVsPtPMcNeg = new TH2D("hDeltaPiVsPtPMcNeg", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC) (q < 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtPMcNeg->Sumw2();
    hDeltaPiVsPtPMcNeg->SetDirectory(0);

    hDeltaPiVsPtPMcPos = new TH2D("hDeltaPiVsPtPMcPos", "dE/dx-<dE/dx>_{#pi} vs p_{T} (MC) (q > 0); p_{T} [GeV/c]; dE/dx - <dE/dx>_{pi}",
				nPtBins, ptBins, nDeltaPiBins, deltaPiLow, deltaPiHigh);
    hDeltaPiVsPtPMcPos->Sumw2();
    hDeltaPiVsPtPMcPos->SetDirectory(0);

    
    const Int_t nPidBins = 7;
    const Double_t pidBinSize = 1;
    Double_t pidBins[nPidBins+1];

    for(Int_t i = 0; i <= nPidBins; i++) {

      pidBins[i] = pidBinSize*i - 0.5;
    }

    const Int_t nPrimaryBins = 2;
    const Double_t primaryBinSize = 1;
    Double_t primaryBins[nPrimaryBins+1];

    for(Int_t i = 0; i <= nPrimaryBins; i++) {
      
      primaryBins[i] = primaryBinSize*i - 0.5;
    }
  
    hPrimaryVsPidVsPt = new TH3D("hPrimaryVsPidVsPt", "primary status vs pid vs Pt; p_{T} [GeV/c]; Pid; Primary status", 
				 nPtBins, ptBins, nPidBins, pidBins, nPrimaryBins, primaryBins);
  }
}

//_________________________________________________________
void AliHighPtDeDxData::FillTrackInfo(Float_t weight) 
{
  AliHighPtDeDxBase::FillTrackInfo(weight);
  
  const Double_t dedxPi  = fDeDxPi->Eval(fTrackP);
  const Double_t sigmaPi = fSigmaDeDx->Eval(dedxPi);
  
  const Double_t dedxK   = fDeDxK->Eval(fTrackP);
  const Double_t sigmaK  = fSigmaDeDx->Eval(dedxK);
  
  const Double_t dedxP   = fDeDxP->Eval(fTrackP);
  const Double_t sigmaP  = fSigmaDeDx->Eval(dedxP);

  const Double_t dedxE   = fDeDxE->Eval(fTrackP);
  const Double_t sigmaE  = fSigmaDeDx->Eval(dedxE);
  
  hDeltaPiVsPt->Fill(fTrackPt, fTrackDeDx-dedxPi);
  if(fTrackCharge<0)
    hDeltaPiVsPtNeg->Fill(fTrackPt, fTrackDeDx-dedxPi);
  else
    hDeltaPiVsPtPos->Fill(fTrackPt, fTrackDeDx-dedxPi);
  
  // Fill MC info
  if(fIsMc) {
    
    hPrimaryVsPidVsPt->Fill(fTrackPt, fTrackPidMc, fTrackPrimaryMc); 
    switch (fTrackPidMc) {
      
    case 1: // pion
      hDeltaPiVsPtPiMc->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      if(fTrackCharge<0)
	hDeltaPiVsPtPiMcNeg->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      else
	hDeltaPiVsPtPiMcPos->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      break;
    case 2: // kaon
      hDeltaPiVsPtKMc->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      if(fTrackCharge<0)
	hDeltaPiVsPtKMcNeg->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      else
	hDeltaPiVsPtKMcPos->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      break;
    case 3: // proton
      hDeltaPiVsPtPMc->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      if(fTrackCharge<0)
	hDeltaPiVsPtPMcNeg->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      else
	hDeltaPiVsPtPMcPos->Fill(fTrackPt, fTrackDeDx-dedxPi, weight);
      break;
    default:
      break;
    }
  }

  
  for(Int_t i = 0; i< 10; i++) {

    const Double_t piShape = gRandom->Gaus(0, sigmaPi);
    const Double_t kShape  = gRandom->Gaus(dedxK - dedxPi, sigmaK);
    const Double_t pShape  = gRandom->Gaus(dedxP - dedxPi, sigmaP);
    const Double_t eShape  = gRandom->Gaus(dedxE - dedxPi, sigmaE);
    hDeltaPiVsPtPiGen->Fill(fTrackPt, piShape);
    hDeltaPiVsPtKGen->Fill(fTrackPt, kShape);
    hDeltaPiVsPtPGen->Fill(fTrackPt, pShape);
    hDeltaPiVsPtEGen->Fill(fTrackPt, eShape);
    if(fTrackCharge<0) {

      hDeltaPiVsPtPiGenNeg->Fill(fTrackPt, piShape);
      hDeltaPiVsPtKGenNeg->Fill(fTrackPt, kShape);
      hDeltaPiVsPtPGenNeg->Fill(fTrackPt, pShape);
      hDeltaPiVsPtEGenNeg->Fill(fTrackPt, eShape);
    } else {
      
      hDeltaPiVsPtPiGenPos->Fill(fTrackPt, piShape);
      hDeltaPiVsPtKGenPos->Fill(fTrackPt, kShape);
      hDeltaPiVsPtPGenPos->Fill(fTrackPt, pShape);
      hDeltaPiVsPtEGenPos->Fill(fTrackPt, eShape);
    }
  }
}

TH2D* AliHighPtDeDxData::GetHistDeltaPiVsPt(Int_t pid, Int_t charge)
{
  switch (pid) {
	
  case 0:
    if(charge==0)
      return hDeltaPiVsPt;
    else if(charge<0)
      return hDeltaPiVsPtNeg;
    else
      return hDeltaPiVsPtPos;
    break;
  case 1:
    if(charge==0)
      return hDeltaPiVsPtPiGen;
    else if(charge<0)
      return hDeltaPiVsPtPiGenNeg;
    else
      return hDeltaPiVsPtPiGenPos;
    break;
  case 2:
    if(charge==0)
      return hDeltaPiVsPtKGen;
    else if(charge<0)
      return hDeltaPiVsPtKGenNeg;
    else
      return hDeltaPiVsPtKGenPos;
    break;
  case 3:
    if(charge==0)
      return hDeltaPiVsPtPGen;
    else if(charge<0)
      return hDeltaPiVsPtPGenNeg;
    else
      return hDeltaPiVsPtPGenPos;
    break;
  case 4:
    if(charge==0)
      return hDeltaPiVsPtEGen;
    else if(charge<0)
      return hDeltaPiVsPtEGenNeg;
    else
      return hDeltaPiVsPtEGenPos;
    break;
  default:
    cout << "PID: " << pid << " not found" << endl;
    break;
  }
  return 0;
}

TH2D* AliHighPtDeDxData::GetHistDeltaPiVsPtMc(Int_t pid, Int_t charge)
{
  switch (pid) {
	
  case 1:
    if(charge==0)
      return hDeltaPiVsPtPiMc;
    else if(charge<0)
      return hDeltaPiVsPtPiMcNeg;
    else
      return hDeltaPiVsPtPiMcPos;
    break;
  case 2:
    if(charge==0)
      return hDeltaPiVsPtKMc;
    else if(charge<0)
      return hDeltaPiVsPtKMcNeg;
    else
      return hDeltaPiVsPtKMcPos;
    break;
  case 3:
    if(charge==0)
      return hDeltaPiVsPtPMc;
    else if(charge<0)
      return hDeltaPiVsPtPMcNeg;
    else
      return hDeltaPiVsPtPMcPos;
    break;
  default:
    cout << "PID: " << pid << " not found" << endl;
    break;
  }
  return 0;
}
