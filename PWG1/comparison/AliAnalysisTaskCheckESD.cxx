/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//------------------------------
//
// Proof-enabled 
// version 
// of CheckESD.C
//
//------------------------------

#include "TChain.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TParticle.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDMuonTrack.h"
#include "AliESDCaloCluster.h"
#include "AliRun.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"

#include "AliAnalysisTaskCheckESD.h"

ClassImp(AliAnalysisTaskCheckESD)

AliAnalysisTaskCheckESD::AliAnalysisTaskCheckESD():
AliAnalysisTaskSE("AliAnalysisTaskCheckESD"),
  fListOfHistos(0),
  fGen(0),
  fRec(0),
  fResPtInv(0),
  fResPhi(0),
  fResTheta(0),
  fDEdxRight(0),
  fDEdxWrong(0),
  fResTOFRight(0),
  fResTOFWrong(0),
  fEPHOS(0),
  fEEMCAL(0),
  fPtMUON(0),
  fMassK0(0),
  fMassLambda(0),
  fMassLambdaBar(0),
  fMassXi(0),
  fMassOmega(0),
  fScalars(0),
  fArrayHist(0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskCheckESD::AliAnalysisTaskCheckESD(const char* name):
AliAnalysisTaskSE(name),
  fListOfHistos(0),
  fGen(0),
  fRec(0),
  fResPtInv(0),
  fResPhi(0),
  fResTheta(0),
  fDEdxRight(0),
  fDEdxWrong(0),
  fResTOFRight(0),
  fResTOFWrong(0),
  fEPHOS(0),
  fEEMCAL(0),
  fPtMUON(0),
  fMassK0(0),
  fMassLambda(0),
  fMassLambdaBar(0),
  fMassXi(0),
  fMassOmega(0),
  fScalars(0),
  fArrayHist(0)
{
  // Constructor
  AliInfo("Constructor AliAnalysisTaskCheckESD");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}

TH1F * AliAnalysisTaskCheckESD::CreateHisto(const char* name, const char* title,Int_t nBins, 
					    Double_t xMin, Double_t xMax,
					    const char* xLabel, const char* yLabel)
{
  // create a histogram
  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  result->SetOption("E");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  result->SetMarkerStyle(kFullCircle);
  return result;
}

TH1F *AliAnalysisTaskCheckESD::CreateEffHisto(const TH1F* hGen, const TH1F* hRec)
{
  // create an efficiency histogram
  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);
  hEff->GetYaxis()->SetTitle("#epsilon [%]");
  
  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGenEff = hGen->GetBinContent(iBin);
    Double_t nRecEff = hRec->GetBinContent(iBin);
    if (nGenEff > 0) {
      Double_t eff = nRecEff/nGenEff;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt(eff*(1.-eff) / nGenEff);
      if (error == 0) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);			
    }
    else {
      hEff->SetBinContent(iBin, -100.);
      hEff->SetBinError(iBin, 0);
    }
  }
  return hEff;
}


Bool_t AliAnalysisTaskCheckESD::FitHisto(TH1* histo, Double_t& res, Double_t& resError)
{
  // fit a gaussian to a histogram
  static TF1* fitFunc = new TF1("fitFunc", "gaus");
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  Double_t maxFitRange = 2;
  
  if (histo->Integral() > 50) {
    Float_t mean = histo->GetMean();
    Float_t rms = histo->GetRMS();
    fitFunc->SetRange(mean - maxFitRange*rms, mean + maxFitRange*rms);
    fitFunc->SetParameters(mean, rms);
    histo->Fit(fitFunc, "QRI0");
    histo->GetFunction("fitFunc")->ResetBit(1<<9);
    res = TMath::Abs(fitFunc->GetParameter(2));
    resError = TMath::Abs(fitFunc->GetParError(2));
    return kTRUE;
  }
  return kFALSE;
}

void AliAnalysisTaskCheckESD::UserCreateOutputObjects()
{
  // Create histograms
  AliInfo("AliAnalysisTaskCheckESD::UserCreateOutputObjects");
  // Create output container
  fListOfHistos = new TList();
  
  // efficiency and resolution histog
  Int_t nBinsPt = 15;
  Float_t minPt = 0.1;
  Float_t maxPt = 3.1;
  
  fGen = CreateHisto("hGen", "generated tracks", nBinsPt, minPt, maxPt, "p_{t} [GeV/c]", "N");
  fRec = CreateHisto("hRec", "reconstructed tracks", nBinsPt, minPt, maxPt, "p_{t} [GeV/c]", "N");
  fResPtInv = CreateHisto("hResPtInv", "", 100, -10, 10, "(p_{t,rec}^{-1}-p_{t,sim}^{-1}) / p_{t,sim}^{-1} [%]", "N");
  fResPhi = CreateHisto("hResPhi", "", 100, -20, 20, "#phi_{rec}-#phi_{sim} [mrad]", "N");
  fResTheta = CreateHisto("hResTheta", "", 100, -20, 20, "#theta_{rec}-#theta_{sim} [mrad]", "N");
  
  // dE/dx and TOF
  fDEdxRight = new TH2F("hDEdxRight", "", 300, 0, 3, 100, 0, 400);
  fDEdxRight->SetStats(kFALSE);
  fDEdxRight->GetXaxis()->SetTitle("p [GeV/c]");
  fDEdxRight->GetYaxis()->SetTitle("dE/dx_{TPC}");
  fDEdxRight->SetMarkerStyle(kFullCircle);
  fDEdxRight->SetMarkerSize(0.4);
  fDEdxWrong = new TH2F("hDEdxWrong", "", 300, 0, 3, 100, 0, 400);
  fDEdxWrong->SetStats(kFALSE);
  fDEdxWrong->GetXaxis()->SetTitle("p [GeV/c]");
  fDEdxWrong->GetYaxis()->SetTitle("dE/dx_{TPC}");
  fDEdxWrong->SetMarkerStyle(kFullCircle);
  fDEdxWrong->SetMarkerSize(0.4);
  fDEdxWrong->SetMarkerColor(kRed);
  
  fResTOFRight = CreateHisto("hResTOFRight", "", 100, -1000, 1000, "t_{TOF}-t_{track} [ps]", "N");
  fResTOFWrong = CreateHisto("hResTOFWrong", "", 100, -1000, 1000, "t_{TOF}-t_{track} [ps]", "N");
  fResTOFWrong->SetLineColor(kRed);
  
  // calorimeters
  fEPHOS = CreateHisto("hEPHOS", "PHOS", 100, 0, 50, "E [GeV]", "N");
  fEEMCAL = CreateHisto("hEEMCAL", "EMCAL", 100, 0, 50, "E [GeV]", "N");
  
  // muons
  fPtMUON = CreateHisto("hPtMUON", "MUON", 100, 0, 20, "p_{t} [GeV/c]", "N");
  
  // V0s and cascades
  fMassK0 = CreateHisto("hMassK0", "K^{0}", 100, 0.4, 0.6, "M(#pi^{+}#pi^{-}) [GeV/c^{2}]", "N");
  fMassLambda = CreateHisto("hMassLambda", "#Lambda", 100, 1.0, 1.2, "M(p#pi^{-}) [GeV/c^{2}]", "N");
  
  fMassLambdaBar = CreateHisto("hMassLambdaBar", "#bar{#Lambda}", 100, 1.0, 1.2, "M(#bar{p}#pi^{+}) [GeV/c^{2}]", "N");
  fMassXi = CreateHisto("hMassXi", "#Xi", 100, 1.2, 1.5, "M(#Lambda#pi) [GeV/c^{2}]", "N");
  fMassOmega = CreateHisto("hMassOmega", "#Omega", 100, 1.5, 1.8, "M(#LambdaK) [GeV/c^{2}]", "N");
  fScalars = new TH1F("hScalars","Container of scalars",8,0,8);
  fArrayHist = new TH1F("hArrayHist","Container for Array",
			(AliPID::kSPECIES+1)*AliPID::kSPECIES,0,(AliPID::kSPECIES+1)*AliPID::kSPECIES);
  
  fListOfHistos->Add(fGen);
  fListOfHistos->Add(fRec);
  fListOfHistos->Add(fResPtInv);
  fListOfHistos->Add(fResPhi);
  fListOfHistos->Add(fResTheta);
  fListOfHistos->Add(fDEdxRight);
  fListOfHistos->Add(fDEdxWrong);
  fListOfHistos->Add(fResTOFRight);
  fListOfHistos->Add(fResTOFWrong);
  fListOfHistos->Add(fEPHOS);
  fListOfHistos->Add(fEEMCAL);
  fListOfHistos->Add(fPtMUON);
  fListOfHistos->Add(fMassK0);
  fListOfHistos->Add(fMassLambda);
  fListOfHistos->Add(fMassLambdaBar);
  fListOfHistos->Add(fMassXi);
  fListOfHistos->Add(fMassOmega);
  fListOfHistos->Add(fScalars);
  fListOfHistos->Add(fArrayHist);
}

void AliAnalysisTaskCheckESD::UserExec(Option_t */*option*/)
{
  // check the content of the ESD
  Double_t cutPtV0 = 0.3;
  Double_t cutPtCascade = 0.5;
  Float_t minPt = 0.1;
  
  // PID
  Int_t partCode[AliPID::kSPECIES] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  Double_t partFrac[AliPID::kSPECIES] = {0.01, 0.01, 0.85, 0.10, 0.05};
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }	
	
  //Primary vertex needed
  TArrayF vertex(3);  
  mcEvent->GenEventHeader()->PrimaryVertex(vertex);
  
  TObjArray selParticles; //Use TClonesArray?
  TObjArray selV0s;
  TObjArray selCascades;
  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {//Check this loop again
    AliMCParticle* track = mcEvent->GetTrack(iTracks);
    TParticle* particle = track->Particle();
    if (!particle) continue;
    if (particle->Pt() < 0.001) continue;
    if (TMath::Abs(particle->Eta()) > 0.9) continue;
    TVector3 dVertex(particle->Vx() - vertex[0], particle->Vy() - vertex[1], particle->Vz() - vertex[2]);
    if (dVertex.Mag() > 0.0001) continue;
    
    switch (TMath::Abs(particle->GetPdgCode())) {
    case kElectron:
    case kMuonMinus:
    case kPiPlus:
    case kKPlus:
    case kProton:
      {
	if (particle->Pt() > minPt) {
	  selParticles.Add(particle);
	  fScalars->Fill(0);
	  fGen->Fill(particle->Pt());
	}
	break;
      }
    case kK0Short:
    case kLambda0:
      {
	if (particle->Pt() > cutPtV0) {
	  fScalars->Fill(3);
	  selV0s.Add(particle);
	}
	break;
      }
    case kXiMinus:
    case kOmegaMinus:
      {
	if (particle->Pt() > cutPtCascade) {
	  fScalars->Fill(6);
	  selCascades.Add(particle);
	}
	break;
      }
    default: break;
    }		
  }
	
  // ESD information
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  // loop over tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = esd->GetTrack(iTrack);

    // select tracks of selected particles
    Int_t label = TMath::Abs(track->GetLabel());
    if (label > mcEvent->GetNumberOfTracks()) continue;     // background
    AliMCParticle* mctrack = mcEvent->GetTrack(label);
    TParticle* particle = mctrack->Particle();
    if (!selParticles.Contains(particle)) continue;
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;
    if (track->GetConstrainedChi2() > 1e9) continue;
    selParticles.Remove(particle);   // don't count multiple tracks
    
    fScalars->Fill(1);
    fRec->Fill(particle->Pt());
    if (track->GetLabel() < 0) { 
      fScalars->Fill(2);
    }

    // resolutions
    fResPtInv->Fill(100. * (TMath::Abs(track->GetSigned1Pt()) - 1./particle->Pt()) *particle->Pt());
    fResPhi->Fill(1000. * (track->Phi() - particle->Phi()));
    fResTheta->Fill(1000. * (track->Theta() - particle->Theta()));
    
    // PID
    if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) continue;
    Int_t iGen = 5;
    
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (TMath::Abs(particle->GetPdgCode()) == partCode[i]) iGen = i;
    }
    Double_t probability[AliPID::kSPECIES];
    track->GetESDpid(probability);
    Double_t pMax = 0;
    Int_t iRec = 0;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      probability[i] *= partFrac[i];
      if (probability[i] > pMax) {
	pMax = probability[i];
	iRec = i;
      }
			
    }
    fArrayHist->Fill(AliPID::kSPECIES*iGen + iRec);
    if (iGen == iRec) {
      fScalars->Fill(5);
    }

    // dE/dx and TOF
    Double_t time[AliPID::kSPECIES];
    track->GetIntegratedTimes(time);
    if (iGen == iRec) {
      fDEdxRight->Fill(particle->P(), track->GetTPCsignal());
      if ((track->GetStatus() & AliESDtrack::kTOFpid) != 0) {
	fResTOFRight->Fill(track->GetTOFsignal() - time[iRec]);
      }
    }
    else {
      fDEdxWrong->Fill(particle->P(), track->GetTPCsignal());
      if ((track->GetStatus() & AliESDtrack::kTOFpid) != 0) {
	fResTOFWrong->Fill(track->GetTOFsignal() - time[iRec]);
      }
    }
  }

  // loop over muon tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfMuonTracks(); iTrack++) {
    AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
    Double_t ptInv = TMath::Abs(muonTrack->GetInverseBendingMomentum());
    if (ptInv > 0.001) {
      fPtMUON->Fill(1./ptInv);
    }
  }

  // loop over V0s
  for (Int_t iV0 = 0; iV0 < esd->GetNumberOfV0s(); iV0++) {
    AliESDv0* v0 = esd->GetV0(iV0);
    if (v0->GetOnFlyStatus()) continue;
    v0->ChangeMassHypothesis(kK0Short);
    fMassK0->Fill(v0->GetEffMass());
    v0->ChangeMassHypothesis(kLambda0);
    fMassLambda->Fill(v0->GetEffMass());
    v0->ChangeMassHypothesis(kLambda0Bar);
    fMassLambdaBar->Fill(v0->GetEffMass());
    
    Int_t negLabel = TMath::Abs(esd->GetTrack(v0->GetNindex())->GetLabel());
    if (negLabel > mcEvent->GetNumberOfTracks()) continue;     // background
    AliMCParticle* track = mcEvent->GetTrack(negLabel);
    TParticle* particle = track->Particle();
    Int_t negMother = particle->GetMother(0);
    if (negMother < 0) continue;
    Int_t posLabel = TMath::Abs(esd->GetTrack(v0->GetPindex())->GetLabel());
    if (posLabel > mcEvent->GetNumberOfTracks()) continue;     // background
    track = mcEvent->GetTrack(posLabel);
    particle = track->Particle();
    Int_t posMother = particle->GetMother(0);
    if (negMother != posMother) continue;
    track = mcEvent->GetTrack(negMother);
    particle = track->Particle();
    if (!selV0s.Contains(particle)) continue;
    selV0s.Remove(particle);
    fScalars->Fill(4);
  }

  // loop over Cascades
  for (Int_t iCascade = 0; iCascade < esd->GetNumberOfCascades(); iCascade++) {
    AliESDcascade* cascade = esd->GetCascade(iCascade);
    Double_t v0q;
    cascade->ChangeMassHypothesis(v0q,kXiMinus);
    fMassXi->Fill(cascade->GetEffMassXi());
    cascade->ChangeMassHypothesis(v0q,kOmegaMinus);
    fMassOmega->Fill(cascade->GetEffMassXi());
    
    Int_t negLabel = TMath::Abs(esd->GetTrack(cascade->GetNindex())->GetLabel());
    if (negLabel > mcEvent->GetNumberOfTracks()) continue;     // background
    AliMCParticle* track = mcEvent->GetTrack(negLabel);
    TParticle* particle = track->Particle();
    Int_t negMother = particle->GetMother(0);
    if (negMother < 0) continue;
    Int_t posLabel = TMath::Abs(esd->GetTrack(cascade->GetPindex())->GetLabel());
    if (posLabel > mcEvent->GetNumberOfTracks()) continue;     // background
    track = mcEvent->GetTrack(posLabel);
    particle = track->Particle();
    Int_t posMother = particle->GetMother(0);
    if (negMother != posMother) continue;
    track = mcEvent->GetTrack(negMother);
    particle = track->Particle();
    Int_t v0Mother = particle->GetMother(0);
    if (v0Mother < 0) continue;
    Int_t bacLabel = TMath::Abs(esd->GetTrack(cascade->GetBindex())->GetLabel());
    if (bacLabel > mcEvent->GetNumberOfTracks()) continue;     // background
    track = mcEvent->GetTrack(bacLabel);
    particle = track->Particle();
    Int_t bacMother = particle->GetMother(0);
    if (v0Mother != bacMother) continue;
    track = mcEvent->GetTrack(v0Mother);
    particle = track->Particle();
    if (!selCascades.Contains(particle)) continue;
    selCascades.Remove(particle);
    fScalars->Fill(7);
  }
  
  // loop over the clusters
  for (Int_t iCluster=0; iCluster<esd->GetNumberOfCaloClusters(); iCluster++) {
    AliESDCaloCluster * clust = esd->GetCaloCluster(iCluster);
    if (clust->IsPHOS()) fEPHOS->Fill(clust->E());
    if (clust->IsEMCAL()) fEEMCAL->Fill(clust->E());
  }
	
  // Post output data.
  PostData(1, fListOfHistos);
}

void AliAnalysisTaskCheckESD::Terminate(Option_t *)
{
  // check values
  Int_t    checkNGenLow = 1;
  Double_t checkEffLow = 0.5;
  Double_t checkEffSigma = 3;
  Double_t checkFakeHigh = 0.5;
  Double_t checkFakeSigma = 3;
  
  Double_t checkResPtInvHigh = 5;
  Double_t checkResPtInvSigma = 3;
  Double_t checkResPhiHigh = 10;
  Double_t checkResPhiSigma = 3;
  Double_t checkResThetaHigh = 10;
  Double_t checkResThetaSigma = 3;
  
  Double_t checkPIDEffLow = 0.5;
  Double_t checkPIDEffSigma = 3;
  Double_t checkResTOFHigh = 500;
  Double_t checkResTOFSigma = 3;
  
  Double_t checkPHOSNLow = 5;
  Double_t checkPHOSEnergyLow = 0.3;
  Double_t checkPHOSEnergyHigh = 1.0;
  Double_t checkEMCALNLow = 50;
  Double_t checkEMCALEnergyLow = 0.05;
  Double_t checkEMCALEnergyHigh = 1.0;
  
  Double_t checkMUONNLow = 1;
  Double_t checkMUONPtLow = 0.5;
  Double_t checkMUONPtHigh = 10.;
  
  Double_t checkV0EffLow = 0.02;
  Double_t checkV0EffSigma = 3;
  Double_t checkCascadeEffLow = 0.01;
  Double_t checkCascadeEffSigma = 3;
  
  const char* partName[AliPID::kSPECIES+1] = {"electron", "muon", "pion", "kaon", "proton", "other"};
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListOfHistos) {
    Printf("ERROR: fListOfHistos not available");
    return;
  }
	
  fGen = dynamic_cast<TH1F*>(fListOfHistos->At(0));
  fRec = dynamic_cast<TH1F*>(fListOfHistos->At(1));
  fResPtInv = dynamic_cast<TH1F*>(fListOfHistos->At(2));
  fResPhi = dynamic_cast<TH1F*>(fListOfHistos->At(3));
  fResTheta = dynamic_cast<TH1F*>(fListOfHistos->At(4));
  fDEdxRight = dynamic_cast<TH2F*>(fListOfHistos->At(5));
  fDEdxWrong = dynamic_cast<TH2F*>(fListOfHistos->At(6));
  fResTOFRight = dynamic_cast<TH1F*>(fListOfHistos->At(7));
  fResTOFWrong = dynamic_cast<TH1F*>(fListOfHistos->At(8));
  fEPHOS = dynamic_cast<TH1F*>(fListOfHistos->At(9));
  fEEMCAL = dynamic_cast<TH1F*>(fListOfHistos->At(10));
  fPtMUON = dynamic_cast<TH1F*>(fListOfHistos->At(11));
  fMassK0 = dynamic_cast<TH1F*>(fListOfHistos->At(12));
  fMassLambda = dynamic_cast<TH1F*>(fListOfHistos->At(13));
  fMassLambdaBar = dynamic_cast<TH1F*>(fListOfHistos->At(14));
  fMassXi = dynamic_cast<TH1F*>(fListOfHistos->At(15));
  fMassOmega = dynamic_cast<TH1F*>(fListOfHistos->At(16));
  fScalars = dynamic_cast<TH1F*>(fListOfHistos->At(17));
  fArrayHist = dynamic_cast<TH1F*>(fListOfHistos->At(18));
  
  Int_t nGen = Int_t(fScalars->GetBinContent(1));
  Int_t nRec = Int_t(fScalars->GetBinContent(2));
  Int_t nFake = Int_t(fScalars->GetBinContent(3));
  Int_t nGenV0s = Int_t(fScalars->GetBinContent(4));
  Int_t nRecV0s = Int_t(fScalars->GetBinContent(5));
  Int_t nIdentified = Int_t(fScalars->GetBinContent(6));
  Int_t nGenCascades = Int_t(fScalars->GetBinContent(7));
  Int_t nRecCascades = Int_t(fScalars->GetBinContent(8));
  
  Int_t k = 1;
  
  Int_t identified[AliPID::kSPECIES+1][AliPID::kSPECIES];
  for(Int_t i = 0; i < (AliPID::kSPECIES+1); i++)
    for(Int_t j = 0; j < AliPID::kSPECIES; j++) {
      identified[i][j] = Int_t(fArrayHist->GetBinContent(k));
      k++;
    }
  

  // perform checks
  if (nGen < checkNGenLow) {
    Warning("CheckESD", "low number of generated particles: %d", Int_t(nGen));
  }
	
  TH1F* hEff = CreateEffHisto(fGen, fRec);
  
  Info("CheckESD", "%d out of %d tracks reconstructed including %d "
       "fake tracks", nRec, nGen, nFake);
  if (nGen > 0) {
    // efficiency
    Double_t eff = nRec*1./nGen;
    Double_t effError = TMath::Sqrt(eff*(1.-eff) / nGen);
    Double_t fake = nFake*1./nGen;
    Double_t fakeError = TMath::Sqrt(fake*(1.-fake) / nGen);
    Info("CheckESD", "eff = (%.1f +- %.1f) %%  fake = (%.1f +- %.1f) %%",
	 100.*eff, 100.*effError, 100.*fake, 100.*fakeError);
    if (eff < checkEffLow - checkEffSigma*effError) {
      Warning("CheckESD", "low efficiency: (%.1f +- %.1f) %%",100.*eff, 100.*effError);
    }
    if (fake > checkFakeHigh + checkFakeSigma*fakeError) {
      Warning("CheckESD", "high fake: (%.1f +- %.1f) %%",100.*fake, 100.*fakeError);
    }
    // resolutions
    Double_t res, resError;
    if (FitHisto(fResPtInv, res, resError)) {
      Info("CheckESD", "relative inverse pt resolution = (%.1f +- %.1f) %%",res, resError);
      if (res > checkResPtInvHigh + checkResPtInvSigma*resError) {
	Warning("CheckESD", "bad pt resolution: (%.1f +- %.1f) %%",res, resError);
      }
    }
    if (FitHisto(fResPhi, res, resError)) {
      Info("CheckESD", "phi resolution = (%.1f +- %.1f) mrad", res, resError);
      if (res > checkResPhiHigh + checkResPhiSigma*resError) {
	Warning("CheckESD", "bad phi resolution: (%.1f +- %.1f) mrad", 
		res, resError);
      }
    }

    if (FitHisto(fResTheta, res, resError)) {
      Info("CheckESD", "theta resolution = (%.1f +- %.1f) mrad", 
	   res, resError);
      if (res > checkResThetaHigh + checkResThetaSigma*resError) {
	Warning("CheckESD", "bad theta resolution: (%.1f +- %.1f) mrad", 
		res, resError);
      }
    }
    // PID
    if (nRec > 0) {
      Double_t eff = nIdentified*1./nRec;
      Double_t effError = TMath::Sqrt(eff*(1.-eff) / nRec);
      Info("CheckESD", "PID eff = (%.1f +- %.1f) %%", 
	   100.*eff, 100.*effError);
      if (eff < checkPIDEffLow - checkPIDEffSigma*effError) {
	Warning("CheckESD", "low PID efficiency: (%.1f +- %.1f) %%", 
		100.*eff, 100.*effError);
      }
    }

    printf("%9s:", "gen\\rec");
    for (Int_t iRec = 0; iRec < AliPID::kSPECIES; iRec++) {
      printf("%9s", partName[iRec]);
    }
    printf("\n");
    for (Int_t iGen = 0; iGen < AliPID::kSPECIES+1; iGen++) {
      printf("%9s:", partName[iGen]);
      for (Int_t iRec = 0; iRec < AliPID::kSPECIES; iRec++) {
	printf("%9d", identified[iGen][iRec]);
      }
      printf("\n");
    }
    
    if (FitHisto(fResTOFRight, res, resError)) {
      Info("CheckESD", "TOF resolution = (%.1f +- %.1f) ps", res, resError);
      if (res > checkResTOFHigh + checkResTOFSigma*resError) {
	Warning("CheckESD", "bad TOF resolution: (%.1f +- %.1f) ps", 
		res, resError);
      }
    }
    
    // calorimeters
    if (fEPHOS->Integral() < checkPHOSNLow) {
      Warning("CheckESD", "low number of PHOS particles: %d", 
	      Int_t(fEPHOS->Integral()));
    } 	
    else {
      Double_t mean = fEPHOS->GetMean();
      if (mean < checkPHOSEnergyLow) {
	Warning("CheckESD", "low mean PHOS energy: %.1f GeV", mean);
      } else if (mean > checkPHOSEnergyHigh) {
	Warning("CheckESD", "high mean PHOS energy: %.1f GeV", mean);
      }
    }
    
    if (fEEMCAL->Integral() < checkEMCALNLow) {
      Warning("CheckESD", "low number of EMCAL particles: %d", 
	      Int_t(fEEMCAL->Integral()));
    } 
    else {
      Double_t mean = fEEMCAL->GetMean();
      if (mean < checkEMCALEnergyLow) {
	Warning("CheckESD", "low mean EMCAL energy: %.1f GeV", mean);
      }
      else if (mean > checkEMCALEnergyHigh) {
	Warning("CheckESD", "high mean EMCAL energy: %.1f GeV", mean);
      }
    }

    // muons
    if (fPtMUON->Integral() < checkMUONNLow) {
      Warning("CheckESD", "low number of MUON particles: %d", 
	      Int_t(fPtMUON->Integral()));
    }
    else {
      Double_t mean = fPtMUON->GetMean();
      if (mean < checkMUONPtLow) {
	Warning("CheckESD", "low mean MUON pt: %.1f GeV/c", mean);
      }
      else if (mean > checkMUONPtHigh) {
	Warning("CheckESD", "high mean MUON pt: %.1f GeV/c", mean);
      }
    }

    // V0s
    if (nGenV0s > 0) {
      Double_t eff = nRecV0s*1./nGenV0s;
      Double_t effError = TMath::Sqrt(eff*(1.-eff) / nGenV0s);
      if (effError == 0) effError = checkV0EffLow / TMath::Sqrt(1.*nGenV0s);
      Info("CheckESD", "V0 eff = (%.1f +- %.1f) %%", 
	   100.*eff, 100.*effError);
      if (eff < checkV0EffLow - checkV0EffSigma*effError) {
	Warning("CheckESD", "low V0 efficiency: (%.1f +- %.1f) %%", 
		100.*eff, 100.*effError);
      }
    }

    // Cascades
    if (nGenCascades > 0) {
      Double_t eff = nRecCascades*1./nGenCascades;
      Double_t effError = TMath::Sqrt(eff*(1.-eff) / nGenCascades);
      if (effError == 0) effError = checkV0EffLow / 
			   TMath::Sqrt(1.*nGenCascades);
      Info("CheckESD", "Cascade eff = (%.1f +- %.1f) %%", 
	   100.*eff, 100.*effError);
      if (eff < checkCascadeEffLow - checkCascadeEffSigma*effError) {
	Warning("CheckESD", "low Cascade efficiency: (%.1f +- %.1f) %%", 
		100.*eff, 100.*effError);
      }
    }
    
  }

  // draw the histograms if not in batch mode
  if (!gROOT->IsBatch()) {
    new TCanvas;
    hEff->DrawCopy();
    new TCanvas;
    fResPtInv->DrawCopy("E");
    new TCanvas;
    fResPhi->DrawCopy("E");
    new TCanvas;
    fResTheta->DrawCopy("E");
    new TCanvas;
    fDEdxRight->DrawCopy();
    fDEdxWrong->DrawCopy("SAME");
    new TCanvas;
    fResTOFRight->DrawCopy("E");
    fResTOFWrong->DrawCopy("SAME");
    new TCanvas;
    fEPHOS->DrawCopy("E");
    new TCanvas;
    fEEMCAL->DrawCopy("E");
    new TCanvas;
    fPtMUON->DrawCopy("E");
    new TCanvas;
    fMassK0->DrawCopy("E");
    new TCanvas;
    fMassLambda->DrawCopy("E");
    new TCanvas;
    fMassLambdaBar->DrawCopy("E");
    new TCanvas;
    fMassXi->DrawCopy("E");
    new TCanvas;
    fMassOmega->DrawCopy("E");
  }

  // write the output histograms to a file
  TFile* outputFile = TFile::Open("check.root", "recreate");
  if (!outputFile || !outputFile->IsOpen())
    {
      Error("CheckESD", "opening output file check.root failed");
      return;
    }
  hEff->Write();
  fResPtInv->Write();
  fResPhi->Write();
  fResTheta->Write();
  fDEdxRight->Write();
  fDEdxWrong->Write();
  fResTOFRight->Write();
  fResTOFWrong->Write();
  fEPHOS->Write();
  fEEMCAL->Write();
  fPtMUON->Write();
  fMassK0->Write();
  fMassLambda->Write();
  fMassLambdaBar->Write();
  fMassXi->Write();
  fMassOmega->Write();
  outputFile->Close();
  delete outputFile;

  // clean up
  /*  	delete hGen;
  	delete hRec;
  	delete hEff;
  	delete hResPtInv;
  	delete hResPhi;
  	delete hResTheta;
  	delete hDEdxRight;
  	delete hDEdxWrong;
  	delete hResTOFRight;
  	delete hResTOFWrong;
  	delete hEPHOS;
  	delete hEEMCAL;
  	delete hPtMUON;
  	delete hMassK0;
  	delete hMassLambda;
  	delete hMassLambdaBar;
  	delete hMassXi;
  	delete hMassOmega;
	delete hScalars;
	delete hArrayHist; */

  //delete esd;
  // result of check
  Info("CheckESD", "check of ESD was successfull");
}






























































