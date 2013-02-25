// $Id$
//
// Jet spectrum task.
//
// Author: R.Reed, M.Connors

#include "AliAnalysisTaskEmcalJetSpectra.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEmcalJet.h"
#include "AliVCluster.h"
#include "AliRhoParameter.h"
#include "AliEmcalParticle.h"

ClassImp(AliAnalysisTaskEmcalJetSpectra)

//________________________________________________________________________
AliAnalysisTaskEmcalJetSpectra::AliAnalysisTaskEmcalJetSpectra() : 
  AliAnalysisTaskEmcalJet("spectra",kFALSE), 
  fHistRhovsCent(0),
  fHistNjetvsCent(0)
{
  // Default constructor.
  for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistRawJetPtvsTrackPt[i]   = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
    fHistJetPtvsdEP[i]          = 0;
    fHistJetPtvsdEPBias[i]      = 0;
    fHistRhovsEP[i]             = 0;

  }
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetSpectra::AliAnalysisTaskEmcalJetSpectra(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fHistRhovsCent(0),
  fHistNjetvsCent(0)
 { 
   for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistRawJetPtvsTrackPt[i]   = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
    fHistJetPtvsdEP[i]          = 0;
    fHistJetPtvsdEPBias[i]      = 0;
    fHistRhovsEP[i]             = 0;
   }
   SetMakeGeneralHistograms(kTRUE);
 }

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSpectra::UserCreateOutputObjects()
{
  if (! fCreateHisto)
    return;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistRhovsCent             = new TH2F("RhovsCent",              "RhovsCent",             100, 0.0, 100.0, 500, 0, 500);
  fHistNjetvsCent            = new TH2F("NjetvsCent",             "NjetvsCent",            100, 0.0, 100.0, 100, 0, 100);

  TString name;
  TString title;
  for (Int_t i = 0;i<6;++i){
    name = TString(Form("JetPtvsTrackPt_%i",i));
    title = TString(Form("Jet pT vs Leading Track pT cent bin %i",i));
    fHistJetPtvsTrackPt[i] = new TH2F(name,title,1000,-500,500,100,0,100);
    fOutput->Add(fHistJetPtvsTrackPt[i]);
    name = TString(Form("RawJetPtvsTrackPt_%i",i));
    title = TString(Form("Raw Jet pT vs Leading Track pT cent bin %i",i));
    fHistRawJetPtvsTrackPt[i] = new TH2F(name,title,1000,-500,500,100,0,100);
    fOutput->Add(fHistRawJetPtvsTrackPt[i]);
    name = TString(Form("TrackPt_%i",i));
    title = TString(Form("Track pT cent bin %i",i));
    fHistTrackPt[i] = new TH1F(name,title,1000,0,200);
    fOutput->Add(fHistTrackPt[i]);
   
    name = TString(Form("EP0_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0[i] = new TH1F(name,title,100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0[i]);
    name = TString(Form("EP0A_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0A[i] = new TH1F(name,title,100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0A[i]);
    name = TString(Form("EP0C_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0C[i] = new TH1F(name,title,100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0C[i]);
    name = TString(Form("EPAvsC_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEPAvsC[i] = new TH2F(name,title,100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEPAvsC[i]);
    name = TString(Form("JetPtvsdEP_%i",i));
    title = TString(Form("Jet pt vs dEP cent bin %i",i));
    fHistJetPtvsdEP[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsdEP[i]);
    name = TString(Form("JetPtvsdEPBias_%i",i));
    title = TString(Form("Bias Jet pt vs dEP cent bin %i",i));
    fHistJetPtvsdEPBias[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsdEPBias[i]);
    name = TString(Form("JetPtvsEP_%i",i));
    title = TString(Form("Jet pt vs EP cent bin %i",i));
    fHistJetPtvsEP[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsEP[i]);
    name = TString(Form("JetPtvsEPBias_%i",i));
    title = TString(Form("Bias Jet pt vs EP cent bin %i",i));
    fHistJetPtvsEPBias[i] = new TH2F(name,title,1000,-500,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsEPBias[i]);
    name = TString(Form("RhovsEP_%i",i));
    title = TString(Form("Rho vs EP cent bin %i",i));
    fHistRhovsEP[i] = new TH2F(name,title,500,0,500,400,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistRhovsEP[i]);
  }

  
  fOutput->Add(fHistRhovsCent);
  fOutput->Add(fHistNjetvsCent);
  
   PostData(1, fOutput);
}

//________________________________________________________________________

Int_t AliAnalysisTaskEmcalJetSpectra::GetCentBin(Double_t cent) const 
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0;
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}

//________________________________________________________________________

Float_t AliAnalysisTaskEmcalJetSpectra:: RelativePhi(Double_t mphi,Double_t vphi) const
{
  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;//dphi in [-Pi, Pi]                                                                                                    
}


//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSpectra::Run()
{
  Int_t centbin = GetCentBin(fCent);
  //for pp analyses we will just use the first centrality bin
  if (centbin == -1)
    centbin = 0;

  if (!fTracks)
    return kTRUE;
  
  const Int_t nTrack = fTracks->GetEntriesFast();
  for (int i = 0;i<nTrack;i++){
    AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(i));
    if (! track)
      continue;
    fHistTrackPt[centbin]->Fill(track->Pt());
  }

  fHistEP0[centbin]->Fill(fEPV0);
  fHistEP0A[centbin]->Fill(fEPV0A);
  fHistEP0C[centbin]->Fill(fEPV0C);
  fHistEPAvsC[centbin]->Fill(fEPV0A,fEPV0C);
  fRho = GetRhoFromEvent(fRhoName);
  fRhoVal = fRho->GetVal();
  fHistRhovsCent->Fill(fCent,fRhoVal);
  fHistRhovsEP[centbin]->Fill(fRhoVal,fEPV0);
  const Int_t Njets = fJets->GetEntriesFast();

  Int_t NjetAcc = 0;
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
     AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
     if (!jet)
       continue; 
     if (jet->Area()==0)
       continue;
     if (jet->Pt()<0.1)
       continue;
     if (jet->MaxTrackPt()>100)
       continue;
     if (! AcceptJet(jet))
       continue;
     //  jets.push_back(jet);
     NjetAcc++;
     Double_t jetPt = -500;
     jetPt = jet->Pt()-jet->Area()*fRhoVal;    
     fHistJetPtvsTrackPt[centbin]->Fill(jetPt,jet->MaxTrackPt());
     fHistRawJetPtvsTrackPt[centbin]->Fill(jet->Pt(),jet->MaxTrackPt());
     fHistJetPtvsdEP[centbin]->Fill(jetPt,RelativePhi((fEPV0+TMath::Pi()),jet->Phi()));
     fHistJetPtvsEP[centbin]->Fill(jetPt,fEPV0);
     if (jet->MaxTrackPt()>5.0){
       fHistJetPtvsdEPBias[centbin]->Fill(jetPt,RelativePhi((fEPV0+TMath::Pi()),jet->Phi()));
       fHistJetPtvsEPBias[centbin]->Fill(jetPt,fEPV0);
     }
  }
  
  fHistNjetvsCent->Fill(fCent,NjetAcc);
  return kTRUE;
}      





