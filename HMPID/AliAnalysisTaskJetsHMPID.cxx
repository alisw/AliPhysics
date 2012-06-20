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

//=================================================================================
// AliAnalysysTaskJetsHMPID - Class performing PID analysis in jets with the HMPID
// A set of histograms is created.
//=================================================================================
//
// By means of AnalysisTrainHMPID.C macro it is possible to use this class
// to perform the analysis on local data, on data on alien using local machine
// and on CAF.

#include <TClonesArray.h>
#include <TRandom2.h>
#include "AliAnalysisTaskJetsHMPID.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODJetEventBackground.h"

ClassImp(AliAnalysisTaskJetsHMPID)

//__________________________________________________________________________
AliAnalysisTaskJetsHMPID::AliAnalysisTaskJetsHMPID() :
  fJetBranch("jets"),
  fBkgBranch(""),
  fJetPtCut(10.),
  fAOD(0x0),
  fHistList(0x0),
  fThetaChJet(0x0),
  fThetaChBkg(0x0),
  fThetaChRndCone(0x0),
  fEvSelDPhi(0x0),
  fJetsPt(0x0),
  fRndConePt(0x0),
  fAwayJetPt(0x0),
  fJetsEtaPhi(0x0),
  fTrksEtaPhiJet(0x0),
  fTrksEtaPhiBkg(0x0),
  fTree(0x0),
  fTrackPt(0),
  fJetPt(0),
  fPionBkg(1.1),
  fKaonBkg(1.1),
  fProtBkg(1.1),
  fPionJet(1.1),
  fKaonJet(1.1),
  fProtJet(1.1)
{
  //
  //Default ctor
  //
}

//__________________________________________________________________________
AliAnalysisTaskJetsHMPID::AliAnalysisTaskJetsHMPID(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fJetBranch("jets"),
  fBkgBranch(""),
  fJetPtCut(10.),
  fAOD(0x0),
  fHistList(0x0),
  fThetaChJet(0x0),
  fThetaChBkg(0x0),
  fThetaChRndCone(0x0),
  fEvSelDPhi(0x0),
  fJetsPt(0x0),
  fRndConePt(0x0),
  fAwayJetPt(0x0),
  fJetsEtaPhi(0x0),
  fTrksEtaPhiJet(0x0),
  fTrksEtaPhiBkg(0x0),
  fTree(0x0),
  fTrackPt(0),
  fJetPt(0),
  fPionBkg(1.1),
  fKaonBkg(1.1),
  fProtBkg(1.1),
  fPionJet(1.1),
  fKaonJet(1.1),
  fProtJet(1.1)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //

  DefineOutput(1,TList::Class());
  DefineOutput(2,TTree::Class());
}

//___________________________________________________________________________
AliAnalysisTaskJetsHMPID& AliAnalysisTaskJetsHMPID::operator=(const AliAnalysisTaskJetsHMPID& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fJetBranch       = c.fJetBranch;
    fBkgBranch       = c.fBkgBranch;
    fJetPtCut        = c.fJetPtCut;
    fAOD             = c.fAOD;
    fHistList        = c.fHistList;
    fThetaChJet      = c.fThetaChJet;
    fThetaChBkg      = c.fThetaChBkg;
    fThetaChRndCone  = c.fThetaChRndCone;
    fEvSelDPhi       = c.fEvSelDPhi;
    fJetsPt          = c.fJetsPt;
    fRndConePt       = c.fRndConePt;
    fAwayJetPt       = c.fAwayJetPt;
    fJetsEtaPhi      = c.fJetsEtaPhi;
    fTrksEtaPhiJet   = c.fTrksEtaPhiJet;
    fTrksEtaPhiBkg   = c.fTrksEtaPhiBkg;
    fTree            = c.fTree;
    fTrackPt         = c.fTrackPt;
    fJetPt           = c.fJetPt;
    fPionBkg         = c.fPionBkg;
    fKaonBkg         = c.fKaonBkg;
    fProtBkg         = c.fProtBkg;
    fPionJet         = c.fPionJet;
    fKaonJet         = c.fKaonJet;
    fProtJet         = c.fProtJet;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskJetsHMPID::AliAnalysisTaskJetsHMPID(const AliAnalysisTaskJetsHMPID& c) :
  AliAnalysisTaskSE(c),
  fJetBranch(c.fJetBranch),
  fBkgBranch(c.fBkgBranch),
  fJetPtCut(c.fJetPtCut),
  fAOD(c.fAOD),
  fHistList(c.fHistList),
  fThetaChJet(c.fThetaChJet),
  fThetaChBkg(c.fThetaChBkg),
  fThetaChRndCone(c.fThetaChRndCone),
  fEvSelDPhi(c.fEvSelDPhi),
  fJetsPt(c.fJetsPt),
  fRndConePt(c.fRndConePt),
  fAwayJetPt(c.fAwayJetPt),
  fJetsEtaPhi(c.fJetsEtaPhi),
  fTrksEtaPhiJet(c.fTrksEtaPhiJet),
  fTrksEtaPhiBkg(c.fTrksEtaPhiBkg),
  fTree(c.fTree),
  fTrackPt(c.fTrackPt),
  fJetPt(c.fJetPt),
  fPionBkg(c.fPionBkg),
  fKaonBkg(c.fKaonBkg),
  fProtBkg(c.fProtBkg),
  fPionJet(c.fPionJet),
  fKaonJet(c.fKaonJet),
  fProtJet(c.fProtJet)
{
  //
  // Copy Constructor
  //
}
 
//___________________________________________________________________________
AliAnalysisTaskJetsHMPID::~AliAnalysisTaskJetsHMPID()
{
  //
  //destructor
  //
  Info("~AliAnalysisTaskJetsHMPID","Calling Destructor");
  if (fHistList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHistList;
}

//__________________________________________________________________________
void AliAnalysisTaskJetsHMPID::ConnectInputData(Option_t *)
{
  if (fDebug) AliInfo("ConnectInputData() ");

  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

  if(handler && handler->InheritsFrom("AliAODInputHandler"))
  { // input AOD
    fAOD = ((AliAODInputHandler*)handler)->GetEvent();
    if (fDebug) AliInfo("Tracks and Jets from AliAODInputHandler");
  }
  else
  {  //output AOD
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if (handler && handler->InheritsFrom("AliAODHandler"))
    {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      if (fDebug > 1) AliInfo("Tracks and Jets from AliAODHandler");
    }
    else
    {  // no AOD
      AliWarning("No AOD event in the input or in the output");
    }
  }

}

//___________________________________________________________________________
void AliAnalysisTaskJetsHMPID::UserCreateOutputObjects()
{
  if (!fHistList) fHistList = new TList();
  fHistList->SetOwner();

  fThetaChJet = new TH2F("ThetaChJet","Theta Cherenkov distribution in the jets",1000,0.,10.,400,0.,0.8);
  fHistList->Add(fThetaChJet);

  fThetaChBkg = new TH2F("ThetaChBkg","Theta Cherenkov distribution out of the jets",1000,0.,10.,400,0.,0.8);
  fHistList->Add(fThetaChBkg);

  fThetaChRndCone = new TH2F("ThetaChRndCone","Theta Cherenkov distribution in random cones",1000,0.,10.,400,0.,0.8);
  fHistList->Add(fThetaChRndCone);

  fEvSelDPhi = new TH1F("PhiJets","Jets phi in different selections",600,0.,6.);
  fHistList->Add(fEvSelDPhi);

  fJetsPt = new TH1F("JetsPt","Jets spectrum in the HMPID",200,0.,200.);
  fHistList->Add(fJetsPt);

  fRndConePt = new TH1F("RndConePt","Pt of random cones in the HMPID",200,0.,200.);
  fHistList->Add(fRndConePt);

  fAwayJetPt = new TH1F("AwayJetPt","Pt of jets opposite to the HMPID",200,0.,200.);
  fHistList->Add(fAwayJetPt);

  fJetsEtaPhi = new TH2F("JetsEtaPhi","Eta-Phi of jets in the HMPID",180,-0.9,0.9,200,-0.25*TMath::Pi(),0.75*TMath::Pi());
  fHistList->Add(fJetsEtaPhi);

  fTrksEtaPhiJet = new TH2F("TrksEtaPhiJet","Eta-Phi of jet tracks in the HMPID",180,-0.9,0.9,200,-0.25*TMath::Pi(),0.75*TMath::Pi());
  fHistList->Add(fTrksEtaPhiJet);

  fTrksEtaPhiBkg = new TH2F("TrksEtaPhiBkg","Eta-Phi of bkg tracks in the HMPID",180,-0.9,0.9,200,-0.25*TMath::Pi(),0.75*TMath::Pi());
  fHistList->Add(fTrksEtaPhiBkg);

  fTree = new TTree("Tree","Tree with data");
  fTree->Branch("Pt",&fTrackPt);
  fTree->Branch("PtJet",&fJetPt);
  fTree->Branch("PionBkg",&fPionBkg);
  fTree->Branch("KaonBkg",&fKaonBkg);
  fTree->Branch("ProtBkg",&fProtBkg);
  fTree->Branch("PionJet",&fPionJet);
  fTree->Branch("KaonJet",&fKaonJet);
  fTree->Branch("ProtJet",&fProtJet);

  PostData(1,fHistList);
  PostData(2,fTree);
}

//___________________________________________________________________________
void AliAnalysisTaskJetsHMPID::UserExec(Option_t *)
{
  if (!fAOD){
    PostData(1,fHistList);
    PostData(2,fTree);
    return;
  }

  const Float_t pi = TMath::Pi(), jetR = 0.4;
  const Float_t hmpPhiMid = 30*TMath::DegToRad();
  const Float_t hmpPhiWid = 25*TMath::DegToRad();

  TClonesArray* jets = (TClonesArray*) fAOD->FindListObject(fJetBranch.Data());
  if (fDebug) printf("Jet Branch: %s\n",fJetBranch.Data());
  if (fDebug) printf("Bkg Branch: %s\n",fBkgBranch.Data());

  if (!jets) AliError("Jet branch not found");

  Int_t nj = jets->GetEntriesFast();
  Float_t *jetPt = new Float_t[nj];
  Int_t *iPerm = new Int_t[nj];
//  Correct jets Pt
  static AliAODJetEventBackground* externalBackground = 0;
  if(fBkgBranch.Length()) externalBackground = (AliAODJetEventBackground*)(fAOD->FindListObject(fBkgBranch.Data()));
  Float_t pTbkg = 0;
  AliAODJet *jet = 0;
  for(Int_t iJet=0; iJet<nj; iJet++){
    jet = (AliAODJet*)jets->At(iJet);
    jetPt[iJet] = jet->Pt();
    pTbkg = (externalBackground) ? externalBackground->GetBackground(2)*jet->EffectiveAreaCharged() : 0;
    jetPt[iJet] -= pTbkg;
  }
  TMath::Sort(nj,jetPt,iPerm,kTRUE);

//  Event characterisation:
//  iEvSel = 0: no jets in the HMPID, consider it as background
//  iEvSel = 1: jet in the HMPID, use it for PID in jets
//  iEvSel = 2: HMPID is in the away side of a leading jet, but no jet found around
//              the detector; dismiss the event as it might alter the analysis

  Int_t iEvSel = -1, iJetSel = -1;
  Float_t phiTrig = 0;
  if (nj && jetPt[iPerm[0]] > fJetPtCut){
    phiTrig = ((AliAODJet*)jets->At(iPerm[0]))->Phi();
    if (phiTrig < pi+hmpPhiMid) phiTrig-=hmpPhiMid;
    else if (phiTrig > pi+hmpPhiMid) phiTrig-=(2*pi+hmpPhiMid);
    if (TMath::Abs(phiTrig) < hmpPhiWid+jetR) iEvSel = 1, iJetSel = 0;
    else if (TMath::Abs(phiTrig) > pi/3 && TMath::Abs(phiTrig) < 2*pi/3) iEvSel = 0;
    else if (TMath::Abs(phiTrig) > 2*pi/3) iEvSel = 2, iJetSel = 0;
  }
  else iEvSel = 0;

  for (Int_t iJet=1; iJet<jets->GetEntriesFast(); iJet++){
    if (jetPt[iPerm[iJet]] < fJetPtCut || iEvSel == 1) break;
    jet = (AliAODJet*)jets->At(iPerm[iJet]);
    phiTrig = jet->Phi();
    if (phiTrig < pi+hmpPhiMid) phiTrig-=hmpPhiMid;
    else if (phiTrig > pi+hmpPhiMid) phiTrig-=(2*pi+hmpPhiMid);
    if (TMath::Abs(phiTrig) < hmpPhiWid+jetR) iEvSel = 1, iJetSel = iJet;
  }

  TRandom2 rndm(0);
  Float_t etaRnd = rndm.Uniform(-0.8,0.8);
  Float_t phiRnd = rndm.Uniform(hmpPhiMid-hmpPhiWid-jetR,hmpPhiMid+hmpPhiWid+jetR);
  AliAODJet *rndJet = new AliAODJet();
  rndJet->SetPtEtaPhiM(10.,etaRnd,phiRnd,0.);
  rndJet->SetEffArea(pi*jetR*jetR,pi*jetR*jetR);
  Float_t ptSum=0., phi=0;
  TRefArray *jetTracks = 0;
  Double_t probs[5]={0, 0, 0, 0, 0};
  switch(iEvSel){
  case 0:
    for (Int_t iTr = 0; iTr < fAOD->GetNumberOfTracks(); iTr++){
      AliAODTrack *tr=(AliAODTrack*)fAOD->GetTrack(iTr);
      if (rndJet->DeltaR(tr) < jetR) ptSum+=tr->Pt();
      AliAODPid *pid = tr->GetDetPid();
      if (!pid) continue;
      if (pid->GetHMPIDsignal() > 0.){
        phi = tr->Phi()>pi ? tr->Phi()-2*pi : tr->Phi();
        fThetaChBkg->Fill(tr->Pt(),pid->GetHMPIDsignal());
        if (rndJet->DeltaR(tr) < jetR) fThetaChRndCone->Fill(tr->Pt(),pid->GetHMPIDsignal());
        fTrksEtaPhiBkg->Fill(tr->Eta(),phi);
        pid->GetHMPIDprobs(probs);
        Short_t charge = tr->Charge();
        fTrackPt = tr->Pt();
        fJetPt   = 0;
        fPionBkg = charge*(probs[0] + probs[1] + probs[2]);
        fKaonBkg = charge*probs[3];
        fProtBkg = charge*probs[4];
        fPionJet = charge*1.2;
        fKaonJet = charge*1.2;
        fProtJet = charge*1.2;
        fTree->Fill();
      }
    }
    pTbkg = (externalBackground) ? externalBackground->GetBackground(2)*jet->EffectiveAreaCharged() : 0;
    fRndConePt->Fill(ptSum-pTbkg);
    rndJet->Delete();
    break;
  case 1:
    jet = (AliAODJet*)jets->At(iPerm[iJetSel]);
    fJetsPt->Fill(jetPt[iPerm[iJetSel]]);
    phi = jet->Phi()>pi ? jet->Phi()-2*pi : jet->Phi();
    fJetsEtaPhi->Fill(jet->Eta(),phi);
    jetTracks = (TRefArray*)jet->GetRefTracks();
    for (Int_t iTr=0; iTr<jetTracks->GetEntriesFast(); iTr++){
      AliAODTrack *tr = (AliAODTrack*)jetTracks->At(iTr);
      AliAODPid *pid = tr->GetDetPid();
      if (!pid) continue;
      if (pid->GetHMPIDsignal() > 0.){
        phi = tr->Phi()>pi ? tr->Phi()-2*pi : tr->Phi();
        fThetaChJet->Fill(tr->Pt(),pid->GetHMPIDsignal());
        fTrksEtaPhiJet->Fill(tr->Eta(),phi);
        pid->GetHMPIDprobs(probs);
        Short_t charge = tr->Charge();
        fTrackPt = tr->Pt();
        fJetPt   = jetPt[iPerm[iJetSel]];
        fPionBkg = charge*1.2;
        fKaonBkg = charge*1.2;
        fProtBkg = charge*1.2;
        fPionJet = charge*(probs[0] + probs[1] + probs[2]);
        fKaonJet = charge*probs[3];
        fProtJet = charge*probs[4];
        fTree->Fill();
      }
    }
    break;
  case 2:
    fAwayJetPt->Fill(jetPt[iPerm[iJetSel]]);
    for (Int_t iTr = 0; iTr < fAOD->GetNumberOfTracks(); iTr++){
      AliAODTrack *tr=(AliAODTrack*)fAOD->GetTrack(iTr);
      AliAODPid *pid = tr->GetDetPid();
      if (!pid) continue;
      if (pid->GetHMPIDsignal() > 0.){
        fThetaChBkg->Fill(tr->Pt(),pid->GetHMPIDsignal());
        pid->GetHMPIDprobs(probs);
        Short_t charge = tr->Charge();
        fTrackPt = tr->Pt();
        fJetPt   = jetPt[iPerm[iJetSel]];
        fPionBkg = charge*(2 + probs[0] + probs[1] + probs[2]);
        fKaonBkg = charge*(2 + probs[3]);
        fProtBkg = charge*(2 + probs[4]);
        fPionJet = charge*1.2;
        fKaonJet = charge*1.2;
        fProtJet = charge*1.2;
        fTree->Fill();
      }
    }
    break;
  }

  if (iEvSel > -1){
    if (iJetSel > -1){
      jet = (AliAODJet*)jets->At(iPerm[iJetSel]);
      fEvSelDPhi->Fill(jet->Phi()/pi+2*iEvSel);
    } else {
      for (Int_t iJet=0; iJet<jets->GetEntriesFast(); iJet++){
        if (jetPt[iPerm[iJet]] < fJetPtCut) break;
        jet = (AliAODJet*)jets->At(iPerm[iJet]);
        fEvSelDPhi->Fill(jet->Phi()/pi+2*iEvSel);
      }
    }
  }

  delete [] jetPt;
  delete [] iPerm;

  PostData(1,fHistList);
  PostData(2,fTree);
  return;
}

//___________________________________________________________________________
void AliAnalysisTaskJetsHMPID::Terminate(Option_t*)
{
  AliAnalysisTaskSE::Terminate();

}
