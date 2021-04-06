/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliAnalysisTaskPbPbTree_SingleMuons.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to compute single muon kinematic distributions
// The output is a tree
// R. Arnaldi
//
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskPbPbTree_SingleMuons_CXX
//#define AliAnalysisTaskPbPbTree_SingleMuons_CXX

#include "TROOT.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TGrid.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"
#include "AliMultSelection.h"

#include "AliAnalysisTaskPbPbTree_SingleMuons.h"
#include "AliTriggerAnalysis.h"


ClassImp(AliAnalysisTaskPbPbTree_SingleMuons)
//__________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons::AliAnalysisTaskPbPbTree_SingleMuons() :
  AliAnalysisTaskSE(),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fPercentV0M(0x0),
  fIsPhysSelected(0x0),
  fAODEvent(0x0),
  fNTracks(0x0)
{
  //
  //Default ctor
  //
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999;
    fPy[i]=999;
    fPz[i]=999;
    fY[i]=999.;
    fEta[i]=999.;
    fPhi[i]=999.;
    fMatchTrig[i]=999.;
    fTrackChi2[i]=999.;
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999;
  }
}

//__________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons::AliAnalysisTaskPbPbTree_SingleMuons(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fPercentV0M(0x0),
  fIsPhysSelected(0x0),
  fAODEvent(0x0),
  fNTracks(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskPbPbTree_SingleMuons","Calling Constructor");

  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "TestStandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999;
    fPy[i]=999;
    fPz[i]=999;
    fY[i]=999.;
    fEta[i]=999.;
    fPhi[i]=999.;
    fMatchTrig[i]=999.;
    fTrackChi2[i]=999.;
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999;
  }

  DefineOutput(1,TTree::Class());
  DefineOutput(2,TH1D::Class());
}

//___________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons& AliAnalysisTaskPbPbTree_SingleMuons::operator=(const AliAnalysisTaskPbPbTree_SingleMuons& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons::AliAnalysisTaskPbPbTree_SingleMuons(const AliAnalysisTaskPbPbTree_SingleMuons& c) :
  AliAnalysisTaskSE(c),
  fOutputTree(c.fOutputTree),
  fNevt(c.fNevt),
  fBeamEnergy(c.fBeamEnergy),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fCountTotEv(c.fCountTotEv),
  fCountTrigger(c.fCountTrigger),
  fCountCINT7(c.fCountCINT7),
  fCountCMUL7(c.fCountCMUL7),
  fCountCMLL7(c.fCountCMLL7),
  fCountCMSL7(c.fCountCMSL7),
  fCountCMSH7(c.fCountCMSH7),
  fNMuons(c.fNMuons),
  fNTracklets(c.fNTracklets),
  fNContributors(c.fNContributors),
  fPercentV0M(c.fPercentV0M),
  fIsPhysSelected(c.fIsPhysSelected),
  fAODEvent(c.fAODEvent),
  fNTracks(c.fNTracks),
  fMuonTrackCuts(c.fMuonTrackCuts)
 
 {
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons::~AliAnalysisTaskPbPbTree_SingleMuons() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskPbPbTree_SingleMuons","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskPbPbTree_SingleMuons::NotifyRun()
{
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetRun((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));
}

//___________________________________________________________________________
void AliAnalysisTaskPbPbTree_SingleMuons::UserCreateOutputObjects(){

  if (fOutputTree) return;

  OpenFile(1,"RECREATE");
  fOutputTree = new TTree("PbPbTree","Data Tree");

  fOutputTree->Branch("FiredTriggerClasses",fTrigClass,"FiredTriggerClasses/C");

  fOutputTree->Branch("NMuons",&fNMuons,"NMuons/I");
  fOutputTree->Branch("Vertex",fVertex,"Vertex[3]/D");
  fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/F");
  fOutputTree->Branch("NTracks",&fNTracks,"NTracks/I");

  fOutputTree->Branch("Pt",fPt,"Pt[NMuons]/D");
  fOutputTree->Branch("E",fE,"E[NMuons]/D");
  fOutputTree->Branch("Px",fPx,"Px[NMuons]/D");
  fOutputTree->Branch("Py",fPy,"Py[NMuons]/D");
  fOutputTree->Branch("Pz",fPz,"Pz[NMuons]/D");
  fOutputTree->Branch("Y",fY,"Y[NMuons]/D");
  fOutputTree->Branch("Eta",fEta,"Eta[NMuons]/D");
  fOutputTree->Branch("Phi",fPhi,"Phi[NMuons]/D");
  fOutputTree->Branch("MatchTrig",fMatchTrig,"MatchTrig[NMuons]/I");
  fOutputTree->Branch("TrackChi2",fTrackChi2,"TrackChi2[NMuons]/D");
  fOutputTree->Branch("MatchTrigChi2",fMatchTrigChi2,"MatchTrigChi2[NMuons]/D");
  fOutputTree->Branch("Charge",fCharge,"Charge[NMuons]/I");
  fOutputTree->Branch("RAtAbsEnd",fRAtAbsEnd,"RAtAbsEnd[NMuons]/D");
  fOutputTree->Branch("pDCA",fpDCA,"pDCA[NMuons]/I");

  fOutputTree->Branch("IsPhysSelected",&fIsPhysSelected,"IsPhysSelected/O");
  fOutputTree->ls();

  PostData(1,fOutputTree);

  fhNEv = new TH1D("fhNEv","hNEv",7,0.,7.);
  TString namelabel1[7]={"TotEv","CINT7","CMUL7","CMLL7","CMSL7","CMSH7","CINT7_CENT"};
  for(int k=0;k<7;k++) fhNEv->GetXaxis()->SetBinLabel(k+1,namelabel1[k]);

  PostData(2,fhNEv);

}

//_________________________________________________
void AliAnalysisTaskPbPbTree_SingleMuons::UserExec(Option_t *)
{

  fNTracks=0;
  fNMuons=0;
  fNTracklets=-1;
  fNContributors=-1;
  fPercentV0M=-1.;
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999;
    fPy[i]=999;
    fPz[i]=999;
    fY[i]=999.;
    fEta[i]=999.;
    fPhi[i]=999.;
    fMatchTrig[i]=999.;
    fTrackChi2[i]=999.;
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999.;
  }

//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }

   AliAODHeader *aodheader=dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
   TString firedtrigger = aodheader->GetFiredTriggerClasses();
   sprintf(fTrigClass,"%s",firedtrigger.Data());

  //   to apply physics selection
    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    fIsPhysSelected = fSelectMask & (AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonSingleLowPt7 |   AliVEvent::kMuonSingleHighPt7  | AliVEvent::kINT7inMUON  | AliVEvent::kINT7);

   Bool_t TriggerSelected=kFALSE;
   Bool_t TriggerSelected_CINT7=kFALSE;
   Bool_t TriggerSelected_CMUL7=kFALSE;
   Bool_t TriggerSelected_CMLL7=kFALSE;
   Bool_t TriggerSelected_CMSL7=kFALSE;
   Bool_t TriggerSelected_CMSH7=kFALSE;
   Bool_t TriggerSelected_CINT7_CENT=kFALSE;

  if(firedtrigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
  else TriggerSelected = kFALSE;
  if(firedtrigger.Contains("CINT7-B-NOPF-MUFAST")) TriggerSelected_CINT7 = kTRUE; 										 
  else TriggerSelected_CINT7 = kFALSE;
  if(firedtrigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected_CMUL7 = kTRUE;
  else TriggerSelected_CMUL7 = kFALSE;
  if(firedtrigger.Contains("CMLL7-B-NOPF-MUFAST")) TriggerSelected_CMLL7 = kTRUE;
  else TriggerSelected_CMLL7 = kFALSE;
  if(firedtrigger.Contains("CMSL7-B-NOPF-MUFAST")) TriggerSelected_CMSL7 = kTRUE;
  else TriggerSelected_CMSL7 = kFALSE;
  if(firedtrigger.Contains("CMSH7-B-NOPF-MUFAST")) TriggerSelected_CMSH7 = kTRUE;
  else TriggerSelected_CMSH7 = kFALSE;
  if(firedtrigger.Contains("CINT7-B-NOPF-CENT")) TriggerSelected_CINT7_CENT = kTRUE;
 
  Double_t DeltaCh=0.;
  fhNEv->Fill(0.+DeltaCh);
  if (TriggerSelected_CINT7) fhNEv->Fill(1.+DeltaCh);
  if (TriggerSelected_CMUL7) fhNEv->Fill(2.+DeltaCh);
  if (TriggerSelected_CMLL7) fhNEv->Fill(3.+DeltaCh);
  if (TriggerSelected_CMSL7) fhNEv->Fill(4.+DeltaCh);
  if (TriggerSelected_CMSH7) fhNEv->Fill(5.+DeltaCh);
  if (TriggerSelected_CINT7_CENT) fhNEv->Fill(6.+DeltaCh);

  // centrality determination
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudiesRun2
  Float_t PercV0M = 300;
 
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) fAODEvent -> FindListObject("MultSelection");
  if( !MultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
    PercV0M = MultSelection -> GetMultiplicityPercentile("V0M"); //second argument in kFALSE
  }
  fPercentV0M = PercV0M;
 
  AliAODVertex *PrimVertex =  fAODEvent->GetPrimaryVertex();
  fVertex[0]=PrimVertex->GetX();
  fVertex[1]=PrimVertex->GetY();
  fVertex[2]=PrimVertex->GetZ();

   Int_t nummu = 0;
   Int_t ntracks = fAODEvent->GetNumberOfTracks();
   if(ntracks!=0) {
   fNTracks = ntracks;

     for (Int_t i=0;i<ntracks;i++){
       AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
       if(!mu0->IsMuonTrack()) continue;
       //if(mu0->Eta() <-4 || mu0->Eta() >-2.5) continue;
       //if(mu0->GetRAtAbsorberEnd()<17.6 || mu0->GetRAtAbsorberEnd()>89.5) continue;
       //if(mu0->GetMatchTrigger()>1){ 
         fCharge[nummu] = mu0->Charge();
         fPt[nummu] = mu0->Pt();
         fPx[nummu] = mu0->Px();
         fPy[nummu] = mu0->Py();
         fPz[nummu] = mu0->Pz();
         fPt[nummu] = mu0->Pt();
         fY[nummu]  = mu0->Y();
         fEta[nummu]= mu0->Eta();
         fPhi[nummu]= mu0->Phi();
         fE[nummu] = mu0->E();
         fMatchTrig[nummu]   = mu0->GetMatchTrigger();
         fMatchTrigChi2[nummu]= mu0->GetChi2MatchTrigger();
         fRAtAbsEnd[nummu]=mu0->GetRAtAbsorberEnd();
         if(fMuonTrackCuts -> IsSelected(mu0)) fpDCA[nummu] = 1;
         nummu++;
       }	 
     //}
     fNMuons = nummu;
   }
   
   if(fNMuons>0){
    fOutputTree->Fill();
    PostData(1,fOutputTree);
  }
  // kepp all events for trigger summary
  PostData(2,fhNEv);

}
//________________________________________________________________________
void AliAnalysisTaskPbPbTree_SingleMuons::Terminate(Option_t *)
{

 }
