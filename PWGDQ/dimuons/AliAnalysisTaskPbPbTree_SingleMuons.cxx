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
// Analysis task to create a tree containing muon tracks
// R. Arnaldi
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskPbPbTree_SingleMuons_CXX
//#define AliAnalysisTaskPbPbTree_SingleMuons_CXX

#include "TROOT.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TGrid.h"
#include "TProcessID.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TParticle.h"

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
using namespace std;


ClassImp(AliAnalysisTaskPbPbTree_SingleMuons)
//__________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons::AliAnalysisTaskPbPbTree_SingleMuons() :
  AliAnalysisTaskSE(),
  fOutputTree(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fPercentV0M(0x0),
  fAODEvent(0x0)
{
   //
  //Default ctor
  //
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTracksCuts", "StandardMuonTracksCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

}

//__________________________________________________________________________
AliAnalysisTaskPbPbTree_SingleMuons::AliAnalysisTaskPbPbTree_SingleMuons(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputTree(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fPercentV0M(0x0),
  fAODEvent(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskPbPbTree_SingleMuons","Calling Constructor");

  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTracksCuts", "TestStandardMuonTracksCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

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
  fCountCINT7(c.fCountCINT7),
  fCountCMUL7(c.fCountCMUL7),
  fCountCMLL7(c.fCountCMLL7),
  fCountCMSL7(c.fCountCMSL7),
  fCountCMSH7(c.fCountCMSH7),
  fNMuons(c.fNMuons),
  fPercentV0M(c.fPercentV0M),
  fAODEvent(c.fAODEvent),
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
  fMuonTracks->Delete();
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

  fMuonTracks = new TObjArray();
  fMuonTracks->SetOwner();

  fOutputTree->Branch("MuonTracks","TObjArray",&fMuonTracks,256000);
  fOutputTree->Branch("NMuons",&fNMuons,"NMuons/I");
  fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/F");
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
  fNMuons=0;
  fPercentV0M=-1.;
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

  //   read physics selection
  UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t IsPhysSelected = fSelectMask & (AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonSingleLowPt7 |   AliVEvent::kMuonSingleHighPt7  | AliVEvent::kINT7inMUON  | AliVEvent::kINT7);

  // read centrality
  Float_t PercV0M = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) fAODEvent -> FindListObject("MultSelection");
  PercV0M = MultSelection -> GetMultiplicityPercentile("V0M"); //second argument in kFALSE
  fPercentV0M = PercV0M;

  Int_t nummu = 0;
  Int_t ntracks = fAODEvent->GetNumberOfTracks();

  // loop on muons - write only muons surviving cuts
  fMuonTracks->Delete();

  if(IsPhysSelected){
    if(TriggerSelected_CMSL7){
      if(ntracks!=0) {
        for (Int_t i=0;i<ntracks;i++){

          AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
          if(!mu0->IsMuonTrack()) continue;
          if(mu0->Eta() <-4 || mu0->Eta() >-2.5) continue;
          if(mu0->GetRAtAbsorberEnd()<17.6 || mu0->GetRAtAbsorberEnd()>89.5) continue;
          if(mu0->GetMatchTrigger()<=1) continue;
          if(!fMuonTrackCuts->IsSelected(mu0)) continue;


          Int_t pdg=999;
          if(mu0->Charge()==-1) pdg=13;
          else if(mu0->Charge()==1) pdg = -13;
          Int_t status=0;
          Int_t mother1=0;
          Int_t mother2=0;
          Int_t daughter1=0;
          Int_t daughter2=0;
          Double_t px = mu0->Px();
          Double_t py = mu0->Py();
          Double_t pz = mu0->Pz();
          Double_t etot = mu0->E();
          Double_t vx=0;
          Double_t vy=0;
          Double_t vz=0;
          Double_t time=0;

          TParticle *MuonTr = new TParticle(pdg,status,mother1,mother2,daughter1,daughter2,px,py,pz,etot,vx,vy,vz,time);

          fMuonTracks-> AddLast(new TParticle (*MuonTr));
          nummu++;
        }
        fNMuons = nummu;
      }
    }
  }
    //MuonTracks.Print();
   // save only events containing muons surviving all standard cuts
   if(fNMuons>0){
    fOutputTree->Fill();
    PostData(1,fOutputTree);
   }

  // keep all events for trigger summary
  PostData(2,fhNEv);

}
//________________________________________________________________________
void AliAnalysisTaskPbPbTree_SingleMuons::Terminate(Option_t *)
{

 }
