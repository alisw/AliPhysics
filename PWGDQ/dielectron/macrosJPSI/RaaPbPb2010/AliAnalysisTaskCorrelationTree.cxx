/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//    Analysis task for creating a reduced data tree                     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
using namespace std;

#include <TChain.h>
#include <TH1D.h>
#include <TFile.h>

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include <AliTriggerAnalysis.h>
#include <AliESDtrackCuts.h>
#include <AliVZDC.h>
#include <AliESDv0.h>
#include <AliESDv0Cuts.h>
#include <AliVCluster.h>
#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowBayesianPID.h"

#include "AliCorrelationReducedEvent.h"
#include "AliAnalysisTaskCorrelationTree.h"

ClassImp(AliAnalysisTaskCorrelationTree)


//_________________________________________________________________________________
AliAnalysisTaskCorrelationTree::AliAnalysisTaskCorrelationTree() :
  AliAnalysisTaskSE(),
  fListDielectron(),
  fListHistos(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fEventFilter(0x0),
  fTrackFilter(0x0),
  fFlowTrackFilter(0x0),
  fK0sCuts(0x0),
  fLambdaCuts(0x0),
  fK0sPionCuts(0x0),
  fLambdaProtonCuts(0x0),
  fLambdaPionCuts(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fV0Histos(0x0),
  fTreeFile(0x0),
  fTree(0x0),
  fFriendTreeFile(0x0),
  fFriendTree(0x0),
  fReducedEvent(0x0),
  fReducedEventFriend(0x0),
  fFlowTrackCuts(0x0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskCorrelationTree::AliAnalysisTaskCorrelationTree(const char *name) :
  AliAnalysisTaskSE(name),
  fListDielectron(),
  fListHistos(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fEventFilter(0x0),
  fTrackFilter(0x0),
  fFlowTrackFilter(0x0),
  fK0sCuts(0x0),
  fLambdaCuts(0x0),
  fK0sPionCuts(0x0),
  fLambdaProtonCuts(0x0),
  fLambdaPionCuts(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fV0Histos(0x0),
  fTreeFile(0x0),
  fTree(0x0),
  fFriendTreeFile(0x0),
  fFriendTree(0x0),
  fReducedEvent(0x0),
  fReducedEventFriend(0x0),
  fFlowTrackCuts(0x0)
{
  //
  // Constructor
  //
  fK0sMassRange[0] = 0.4; fK0sMassRange[1] = 0.6;
  fLambdaMassRange[0] = 1.08; fLambdaMassRange[1] = 1.15;
  
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());   // QA histograms
  //DefineOutput(2, TTree::Class());   // reduced information tree
  //DefineOutput(3, TTree::Class());   // reduced information tree with friends
  DefineOutput(2, TTree::Class());   // reduced information tree with friends
  
  fListHistos.SetName("QAhistograms");
  fListDielectron.SetOwner();
  fListHistos.SetOwner(kFALSE);
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  if (!fListHistos.IsEmpty() || fTree || fFriendTree) return; //already initialised

  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->Init();
    if (die->GetHistogramList()) fListHistos.Add(const_cast<THashList*>(die->GetHistogramList()));
  }  
  if(fV0Histos) fListHistos.Add(const_cast<THashList*>(fV0Histos->GetHistogramList()));

  fFriendTree = new TTree("DstFriendTree","Reduced ESD information");
  fReducedEventFriend = new AliCorrelationReducedEventFriend();
  fFriendTree->Branch("Event",&fReducedEventFriend,16000,99);
  
  fTreeFile = new TFile("dstTree.root", "RECREATE");
  fTree = new TTree("DstTree","Reduced ESD information");
  fReducedEvent = new AliCorrelationReducedEvent();
  fTree->Branch("Event",&fReducedEvent,16000,99);
    
  fFlowTrackCuts = new AliFlowTrackCuts("flow cuts");
  fFlowTrackCuts->SetPID(AliPID::kPion, AliFlowTrackCuts::kTOFbayesian);
  
  PostData(1, &fListHistos);
  //PostData(2, fTree);
  //PostData(3, fFriendTree);
  PostData(2, fFriendTree);
}

//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::UserExec(Option_t *option)
{
  //
  // Main loop. Called for every event
  //  
  //cout << "Event" << endl;
  option = option;
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;
  
  if ( inputHandler->GetPIDResponse() ){
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    AliFatal("This task needs the PID response attached to the input event handler!");
  }
  
  fFlowTrackCuts->SetEvent(InputEvent());
  
  // Was event selected ?
  UInt_t isSelected = AliVEvent::kAny;
  if(fSelectPhysics && inputHandler){
    if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      isSelected&=fTriggerMask;
    }
  }

  if(isSelected==0) {
    return;
  }

  // fill event histograms before event filter
  Double_t values[AliDielectronVarManager::kNMaxValues]={0};
  AliDielectronVarManager::Fill(InputEvent(),values);
  
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    AliDielectronHistos *h=die->GetHistoManager();
    if (h){
      if (h->GetHistogramList()->FindObject("Event_noCuts"))
        h->FillClass("Event_noCuts",AliDielectronVarManager::kNMaxValues,values);
    }
  }
  nextDie.Reset();

  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
  }

  //cout << "Event selected" << endl;
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  //Process event in all AliDielectron instances
  fReducedEvent->ClearEvent();
  fReducedEventFriend->ClearEvent();
  FillEventInfo();
  FillV0PairInfo();
  
  Short_t idie=0;
  while((die=static_cast<AliDielectron*>(nextDie()))){
    die->Process(InputEvent());
    FillDielectronPairInfo(die, idie);
    ++idie;
  }
  nextDie.Reset();
  
  FillTrackInfo();
  FillFriendEventInfo();              // Q-vector calculation
  
  fTree->Fill();
  fFriendTree->Fill();
      
  // if there are candidate pairs, add the information to the reduced tree
  PostData(1, &fListHistos);
  //PostData(2, fTree);
  //PostData(3, fFriendTree);
  PostData(2, fFriendTree);
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FillEventInfo() 
{
  //
  // fill reduced event information
  //
  AliVEvent* event = InputEvent();
  // Was event selected ?
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  UInt_t isSelected = AliVEvent::kAny;
  if(inputHandler){
    if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      isSelected&=fTriggerMask;
    }
  }

  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(event, values);
  
  fReducedEvent->fRunNo       = event->GetRunNumber();
  fReducedEvent->fBC          = event->GetBunchCrossNumber();
  fReducedEvent->fTriggerMask = event->GetTriggerMask();
  fReducedEvent->fIsPhysicsSelection = (isSelected!=0 ? kTRUE : kFALSE);
  AliVVertex* eventVtx = 0x0;
  if(isESD) eventVtx = const_cast<AliESDVertex*>((static_cast<AliESDEvent*>(event))->GetPrimaryVertexTracks());
  if(isAOD) eventVtx = const_cast<AliAODVertex*>((static_cast<AliAODEvent*>(event))->GetPrimaryVertex());
  if(eventVtx) {
    fReducedEvent->fVtx[0] = (isESD ? ((AliESDVertex*)eventVtx)->GetXv() : ((AliAODVertex*)eventVtx)->GetX());
    fReducedEvent->fVtx[1] = (isESD ? ((AliESDVertex*)eventVtx)->GetYv() : ((AliAODVertex*)eventVtx)->GetY());
    fReducedEvent->fVtx[2] = (isESD ? ((AliESDVertex*)eventVtx)->GetZv() : ((AliAODVertex*)eventVtx)->GetZ());
    fReducedEvent->fNVtxContributors = eventVtx->GetNContributors();
  }
  if(isESD) {
    eventVtx = const_cast<AliESDVertex*>((static_cast<AliESDEvent*>(event))->GetPrimaryVertexTPC());
    if(eventVtx) {
      fReducedEvent->fVtxTPC[0] = ((AliESDVertex*)eventVtx)->GetXv();
      fReducedEvent->fVtxTPC[1] = ((AliESDVertex*)eventVtx)->GetYv();
      fReducedEvent->fVtxTPC[2] = ((AliESDVertex*)eventVtx)->GetZv();
      fReducedEvent->fNVtxTPCContributors = eventVtx->GetNContributors();
    }
  }
  
  AliCentrality *centrality = event->GetCentrality();
  if(centrality) {
    fReducedEvent->fCentrality[0] = centrality->GetCentralityPercentile("V0M");
    fReducedEvent->fCentrality[1] = centrality->GetCentralityPercentile("CL1");
    fReducedEvent->fCentrality[2] = centrality->GetCentralityPercentile("TRK");
    fReducedEvent->fCentrality[3] = centrality->GetCentralityPercentile("ZEMvsZDC");    
    fReducedEvent->fCentQuality   = centrality->GetQuality();
  }
  
  fReducedEvent->fNtracks[0] = event->GetNumberOfTracks();
  fReducedEvent->fSPDntracklets = values[AliDielectronVarManager::kNaccTrckltsEsd10Corr];

  AliVVZERO* vzero = event->GetVZEROData();
  for(Int_t i=0;i<64;++i) 
    fReducedEvent->fVZEROMult[i] = vzero->GetMultiplicity(i);  

  AliESDZDC* zdc = (isESD ? (static_cast<AliESDEvent*>(event))->GetESDZDC() : 0x0);
  if(zdc) {
    for(Int_t i=0; i<4; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZN1TowerEnergy()[i];
    for(Int_t i=4; i<8; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZN2TowerEnergy()[i-5];
  }
  
  // EMCAL/PHOS clusters
  FillCaloClusters();
  
  // TODO FMD multiplicities
  
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FillCaloClusters() {
  //
  // Fill info about the calorimeter clusters
  //
  AliVEvent* event = InputEvent();
  Int_t nclusters = event->GetNumberOfCaloClusters();

  fReducedEvent->fNCaloClusters = 0;
  for(Int_t iclus=0; iclus<nclusters; ++iclus) {
    AliVCluster* cluster = event->GetCaloCluster(iclus);
    
    TClonesArray& clusters = *(fReducedEvent->fCaloClusters);
    AliCorrelationReducedCaloCluster *reducedCluster=new(clusters[fReducedEvent->fNCaloClusters]) AliCorrelationReducedCaloCluster();
    
    reducedCluster->fType    = (cluster->IsEMCAL() ? AliCorrelationReducedCaloCluster::kEMCAL : AliCorrelationReducedCaloCluster::kPHOS);
    reducedCluster->fEnergy  = cluster->E();
    reducedCluster->fTrackDx = cluster->GetTrackDx();
    reducedCluster->fTrackDz = cluster->GetTrackDz();
    fReducedEvent->fNCaloClusters += 1;
  }  // end loop over clusters
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FillFriendEventInfo() {
  //
  // Fill event info into the friend tree
  //
  // Add here calculated Q-vector components from all detectors
  for(Int_t idet=0; idet<AliCorrelationReducedEventFriend::kNdetectors; ++idet) {
    fReducedEvent->GetQvector(fReducedEventFriend->fQvector[idet], idet);
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih)
      fReducedEventFriend->fEventPlaneStatus[idet][ih] = AliCorrelationReducedEventFriend::kRaw;
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FillTrackInfo() 
{
  //
  // fill reduced track information
  //
  AliVEvent* event = InputEvent();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());

  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack=0; itrack<ntracks; ++itrack){
    AliVParticle *particle=event->GetTrack(itrack);
    //apply track cuts
    if(fTrackFilter && !fTrackFilter->IsSelected(particle)) continue;
    
    TClonesArray& tracks = *(fReducedEvent->fTracks);
    AliCorrelationReducedTrack *reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliCorrelationReducedTrack();
        
    Double_t values[AliDielectronVarManager::kNMaxValues];
    AliDielectronVarManager::Fill(particle, values);
    reducedParticle->fStatus        = (ULong_t)values[AliDielectronVarManager::kTrackStatus];
    reducedParticle->fGlobalPhi     = values[AliDielectronVarManager::kPhi];
    reducedParticle->fGlobalPt      = values[AliDielectronVarManager::kPt]*values[AliDielectronVarManager::kCharge];
    reducedParticle->fGlobalEta     = values[AliDielectronVarManager::kEta];    
    reducedParticle->fMomentumInner = values[AliDielectronVarManager::kPIn];
    reducedParticle->fDCA[0]        = values[AliDielectronVarManager::kImpactParXY];
    reducedParticle->fDCA[1]        = values[AliDielectronVarManager::kImpactParZ];
    
    reducedParticle->fITSclusterMap = values[AliDielectronVarManager::kITSclusterMap];
    reducedParticle->fITSsignal     = values[AliDielectronVarManager::kITSsignal];
    
    reducedParticle->fTPCNcls      = (UChar_t)values[AliDielectronVarManager::kNclsTPC];
    reducedParticle->fTPCNclsF     = (UChar_t)values[AliDielectronVarManager::kNFclsTPC];
    reducedParticle->fTPCNclsIter1 = (UChar_t)values[AliDielectronVarManager::kNclsTPCiter1];
    reducedParticle->fTPCsignal    = values[AliDielectronVarManager::kTPCsignal];
    reducedParticle->fTPCnSig[0]   = values[AliDielectronVarManager::kTPCnSigmaEle];
    reducedParticle->fTPCnSig[1]   = values[AliDielectronVarManager::kTPCnSigmaPio];
    reducedParticle->fTPCnSig[2]   = values[AliDielectronVarManager::kTPCnSigmaKao];
    reducedParticle->fTPCnSig[3]   = values[AliDielectronVarManager::kTPCnSigmaPro];
    
    reducedParticle->fTOFbeta      = values[AliDielectronVarManager::kTOFbeta];
    reducedParticle->fTOFnSig[0]   = values[AliDielectronVarManager::kTOFnSigmaEle];
    reducedParticle->fTOFnSig[1]   = values[AliDielectronVarManager::kTOFnSigmaPio];
    reducedParticle->fTOFnSig[2]   = values[AliDielectronVarManager::kTOFnSigmaKao];
    reducedParticle->fTOFnSig[3]   = values[AliDielectronVarManager::kTOFnSigmaPro];

    reducedParticle->fTRDpid[0]    = values[AliDielectronVarManager::kTRDprobEle];
    reducedParticle->fTRDpid[1]    = values[AliDielectronVarManager::kTRDprobPio];
    
    if(fFlowTrackFilter) {
      // switch on the first bit if this particle should be used for the event plane
      if(fFlowTrackFilter->IsSelected(particle)) reducedParticle->fFlags |= (1<<0);
    }
    
    if(isESD){
      AliESDtrack *track=static_cast<AliESDtrack*>(particle);
      reducedParticle->fTrackId          = (UShort_t)track->GetID();
      reducedParticle->fTPCCrossedRows   = (UChar_t)track->GetTPCCrossedRows();
      reducedParticle->fTPCClusterMap    = EncodeTPCClusterMap(track);
      const AliExternalTrackParam* tpcInner = track->GetTPCInnerParam();
      reducedParticle->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
      reducedParticle->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
      reducedParticle->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
      reducedParticle->fTRDntracklets[0] = track->GetTRDntracklets();
      reducedParticle->fTRDntracklets[1] = track->GetTRDntrackletsPID();
      fFlowTrackCuts->IsSelected(track);
      AliFlowBayesianPID* bayesResponse = fFlowTrackCuts->GetBayesianResponse();
      reducedParticle->fBayesPID[0] = bayesResponse->GetProb()[AliPID::kPion];
      reducedParticle->fBayesPID[1] = bayesResponse->GetProb()[AliPID::kKaon];
      reducedParticle->fBayesPID[2] = bayesResponse->GetProb()[AliPID::kProton];
      if(track->IsEMCAL()) reducedParticle->fCaloClusterId = track->GetEMCALcluster();
      if(track->IsPHOS()) reducedParticle->fCaloClusterId = track->GetPHOScluster();
    }
    if(isAOD) {
      AliAODTrack *track=static_cast<AliAODTrack*>(particle);
      reducedParticle->fTrackId = track->GetID(); 
      if(track->IsEMCAL()) reducedParticle->fCaloClusterId = track->GetEMCALcluster();
      if(track->IsPHOS()) reducedParticle->fCaloClusterId = track->GetPHOScluster();
    }

    fReducedEvent->fNtracks[1] += 1;
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FillDielectronPairInfo(AliDielectron* die, Short_t iDie) 
{
  //
  // fill reduced pair information
  //
  Bool_t hasMC=AliDielectronMC::Instance()->HasMC();

  for(Int_t iType=0; iType<3; ++iType) {
    
    const TObjArray* array = die->GetPairArray(iType);
    if(!array || array->GetEntriesFast()==0) continue;
    cout << "dielectron tracks (idie/type/n): " << iDie << " / 0 / " << die->GetTrackArray(0)->GetEntriesFast() << endl;
    cout << "dielectron tracks (idie/type/n): " << iDie << " / 1 / " << die->GetTrackArray(1)->GetEntriesFast() << endl;

    cout << "diele idie/type/n-candidates: " << iDie << " / " << iType << " / " << array->GetEntriesFast() << endl;
    for(Int_t iCandidate=0; iCandidate<array->GetEntriesFast(); ++iCandidate) {
      if(iCandidate%100==0) cout << "iCandidate = " << iCandidate << endl;
      AliDielectronPair* pair = (AliDielectronPair*)array->At(iCandidate);
      Double_t values[AliDielectronVarManager::kNMaxValues];
      AliDielectronVarManager::Fill(pair, values);
      
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliCorrelationReducedPair *reducedParticle= 
         new (tracks[fReducedEvent->fNV0candidates[1]+fReducedEvent->fNDielectronCandidates]) AliCorrelationReducedPair();
      // !!! hardcoded flag for dielectron id 
      reducedParticle->fCandidateId  = (iDie==0 ? AliCorrelationReducedPair::kJpsiToEE : AliCorrelationReducedPair::kPhiToKK);
      reducedParticle->fPairType     = values[AliDielectronVarManager::kPairType];
      reducedParticle->fLegIds[0]    = (UShort_t)(static_cast<AliVTrack*>(pair->GetFirstDaughter()))->GetID();
      reducedParticle->fLegIds[1]    = (UShort_t)(static_cast<AliVTrack*>(pair->GetSecondDaughter()))->GetID();
      reducedParticle->fMass[0]      = values[AliDielectronVarManager::kM];
      reducedParticle->fMass[1]      = -999.;
      reducedParticle->fMass[2]      = -999.;
      reducedParticle->fPhi          = values[AliDielectronVarManager::kPhi];  // in the [-pi,pi] interval
      if(reducedParticle->fPhi<0.0) reducedParticle->fPhi = 2.0*TMath::Pi() + reducedParticle->fPhi;  // converted to [0,2pi]
      reducedParticle->fPt           = values[AliDielectronVarManager::kPt];
      reducedParticle->fEta          = values[AliDielectronVarManager::kEta];
      reducedParticle->fLxy          = values[AliDielectronVarManager::kPseudoProperTime];
      reducedParticle->fOpeningAngle = values[AliDielectronVarManager::kOpeningAngle];     
     
      reducedParticle->fMCid         = 0;
      if(hasMC) {
	AliDielectronMC::Instance()->ConnectMCEvent();
	const TObjArray* mcSignals = die->GetMCSignals();
	for(Int_t iSig=0; iSig<mcSignals->GetEntries(); ++iSig) {
	  if(iSig>31) break;
	  AliDielectronMC *mc=AliDielectronMC::Instance();
	  if(mc->IsMCTruth(pair, (AliDielectronSignalMC*)mcSignals->At(iSig))) {
	    reducedParticle->fMCid = reducedParticle->fMCid | (1<<iSig);
	  }
	}
      }   // end if has MC
      fReducedEvent->fNDielectronCandidates += 1;
    }    // end loop over candidates
  }    // end loop over pair type
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FillV0PairInfo() 
{
  //
  // fill reduced pair information
  //
  AliESDEvent* esd = (AliESDEvent*)InputEvent();
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*primaryVertex);
  
  fReducedEvent->fNV0candidates[0] = InputEvent()->GetNumberOfV0s();
  
  Double_t valuesPos[AliDielectronVarManager::kNMaxValues];
  Double_t valuesNeg[AliDielectronVarManager::kNMaxValues];
   
  cout << "n total V0s: " << InputEvent()->GetNumberOfV0s() << endl;
  for(Int_t iV0=0; iV0<InputEvent()->GetNumberOfV0s(); ++iV0) {   // loop over V0s
    if(iV0%1000==0) cout << "iV0 = " << iV0 << endl;
    AliESDv0 *v0 = esd->GetV0(iV0);
       
    AliESDtrack* legPos = esd->GetTrack(v0->GetPindex());
    AliESDtrack* legNeg = esd->GetTrack(v0->GetNindex());
 
    if(legPos->GetSign() == legNeg->GetSign()) {
      //cout << "V0 rejected: same sign legs" << endl;
      continue;
    }

    Bool_t v0ChargesAreCorrect = (legPos->GetSign()==+1 ? kTRUE : kFALSE);
    legPos = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetNindex()) : legPos);
    legNeg = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetPindex()) : legNeg); 
    
    // apply the K0s filter
    Bool_t goodK0s = kTRUE;
    if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(legPos) || !fK0sPionCuts->IsSelected(legNeg))) goodK0s = kFALSE;
    if(fK0sCuts) {
      TList k0sList;
      k0sList.Add(v0);
      k0sList.Add(legPos); k0sList.Add(legNeg);
      k0sList.Add(const_cast<AliESDVertex*>(primaryVertex));
      if(!fK0sCuts->IsSelected(&k0sList)) goodK0s = kFALSE;
    }
    
    // apply the lambda filter
    Bool_t goodLambda = kTRUE;
    if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legPos)) goodLambda = kFALSE;
    if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legNeg)) goodLambda = kFALSE;
    Bool_t goodALambda = kTRUE;
    if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legNeg)) goodALambda = kFALSE;
    if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legPos)) goodALambda = kFALSE;
    if(fLambdaCuts) {
      TList lambdaList;
      lambdaList.Add(v0);
      lambdaList.Add(legPos); lambdaList.Add(legNeg);
      lambdaList.Add(const_cast<AliESDVertex*>(primaryVertex));
      if(!fLambdaCuts->IsSelected(&lambdaList)) {goodLambda = kFALSE; goodALambda = kFALSE;}
    }
    
    if(!(goodK0s || goodLambda || goodALambda)) continue;
    
    // Fill the V0 information into the tree for 3 hypothesis: K0s, Lambda, Anti-Lambda
    AliCorrelationReducedPair* k0sReducedPair     = FillV0PairInfo(v0, AliCorrelationReducedPair::kK0sToPiPi,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliCorrelationReducedPair* lambdaReducedPair  = FillV0PairInfo(v0, AliCorrelationReducedPair::kLambda0ToPPi,  legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliCorrelationReducedPair* alambdaReducedPair = FillV0PairInfo(v0, AliCorrelationReducedPair::kALambda0ToPPi, legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);

    if(goodK0s && k0sReducedPair->fMass[0]>fK0sMassRange[0] && k0sReducedPair->fMass[0]<fK0sMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliCorrelationReducedPair *goodK0sPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliCorrelationReducedPair(*k0sReducedPair);
      goodK0sPair->fMass[0] = k0sReducedPair->fMass[0];
      goodK0sPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodK0sPair->fMass[2] = alambdaReducedPair->fMass[0];
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodK0s=kFALSE;}
    if(goodLambda && lambdaReducedPair->fMass[0]>fLambdaMassRange[0] && lambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliCorrelationReducedPair *goodLambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliCorrelationReducedPair(*lambdaReducedPair);
      fReducedEvent->fNV0candidates[1] += 1;
      goodLambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodLambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
    } else {goodLambda=kFALSE;}
    if(goodALambda && alambdaReducedPair->fMass[0]>fLambdaMassRange[0] && alambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliCorrelationReducedPair *goodALambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliCorrelationReducedPair(*alambdaReducedPair);
      fReducedEvent->fNV0candidates[1] += 1;  
      goodALambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodALambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
    } else {goodALambda = kFALSE;}
    delete k0sReducedPair;
    delete lambdaReducedPair;
    delete alambdaReducedPair;

    if(!(goodK0s || goodLambda || goodALambda)) continue;
    
    //  Fill histograms and the CF container
    AliDielectronVarManager::Fill(legPos, valuesPos);
    AliDielectronVarManager::Fill(legNeg, valuesNeg);
    
    if(fV0Histos && fV0Histos->GetHistogramList()->FindObject("V0Track_Pos"))
      fV0Histos->FillClass("V0Track_Pos", AliDielectronVarManager::kNMaxValues, valuesPos);
    if(fV0Histos && fV0Histos->GetHistogramList()->FindObject("V0Track_Neg"))
      fV0Histos->FillClass("V0Track_Neg", AliDielectronVarManager::kNMaxValues, valuesNeg);    
  }   // end loop over V0s
}


//_________________________________________________________________________________
AliCorrelationReducedPair* AliAnalysisTaskCorrelationTree::FillV0PairInfo(AliESDv0* v0, Int_t id, 
						    AliESDtrack* legPos, AliESDtrack* legNeg,
						    AliKFVertex* vtxKF, Bool_t chargesAreCorrect) {
  //
  // Create a reduced V0 object and fill it
  //
  AliCorrelationReducedPair* reducedPair=new AliCorrelationReducedPair();  
  reducedPair->fCandidateId = id;
  reducedPair->fPairType    = v0->GetOnFlyStatus();    // on the fly status
  //reducedPair->fOnTheFly    = v0->GetOnFlyStatus();
  reducedPair->fLegIds[0]   = legPos->GetID();
  reducedPair->fLegIds[1]   = legNeg->GetID();
  if(!reducedPair->fPairType) {    // offline
    UInt_t pidPos = AliPID::kPion;
    if(id==AliCorrelationReducedPair::kLambda0ToPPi) pidPos = AliPID::kProton;
    UInt_t pidNeg = AliPID::kPion;
    if(id==AliCorrelationReducedPair::kALambda0ToPPi) pidNeg = AliPID::kProton;
    reducedPair->fMass[0]      = v0->GetEffMass(pidPos, pidNeg);
    reducedPair->fPhi          = v0->Phi();
    if(reducedPair->fPhi<0.0) reducedPair->fPhi = 2.0*TMath::Pi() + reducedPair->fPhi;  // converted to [0,2pi]
    reducedPair->fPt           = v0->Pt();
    reducedPair->fEta          = v0->Eta();
    reducedPair->fLxy          = 0.0;           //TODO
    reducedPair->fOpeningAngle = 0.0;
  }
  else {
    const AliExternalTrackParam *negHelix=v0->GetParamN();
    const AliExternalTrackParam *posHelix=v0->GetParamP();
    if(!chargesAreCorrect) {
      negHelix = v0->GetParamP();
      posHelix = v0->GetParamN();
    }
    AliKFParticle negKF(*(negHelix),(id==AliCorrelationReducedPair::kALambda0ToPPi ? -2212 : -211));
    AliKFParticle posKF(*(posHelix),(id==AliCorrelationReducedPair::kLambda0ToPPi ? +2212 : +211));
    AliKFParticle v0Refit;
    v0Refit += negKF;
    v0Refit += posKF;
    Double_t massFit=0.0, massErrFit=0.0;
    v0Refit.GetMass(massFit,massErrFit);
    reducedPair->fMass[0] = massFit;
    reducedPair->fPhi          = v0Refit.GetPhi();
    if(reducedPair->fPhi<0.0) reducedPair->fPhi = 2.0*TMath::Pi() + reducedPair->fPhi;  // converted to [0,2pi]
    reducedPair->fPt           = v0Refit.GetPt();
    reducedPair->fEta          = v0Refit.GetEta();
    reducedPair->fLxy          = v0Refit.GetPseudoProperDecayTime(*vtxKF, massFit);
    reducedPair->fOpeningAngle = negKF.GetAngle(posKF);
  }
  return reducedPair;
}


//_________________________________________________________________________________
UChar_t AliAnalysisTaskCorrelationTree::EncodeTPCClusterMap(AliESDtrack* track) {
  //
  // Encode the TPC cluster map into an UChar_t
  // Divide the 159 bits from the bit map into 8 groups of adiacent clusters
  // For each group enable its corresponding bit if in that group there are more clusters compared to
  // a threshold.
  //
  const UChar_t threshold=5;
  TBits tpcClusterMap = track->GetTPCClusterMap();
  UChar_t map=0;
  UChar_t n=0;
  UChar_t j=0;
  for(UChar_t i=0; i<8; ++i) {
    n=0;
    for(j=i*20; j<(i+1)*20 && j<159; ++j) n+=tpcClusterMap.TestBitNumber(j);
    if(n>=threshold) map |= (1<<i);
  }
  return map;
}


//_________________________________________________________________________________
void AliAnalysisTaskCorrelationTree::FinishTaskOutput()
{
  //
  // Finish Task 
  //
  fTreeFile->Write();
  fTreeFile->Close();
}
