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
#include <AliESDv0KineCuts.h>
#include <AliVCluster.h>
#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronMC.h"
#include "AliDielectronVarManager.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowBayesianPID.h"

#include "AliReducedEvent.h"
#include "AliAnalysisTaskReducedTree.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskReducedTree)


//_________________________________________________________________________________
AliAnalysisTaskReducedTree::AliAnalysisTaskReducedTree() :
  AliAnalysisTaskSE(),
  fListDielectron(),
  fListHistos(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fFillTrackInfo(kTRUE),
  fFillDielectronInfo(kTRUE),
  fFillV0Info(kTRUE),
  fFillGammaConversions(kTRUE),
  fFillK0s(kTRUE),
  fFillLambda(kTRUE),
  fFillALambda(kTRUE),
  fFillCaloClusterInfo(kTRUE),
  fFillFriendInfo(kTRUE),
  fEventFilter(0x0),
  fTrackFilter(0x0),
  fFlowTrackFilter(0x0),
  fK0sCuts(0x0),
  fLambdaCuts(0x0),
  fGammaConvCuts(0x0),
  fK0sPionCuts(0x0),
  fLambdaProtonCuts(0x0),
  fLambdaPionCuts(0x0),
  fGammaElectronCuts(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fGammaMassRange(),
  fV0Histos(0x0),
  fTreeFile(0x0),
  fTree(0x0),
  fFriendTreeFile(0x0),
  fFriendTree(0x0),
  fReducedEvent(0x0),
  fReducedEventFriend(0x0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskReducedTree::AliAnalysisTaskReducedTree(const char *name) :
  AliAnalysisTaskSE(name),
  fListDielectron(),
  fListHistos(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fFillTrackInfo(kTRUE),
  fFillDielectronInfo(kTRUE),
  fFillV0Info(kTRUE),
  fFillGammaConversions(kTRUE),
  fFillK0s(kTRUE),
  fFillLambda(kTRUE),
  fFillALambda(kTRUE),
  fFillCaloClusterInfo(kTRUE),
  fFillFriendInfo(kTRUE),
  fEventFilter(0x0),
  fTrackFilter(0x0),
  fFlowTrackFilter(0x0),
  fK0sCuts(0x0),
  fLambdaCuts(0x0),
  fGammaConvCuts(0x0),
  fK0sPionCuts(0x0),
  fLambdaProtonCuts(0x0),
  fLambdaPionCuts(0x0),
  fGammaElectronCuts(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fGammaMassRange(),
  fV0Histos(0x0),
  fTreeFile(0x0),
  fTree(0x0),
  fFriendTreeFile(0x0),
  fFriendTree(0x0),
  fReducedEvent(0x0),
  fReducedEventFriend(0x0)
{
  //
  // Constructor
  //
  fK0sMassRange[0] = 0.4; fK0sMassRange[1] = 0.6;
  fLambdaMassRange[0] = 1.08; fLambdaMassRange[1] = 1.15;
  fGammaMassRange[0] = 0.0; fGammaMassRange[1] = 0.1;
  
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());   // QA histograms
  DefineOutput(2, TTree::Class());   // reduced information tree
  //if(fFillFriendInfo) DefineOutput(3, TTree::Class());   // reduced information tree with friends
  //DefineOutput(2, TTree::Class());   // reduced information tree with friends
  //DefineOutput(2, TTree::Class());   // reduced information tree  

  fListHistos.SetName("QAhistograms");
  fListDielectron.SetOwner();
  fListHistos.SetOwner(kFALSE);
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::UserCreateOutputObjects()
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

  if(fFillFriendInfo) {
    fFriendTree = new TTree("DstFriendTree","Reduced ESD information");
    fReducedEventFriend = new AliReducedEventFriend();
    fFriendTree->Branch("Event",&fReducedEventFriend,16000,99);
  }
  
  //fTreeFile = new TFile("dstTree.root", "RECREATE");
  fTree = new TTree("DstTree","Reduced ESD information");
  fReducedEvent = new AliReducedEvent("DstEvent");
  fTree->Branch("Event",&fReducedEvent,16000,99);
    
  PostData(1, &fListHistos);
  PostData(2, fTree);
  //if(fFillFriendInfo) PostData(3, fFriendTree);
  //PostData(2, fFriendTree);
  //PostData(1, fTree);
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::UserExec(Option_t *option)
{
  //
  // Main loop. Called for every event
  //  
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

  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  //Process event in all AliDielectron instances
  fReducedEvent->ClearEvent();
  if(fFillFriendInfo) fReducedEventFriend->ClearEvent();
  FillEventInfo();
  if(fFillV0Info) FillV0PairInfo();
  
  Short_t idie=0;
  if(fFillDielectronInfo) {
    while((die=static_cast<AliDielectron*>(nextDie()))){
      die->Process(InputEvent());
      FillDielectronPairInfo(die, idie);
      ++idie;
    }
  }
  nextDie.Reset();
  
  if(fFillTrackInfo) FillTrackInfo();
  if(fFillFriendInfo) FillFriendEventInfo();              // Q-vector calculation
  
  fTree->Fill();
  if(fFillFriendInfo) fFriendTree->Fill();
      
  // if there are candidate pairs, add the information to the reduced tree
  PostData(1, &fListHistos);
  PostData(2, fTree);
  //if(fFillFriendInfo) PostData(3, fFriendTree);
  //PostData(2, fFriendTree);
  //PostData(2, fTree);
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::FillEventInfo() 
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
  fReducedEvent->fSPDntracklets = (UChar_t)values[AliDielectronVarManager::kNaccTrckltsEsd10Corr];

  AliVVZERO* vzero = event->GetVZEROData();
  for(Int_t i=0;i<64;++i) 
    fReducedEvent->fVZEROMult[i] = vzero->GetMultiplicity(i);  

  AliESDZDC* zdc = (isESD ? (static_cast<AliESDEvent*>(event))->GetESDZDC() : 0x0);
  if(zdc) {
    for(Int_t i=0; i<4; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZN1TowerEnergy()[i];
    for(Int_t i=4; i<8; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZN2TowerEnergy()[i-4];
  }
  
  // EMCAL/PHOS clusters
  if(fFillCaloClusterInfo) FillCaloClusters();
  
  // TODO FMD multiplicities
  
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::FillCaloClusters() {
  //
  // Fill info about the calorimeter clusters
  //
  AliVEvent* event = InputEvent();
  Int_t nclusters = event->GetNumberOfCaloClusters();

  fReducedEvent->fNCaloClusters = 0;
  for(Int_t iclus=0; iclus<nclusters; ++iclus) {
    AliVCluster* cluster = event->GetCaloCluster(iclus);
    
    TClonesArray& clusters = *(fReducedEvent->fCaloClusters);
    AliReducedCaloCluster *reducedCluster=new(clusters[fReducedEvent->fNCaloClusters]) AliReducedCaloCluster();
    
    reducedCluster->fType    = (cluster->IsEMCAL() ? AliReducedCaloCluster::kEMCAL : AliReducedCaloCluster::kPHOS);
    reducedCluster->fEnergy  = cluster->E();
    reducedCluster->fTrackDx = cluster->GetTrackDx();
    reducedCluster->fTrackDz = cluster->GetTrackDz();
    fReducedEvent->fNCaloClusters += 1;
  }  // end loop over clusters
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::FillFriendEventInfo() {
  //
  // Fill event info into the friend tree
  //
  // Add here calculated Q-vector components from all detectors
  for(Int_t idet=0; idet<AliReducedEventFriend::kNdetectors; ++idet) {
    fReducedEvent->GetQvector(fReducedEventFriend->fQvector[idet], idet);
    for(Int_t ih=0; ih<fgkNMaxHarmonics; ++ih)
      fReducedEventFriend->fEventPlaneStatus[idet][ih] = AliReducedEventFriend::kRaw;
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::FillTrackInfo() 
{
  //
  // fill reduced track information
  //
  AliVEvent* event = InputEvent();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());

  // find all the tracks which belong to a V0 stored in the reduced event
  UShort_t trackIdsV0[4][20000]={{0}};
  Int_t nV0LegsTagged[4] = {0};
  Bool_t leg1Found[4]; Bool_t leg2Found[4];
  for(Int_t iv0=0;iv0<fReducedEvent->fNV0candidates[1];++iv0) {
    AliReducedPair* pair = fReducedEvent->GetV0Pair(iv0);
    Int_t pairId = 0;
    if(pair->fCandidateId==AliReducedPair::kGammaConv) pairId=0;
    if(pair->fCandidateId==AliReducedPair::kK0sToPiPi) pairId=1;
    if(pair->fCandidateId==AliReducedPair::kLambda0ToPPi) pairId=2;
    if(pair->fCandidateId==AliReducedPair::kALambda0ToPPi) pairId=3;
    leg1Found[pairId] = kFALSE; leg2Found[pairId] = kFALSE;
    for(Int_t it=0;it<nV0LegsTagged[pairId];++it) {
      if(trackIdsV0[pairId][it]==pair->fLegIds[0]) leg1Found[pairId]=kTRUE;
      if(trackIdsV0[pairId][it]==pair->fLegIds[1]) leg2Found[pairId]=kTRUE;
    }
    // if the legs of this V0 were not already stored then add them now to the list
    if(!leg1Found[pairId]) {trackIdsV0[pairId][nV0LegsTagged[pairId]] = pair->fLegIds[0]; ++nV0LegsTagged[pairId];}
    if(!leg2Found[pairId]) {trackIdsV0[pairId][nV0LegsTagged[pairId]] = pair->fLegIds[1]; ++nV0LegsTagged[pairId];}
  }
  
  AliESDtrack* esdTrack=0;
  AliAODTrack* aodTrack=0;
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t trackId = 0; Bool_t usedForV0[4] = {kFALSE}; Bool_t usedForV0Or = kFALSE;
  for(Int_t itrack=0; itrack<ntracks; ++itrack){
    AliVParticle *particle=event->GetTrack(itrack);
    if(isESD) {
      esdTrack=static_cast<AliESDtrack*>(particle);
      trackId = esdTrack->GetID();
    }
    if(isAOD) {
      aodTrack=static_cast<AliAODTrack*>(particle);
      trackId = aodTrack->GetID();
    }
    // check whether this track belongs to a V0 stored in the reduced event
    usedForV0Or = kFALSE;
    for(Int_t i=0; i<4; ++i) {
      usedForV0[i] = kFALSE;
      for(Int_t ii=0; ii<nV0LegsTagged[i]; ++ii) {
        if(UShort_t(trackId)==trackIdsV0[i][ii]) {
          usedForV0[i] = kTRUE;
          break;
        }
      }
      usedForV0Or |= usedForV0[i];
    }
    
    //apply track cuts
    if(!usedForV0Or && fTrackFilter && !fTrackFilter->IsSelected(particle)) continue;
       
    TClonesArray& tracks = *(fReducedEvent->fTracks);
    AliReducedTrack *reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliReducedTrack();
        
    Double_t values[AliDielectronVarManager::kNMaxValues];
    AliDielectronVarManager::Fill(particle, values);
    reducedParticle->fStatus        = (ULong_t)values[AliDielectronVarManager::kTrackStatus];
    reducedParticle->fGlobalPhi     = values[AliDielectronVarManager::kPhi];
    reducedParticle->fGlobalPt      = values[AliDielectronVarManager::kPt]*values[AliDielectronVarManager::kCharge];
    reducedParticle->fGlobalEta     = values[AliDielectronVarManager::kEta];    
    reducedParticle->fMomentumInner = values[AliDielectronVarManager::kPIn];
    reducedParticle->fDCA[0]        = values[AliDielectronVarManager::kImpactParXY];
    reducedParticle->fDCA[1]        = values[AliDielectronVarManager::kImpactParZ];
    
    reducedParticle->fITSclusterMap = (UChar_t)values[AliDielectronVarManager::kITSclusterMap];
    reducedParticle->fITSsignal     = values[AliDielectronVarManager::kITSsignal];
    reducedParticle->fITSnSig[0]    = values[AliDielectronVarManager::kITSnSigmaEle];
    reducedParticle->fITSnSig[1]    = values[AliDielectronVarManager::kITSnSigmaPio];
    reducedParticle->fITSnSig[2]    = values[AliDielectronVarManager::kITSnSigmaKao];
    reducedParticle->fITSnSig[3]    = values[AliDielectronVarManager::kITSnSigmaPro];
    
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
    // switch bits to show wheter the track is used in a V0
    for(Int_t iV0type=0;iV0type<4;++iV0type) {
      if(usedForV0[iV0type]) reducedParticle->fFlags |= (1<<(iV0type+1));
    }
    
    if(isESD){
      reducedParticle->fTrackId          = (UShort_t)esdTrack->GetID();
      reducedParticle->fTPCCrossedRows   = (UChar_t)esdTrack->GetTPCCrossedRows();
      reducedParticle->fTPCClusterMap    = EncodeTPCClusterMap(esdTrack);
      const AliExternalTrackParam* tpcInner = esdTrack->GetTPCInnerParam();
      reducedParticle->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
      reducedParticle->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
      reducedParticle->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
      reducedParticle->fTRDntracklets[0] = esdTrack->GetTRDntracklets();
      reducedParticle->fTRDntracklets[1] = esdTrack->GetTRDntrackletsPID();
      for(Int_t idx=0; idx<3; ++idx) if(esdTrack->GetKinkIndex(idx)>0) reducedParticle->fFlags |= (1<<(5+idx));
      if(esdTrack->IsEMCAL()) reducedParticle->fCaloClusterId = esdTrack->GetEMCALcluster();
      if(esdTrack->IsPHOS()) reducedParticle->fCaloClusterId = esdTrack->GetPHOScluster();
    }
    if(isAOD) {
      reducedParticle->fTrackId = aodTrack->GetID(); 
      if(aodTrack->IsEMCAL()) reducedParticle->fCaloClusterId = aodTrack->GetEMCALcluster();
      if(aodTrack->IsPHOS()) reducedParticle->fCaloClusterId = aodTrack->GetPHOScluster();
      if(values[AliDielectronVarManager::kKinkIndex0]>0.0) reducedParticle->fFlags |= (1<<5);
    }

    fReducedEvent->fNtracks[1] += 1;
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::FillDielectronPairInfo(AliDielectron* die, Short_t iDie) 
{
  //
  // fill reduced pair information
  //
  Bool_t hasMC=AliDielectronMC::Instance()->HasMC();

  for(Int_t iType=0; iType<3; ++iType) {
    
    const TObjArray* array = die->GetPairArray(iType);
    if(!array || array->GetEntriesFast()==0) continue;
    
    for(Int_t iCandidate=0; iCandidate<array->GetEntriesFast(); ++iCandidate) {
      AliDielectronPair* pair = (AliDielectronPair*)array->At(iCandidate);
      Double_t values[AliDielectronVarManager::kNMaxValues];
      AliDielectronVarManager::Fill(pair, values);
      
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *reducedParticle= 
         new (tracks[fReducedEvent->fNV0candidates[1]+fReducedEvent->fNDielectronCandidates]) AliReducedPair();
      // !!! hardcoded flag for dielectron id 
      reducedParticle->fCandidateId  = (iDie==0 ? AliReducedPair::kJpsiToEE : AliReducedPair::kPhiToKK);
      reducedParticle->fPairType     = (Char_t)values[AliDielectronVarManager::kPairType];
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
      reducedParticle->fLxyErr       = values[AliDielectronVarManager::kPseudoProperTimeErr];
      reducedParticle->fPointingAngle = values[AliDielectronVarManager::kCosPointingAngle];
      
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
void AliAnalysisTaskReducedTree::FillV0PairInfo() 
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
  
  fGammaConvCuts->SetEvent(esd);
  fGammaConvCuts->SetPrimaryVertex(&primaryVertexKF);
  
  Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
  for(Int_t iV0=0; iV0<InputEvent()->GetNumberOfV0s(); ++iV0) {   // loop over V0s
    AliESDv0 *v0 = esd->GetV0(iV0);
       
    AliESDtrack* legPos = esd->GetTrack(v0->GetPindex());
    AliESDtrack* legNeg = esd->GetTrack(v0->GetNindex());
 
    if(legPos->GetSign() == legNeg->GetSign()) {
      continue;
    }

    Bool_t v0ChargesAreCorrect = (legPos->GetSign()==+1 ? kTRUE : kFALSE);
    legPos = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetNindex()) : legPos);
    legNeg = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetPindex()) : legNeg); 
    
    // apply the K0s filter
    Bool_t goodK0s = kTRUE;
    if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(legPos) || !fK0sPionCuts->IsSelected(legNeg))) goodK0s = kFALSE;
    if(goodK0s && fK0sCuts) {
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
    
    // apply the gamma conversion filter
    Bool_t goodGamma = kTRUE;
    //cout << "fGammaElectronCuts " << fGammaElectronCuts << endl;
    if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(legPos) || !fGammaElectronCuts->IsSelected(legNeg))) goodGamma = kFALSE;
    //cout << "goodGamma1 " << goodGamma << endl;
    goodGamma = goodGamma && fGammaConvCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    //cout << "goodGamma2 " << goodGamma << endl;
    //cout << "pdg V0/p/n: " << pdgV0 << "/" << pdgP << "/" << pdgN << endl;
    if(pdgV0!=22 || TMath::Abs(pdgP)!=11 || TMath::Abs(pdgN)!=11) goodGamma = kFALSE;
    //cout << "goodGamma3 " << goodGamma << endl;
    
    if(!(goodK0s || goodLambda || goodALambda || goodGamma)) continue;
    
    // Fill the V0 information into the tree for 4 hypothesis: K0s, Lambda, Anti-Lambda and gamma conversion
    AliReducedPair* k0sReducedPair     = FillV0PairInfo(v0, AliReducedPair::kK0sToPiPi,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPair* lambdaReducedPair  = FillV0PairInfo(v0, AliReducedPair::kLambda0ToPPi,  legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPair* alambdaReducedPair = FillV0PairInfo(v0, AliReducedPair::kALambda0ToPPi, legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPair* gammaReducedPair   = FillV0PairInfo(v0, AliReducedPair::kGammaConv,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    
    if(fFillK0s && goodK0s && k0sReducedPair->fMass[0]>fK0sMassRange[0] && k0sReducedPair->fMass[0]<fK0sMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodK0sPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*k0sReducedPair);
      goodK0sPair->fMass[0] = k0sReducedPair->fMass[0];
      goodK0sPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodK0sPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodK0sPair->fMass[3] = gammaReducedPair->fMass[0];
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodK0s=kFALSE;}
    if(fFillLambda && goodLambda && lambdaReducedPair->fMass[0]>fLambdaMassRange[0] && lambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodLambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*lambdaReducedPair);
      fReducedEvent->fNV0candidates[1] += 1;
      goodLambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodLambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[3] = gammaReducedPair->fMass[0];
    } else {goodLambda=kFALSE;}
    if(fFillALambda && goodALambda && alambdaReducedPair->fMass[0]>fLambdaMassRange[0] && alambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodALambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*alambdaReducedPair);
      fReducedEvent->fNV0candidates[1] += 1;  
      goodALambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodALambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[3] = gammaReducedPair->fMass[0];
    } else {goodALambda = kFALSE;}
    //cout << "gamma mass: " << gammaReducedPair->fMass[0] << endl;
    if(fFillGammaConversions && goodGamma && gammaReducedPair->fMass[0]>fGammaMassRange[0] && gammaReducedPair->fMass[0]<fGammaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodGammaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*gammaReducedPair);
      goodGammaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodGammaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodGammaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodGammaPair->fMass[3] = gammaReducedPair->fMass[0];
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodGamma=kFALSE;}
    delete k0sReducedPair;
    delete lambdaReducedPair;
    delete alambdaReducedPair;
    delete gammaReducedPair;
    
    if(!(goodK0s || goodLambda || goodALambda || goodGamma)) continue;
    
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
AliReducedPair* AliAnalysisTaskReducedTree::FillV0PairInfo(AliESDv0* v0, Int_t id, 
						    AliESDtrack* legPos, AliESDtrack* legNeg,
						    AliKFVertex* vtxKF, Bool_t chargesAreCorrect) {
  //
  // Create a reduced V0 object and fill it
  //
  AliReducedPair* reducedPair=new AliReducedPair();  
  reducedPair->fCandidateId = id;
  reducedPair->fPairType    = v0->GetOnFlyStatus();    // on the fly status
  //reducedPair->fOnTheFly    = v0->GetOnFlyStatus();
  reducedPair->fLegIds[0]   = legPos->GetID();
  reducedPair->fLegIds[1]   = legNeg->GetID();
  if(!reducedPair->fPairType) {    // offline
    UInt_t pidPos = AliPID::kPion;
    if(id==AliReducedPair::kLambda0ToPPi) pidPos = AliPID::kProton;
    if(id==AliReducedPair::kGammaConv) pidPos = AliPID::kElectron;
    UInt_t pidNeg = AliPID::kPion;
    if(id==AliReducedPair::kALambda0ToPPi) pidNeg = AliPID::kProton;
    if(id==AliReducedPair::kGammaConv) pidNeg = AliPID::kElectron;
    reducedPair->fMass[0]      = v0->GetEffMass(pidPos, pidNeg);
    reducedPair->fPhi          = v0->Phi();
    if(reducedPair->fPhi<0.0) reducedPair->fPhi = 2.0*TMath::Pi() + reducedPair->fPhi;  // converted to [0,2pi]
    reducedPair->fPt           = v0->Pt();
    reducedPair->fEta          = v0->Eta();
    reducedPair->fLxy          = v0->GetRr();
    reducedPair->fPointingAngle = v0->GetV0CosineOfPointingAngle(vtxKF->GetX(), vtxKF->GetY(), vtxKF->GetZ());
  }
  else {
    const AliExternalTrackParam *negHelix=v0->GetParamN();
    const AliExternalTrackParam *posHelix=v0->GetParamP();
    if(!chargesAreCorrect) {
      negHelix = v0->GetParamP();
      posHelix = v0->GetParamN();
    }
    Int_t pdgPos = 211;
    if(id==AliReducedPair::kLambda0ToPPi) pdgPos = 2212;
    if(id==AliReducedPair::kGammaConv) pdgPos = -11;
    Int_t pdgNeg = -211;
    if(id==AliReducedPair::kALambda0ToPPi) pdgNeg = -2212;
    if(id==AliReducedPair::kGammaConv) pdgNeg = 11;
    AliKFParticle negKF(*(negHelix), pdgPos);
    AliKFParticle posKF(*(posHelix), pdgNeg);
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
    Double_t deltaPos[3];
    deltaPos[0] = v0Refit.GetX() - vtxKF->GetX(); deltaPos[1] = v0Refit.GetY() - vtxKF->GetY(); deltaPos[2] = v0Refit.GetZ() - vtxKF->GetZ();
    Double_t momV02 = v0Refit.GetPx()*v0Refit.GetPx() + v0Refit.GetPy()*v0Refit.GetPy() + v0Refit.GetPz()*v0Refit.GetPz();
    Double_t deltaPos2 = deltaPos[0]*deltaPos[0] + deltaPos[1]*deltaPos[1] + deltaPos[2]*deltaPos[2];
    reducedPair->fPointingAngle = (deltaPos[0]*v0Refit.GetPx() + deltaPos[1]*v0Refit.GetPy() + deltaPos[2]*v0Refit.GetPz()) / 
                                  TMath::Sqrt(momV02*deltaPos2);
  }
  return reducedPair;
}


//_________________________________________________________________________________
UChar_t AliAnalysisTaskReducedTree::EncodeTPCClusterMap(AliESDtrack* track) {
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
void AliAnalysisTaskReducedTree::FinishTaskOutput()
{
  //
  // Finish Task 
  //
  //fTreeFile->Write();
  //fTreeFile->Close();
}
