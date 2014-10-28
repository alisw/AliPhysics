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
#include <AliESDHeader.h>
#include <AliAODHeader.h>
#include <AliAODTrack.h>
//#include <AliAODForwardMult.h>
//#include <AliForwardUtil.h>
#include <AliTriggerAnalysis.h>
#include <AliESDtrackCuts.h>
#include <AliVZDC.h>
#include <AliESDv0.h>
#include <AliESDv0Cuts.h>
#include <AliESDv0KineCuts.h>
#include <AliESDFMD.h>
#include <AliVCluster.h>
#include <AliAODTracklets.h>
#include <AliMultiplicity.h>
#include <AliPIDResponse.h>
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
  fFillFMDSectorInfo(kFALSE),
  fFillFMDChannelInfo(kFALSE),
  //fFillCorrectedFMDInfo(kTRUE),
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
  fV0OpenCuts(0x0),
  fV0StrongCuts(0x0),
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
  fFillFMDSectorInfo(kFALSE),
  fFillFMDChannelInfo(kFALSE),
  //fFillCorrectedFMDInfo(kTRUE),
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
  fV0OpenCuts(0x0),
  fV0StrongCuts(0x0),
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
  //DefineInput(2,AliAODForwardMult::Class());
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
  OpenFile(2);
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
void AliAnalysisTaskReducedTree::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //
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
/*
  cout<<"get AOD"<<endl;
  
  AliAODEvent* aodEvent = AliForwardUtil::GetAODEvent(this);
  if (!aodEvent) return;
  cout<<"got AOD"<<endl;

  TObject* obj = aodEvent->FindListObject("Forward");  
  if (!obj) return;
  cout<<"got AOD forward"<<endl;

  AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(obj);

   //if (!aodForward->CheckEvent(mask,ipZmin,ipZmax,cMin,cMax)) return 0;

  Double_t ret = 0;
  const TH2D& d2Ndetadphi = aodForward->GetHistogram();


  cout<<d2Ndetadphi.GetXaxis()->GetNbins()<<endl;
*/

  // FMD corrections
//  AliAODEvent* aodEvent;
//  TObject* obj;
//  AliAODForwardMult* aodForward=NULL;
//
//  if(fFillCorrectedFMDInfo){
//  aodEvent = AliForwardUtil::GetAODEvent(this);
//  if (!aodEvent) return;
//
//  obj = aodEvent->FindListObject("Forward");  
//  if (!obj) return;}
//
//  aodForward = static_cast<AliAODForwardMult*>(obj);
//
//   //if (!aodForward->CheckE	vent(fTriggerMask,ipZmin,ipZmax,cMin,cMax)) return 0;
//
//  const TH2D& d2Ndetadphi = aodForward->GetHistogram();
  
  
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
  //if(fFillCorrectedFMDInfo) FillCorrectedFMDInfo(d2Ndetadphi);   //Fill corrected FMD info
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
  
  AliESDEvent* esdEvent = 0x0;
  if(isESD) esdEvent = static_cast<AliESDEvent*>(event);
  AliAODEvent* aodEvent = 0x0;
  if(isAOD) aodEvent = static_cast<AliAODEvent*>(event);
  
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
  fReducedEvent->fEventType   = event->GetEventType();
  fReducedEvent->fTriggerMask = event->GetTriggerMask();
  fReducedEvent->fIsPhysicsSelection = (isSelected!=0 ? kTRUE : kFALSE);
  fReducedEvent->fIsSPDPileup = event->IsPileupFromSPD(3,0.8,3.,2.,5.);
  fReducedEvent->fIsSPDPileupMultBins = event->IsPileupFromSPDInMultBins();
  AliVVertex* eventVtx = 0x0;
  if(isESD) eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTracks());
  if(isAOD) eventVtx = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertex());
  if(eventVtx) {
    fReducedEvent->fVtx[0] = (isESD ? ((AliESDVertex*)eventVtx)->GetX() : ((AliAODVertex*)eventVtx)->GetX());
    fReducedEvent->fVtx[1] = (isESD ? ((AliESDVertex*)eventVtx)->GetY() : ((AliAODVertex*)eventVtx)->GetY());
    fReducedEvent->fVtx[2] = (isESD ? ((AliESDVertex*)eventVtx)->GetZ() : ((AliAODVertex*)eventVtx)->GetZ());
    fReducedEvent->fNVtxContributors = eventVtx->GetNContributors();
  }
  if(isESD) {
    eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTPC());
    fReducedEvent->fEventNumberInFile = esdEvent->GetEventNumberInFile();
    fReducedEvent->fL0TriggerInputs = esdEvent->GetHeader()->GetL0TriggerInputs();
    fReducedEvent->fL1TriggerInputs = esdEvent->GetHeader()->GetL1TriggerInputs();
    fReducedEvent->fL2TriggerInputs = esdEvent->GetHeader()->GetL2TriggerInputs();
    fReducedEvent->fIRIntClosestIntMap[0] = esdEvent->GetHeader()->GetIRInt1ClosestInteractionMap();
    fReducedEvent->fIRIntClosestIntMap[1] = esdEvent->GetHeader()->GetIRInt2ClosestInteractionMap();
    if(eventVtx) {
      fReducedEvent->fVtxTPC[0] = ((AliESDVertex*)eventVtx)->GetX();
      fReducedEvent->fVtxTPC[1] = ((AliESDVertex*)eventVtx)->GetY();
      fReducedEvent->fVtxTPC[2] = ((AliESDVertex*)eventVtx)->GetZ();
      fReducedEvent->fNVtxTPCContributors = eventVtx->GetNContributors();
    }
    fReducedEvent->fTimeStamp     = esdEvent->GetTimeStamp();
    fReducedEvent->fNpileupSPD    = esdEvent->GetNumberOfPileupVerticesSPD();
    fReducedEvent->fNpileupTracks = esdEvent->GetNumberOfPileupVerticesTracks();
    fReducedEvent->fNPMDtracks    = esdEvent->GetNumberOfPmdTracks();
    fReducedEvent->fNTRDtracks    = esdEvent->GetNumberOfTrdTracks();
    fReducedEvent->fNTRDtracklets = esdEvent->GetNumberOfTrdTracklets();
    
    AliESDZDC* zdc = esdEvent->GetESDZDC();
    if(zdc) {
      for(Int_t i=0; i<5; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZN1TowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZN2TowerEnergy()[i-5];
      for(Int_t i=0; i<5; ++i)  fReducedEvent->fZDCpEnergy[i]   = zdc->GetZP1TowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  fReducedEvent->fZDCpEnergy[i]   = zdc->GetZP2TowerEnergy()[i-5];
    }
  }
  if(isAOD) {
    AliAODHeader * header = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
    if(!header) AliFatal("Not a standard AOD");


    fReducedEvent->fIRIntClosestIntMap[0] = header->GetIRInt1ClosestInteractionMap();
    fReducedEvent->fIRIntClosestIntMap[1] = header->GetIRInt2ClosestInteractionMap();
    fReducedEvent->fEventNumberInFile = header->GetEventNumberESDFile();
    fReducedEvent->fL0TriggerInputs = header->GetL0TriggerInputs();
    fReducedEvent->fL1TriggerInputs = header->GetL1TriggerInputs();
    fReducedEvent->fL2TriggerInputs = header->GetL2TriggerInputs();
    fReducedEvent->fTimeStamp     = 0;
    fReducedEvent->fNpileupSPD    = aodEvent->GetNumberOfPileupVerticesSPD();
    fReducedEvent->fNpileupTracks = aodEvent->GetNumberOfPileupVerticesTracks();
    fReducedEvent->fNPMDtracks    = aodEvent->GetNPmdClusters();
    fReducedEvent->fNTRDtracks    = 0;
    fReducedEvent->fNTRDtracklets = 0;
    
    AliAODZDC* zdc = aodEvent->GetZDCData();
    if(zdc) {
      for(Int_t i=0; i<5; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZNATowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  fReducedEvent->fZDCnEnergy[i]   = zdc->GetZNCTowerEnergy()[i-5];
      for(Int_t i=0; i<5; ++i)  fReducedEvent->fZDCpEnergy[i]   = zdc->GetZPATowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  fReducedEvent->fZDCpEnergy[i]   = zdc->GetZPCTowerEnergy()[i-5];
    }
  }
  
  // Fill TZERO information
  if(isESD) {
    const AliESDTZERO* tzero = esdEvent->GetESDTZERO();
    if(tzero) {
      fReducedEvent->fT0start = tzero->GetT0();
      fReducedEvent->fT0zVertex = tzero->GetT0zVertex();
      for(Int_t i = 0;i<24;i++)
        fReducedEvent->fT0amplitude[i] = tzero->GetT0amplitude()[i];
      for(Int_t i = 0;i<3;i++)
        fReducedEvent->fT0TOF[i] = tzero->GetT0TOF()[i];
      for(Int_t i = 0;i<3;i++)
        fReducedEvent->fT0TOFbest[i] = tzero->GetT0TOFbest()[i];
      fReducedEvent->fT0pileup = tzero->GetPileupFlag();
      fReducedEvent->fT0sattelite = tzero->GetSatellite();
    }
  }
  if(isAOD) {
    AliAODTZERO* tzero = aodEvent->GetTZEROData();
    if(tzero) {
      fReducedEvent->fT0start = -999.;   // not available
      fReducedEvent->fT0zVertex = tzero->GetT0zVertex();
      for(Int_t i = 0;i<26;i++)
        fReducedEvent->fT0amplitude[i] = tzero->GetAmp(i);
      for(Int_t i = 0;i<3;i++)
        fReducedEvent->fT0TOF[i] = tzero->GetT0TOF()[i];
      for(Int_t i = 0;i<3;i++)
        fReducedEvent->fT0TOFbest[i] = tzero->GetT0TOFbest()[i];
      fReducedEvent->fT0pileup = tzero->GetPileupFlag();
      fReducedEvent->fT0sattelite = tzero->GetSatellite();
    }
  }
  
  if(fFillFMDChannelInfo&&isESD) fReducedEvent->fIsFMDReduced = kFALSE;
  if((fFillFMDSectorInfo||fFillFMDChannelInfo)&&isESD) FillFMDInfo();
  
  AliCentrality *centrality = event->GetCentrality();
  if(centrality) {
    fReducedEvent->fCentrality[0] = centrality->GetCentralityPercentile("V0M");
    fReducedEvent->fCentrality[1] = centrality->GetCentralityPercentile("CL1");
    fReducedEvent->fCentrality[2] = centrality->GetCentralityPercentile("TRK");
    fReducedEvent->fCentrality[3] = centrality->GetCentralityPercentile("ZEMvsZDC");    
    fReducedEvent->fCentQuality   = centrality->GetQuality();
  }
  
  //cout << "event vtxZ/cent: " << fReducedEvent->fVtx[2] << "/" << fReducedEvent->fCentrality[0] << endl;
  
  fReducedEvent->fNtracks[0] = event->GetNumberOfTracks();
  fReducedEvent->fSPDntracklets = GetSPDTrackletMultiplicity(event, -1.0, 1.0);
  for(Int_t ieta=0; ieta<32; ++ieta)
    fReducedEvent->fSPDntrackletsEta[ieta] = GetSPDTrackletMultiplicity(event, -1.6+0.1*ieta, -1.6+0.1*(ieta+1));
  
  AliVVZERO* vzero = event->GetVZEROData();
  for(Int_t i=0;i<64;++i) 
    fReducedEvent->fVZEROMult[i] = vzero->GetMultiplicity(i);  
  
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
    reducedCluster->fM20     = cluster->GetM20();
    reducedCluster->fM02     = cluster->GetM02();
    reducedCluster->fDispersion = cluster->GetDispersion();
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
void AliAnalysisTaskReducedTree::FillFMDInfo()
{
  //
  // fill reduced FMD information
  //
  AliVEvent* event = InputEvent();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());

  if(!isESD) return;

  AliESDEvent* esdEvent = 0x0;
  if(isESD) esdEvent = static_cast<AliESDEvent*>(event);

  AliESDFMD* esdFmd = esdEvent->GetFMDData();
  //const AliFMDFloatMap multMap = esdFmd->MultiplicityMap();
  //const AliFMDFloatMap etaMap  = esdFmd->EtaMap();

  //esdFmd->Print();
  Int_t nFMD=0;
  Int_t id=-1;
  Int_t maxDet=3;
  Int_t maxRing=2;
  Int_t maxSector;
  Int_t maxStrip;
  Float_t m=0.0;
  Double_t phi;
  Char_t ring;
  Float_t fmdMult;
  Float_t msum=0;
  UShort_t fmdDet=0;
  Int_t phiBin=0;
    
  for(UShort_t det = 1; det <= maxDet; ++det) {
    (det == 1 ? maxRing=1 : maxRing=2);
    for(UShort_t ir = 0; ir < maxRing; ++ir) {
      ring = (ir == 0 ? 'I' : 'O');
      (ir == 0 ? maxSector=20 : maxSector=40);
      (ir == 0 ? maxStrip=512 : maxStrip=256);
      nFMD=-1;
      for(UShort_t sec = 0; sec < maxSector; ++sec) {
        phi  =  esdFmd->Phi(det, ring, sec, 0)/180.*TMath::Pi();
        phiBin = Int_t (phi/2/TMath::Pi()*maxSector);
	fmdMult = 0;
        for(UShort_t str = 0; str < maxStrip; ++str) {
          ++id;
          m  =  esdFmd->Multiplicity(det, ring, sec, str);
	  //cout << "det/ir/sec/str/m :: " << det << "/" << ir << "/" << sec << "/" << str << "/" << m << endl;
          if(fFillFMDChannelInfo) 
	    if(m<1.e-6) continue;
          if(m ==  AliESDFMD::kInvalidMult) m=0;
          fmdMult += m;
          msum+=m;
          ++nFMD;
            
	  if(fFillFMDChannelInfo){
            if(m>15.) m=15.;
            m = UShort_t(m*4369+0.5);
            TClonesArray& fmd = *(fReducedEvent->GetFMD(fmdDet));
            AliReducedFMD   *reducedFMD=new(fmd[nFMD]) AliReducedFMD();
	    fReducedEvent->fNFMDchannels[fmdDet] += 1;
            reducedFMD->fMultiplicity     =  m;    
            //reducedFMD->fEta              =  esdFmd->Eta(det, ring, 0, str);
            reducedFMD->fId               =  id;
          }  
        }  // end loop over strips
        
        if(fFillFMDSectorInfo) {
          TClonesArray& fmd = *(fReducedEvent->GetFMD(fmdDet));
          AliReducedFMD   *reducedFMD=new(fmd[phiBin]) AliReducedFMD();
          reducedFMD->fMultiplicity     = fmdMult;
	  fReducedEvent->fNFMDchannels[fmdDet] += 1;
          //cout<<sec<<"  "<<fmdMult<<endl;
          fmdMult=0.0;
        }
      }  // end loop over sectors      
      
      fReducedEvent->fFMDtotalMult[fmdDet] = msum;
      msum=0.0;
      ++fmdDet;
      id=-1;
    }  // end loop over rings
  } // end loop over detectors
}

////_________________________________________________________________________________
//void AliAnalysisTaskReducedTree::FillCorrectedFMDInfo(const TH2D& fmdhist) 
//{
//
//
//Int_t nEta = fmdhist.GetXaxis()->GetNbins();
//Int_t nPhi = fmdhist.GetYaxis()->GetNbins();
////Int_t nBins= fmdhist.GetNbins();
//Float_t eta=0.0;
//Float_t phi=0.0;
//Float_t mult=0.0;
//
//
//Int_t nFMD=0;
////fReducedEvent->fNCorFmdChannels = 0;
//for (Int_t e = 1; e <= nEta; e++) {
//    eta = fmdhist.GetXaxis()->GetBinCenter(e);
//    for (Int_t p = 1; p <= nPhi; p++) {
//         phi = fmdhist.GetYaxis()->GetBinCenter(p);
//         mult = fmdhist.GetBinContent(e, p);
//  	 //TClonesArray& Corfmd = *(fReducedEvent->fCorFMD);
//  	 //AliReducedFMD *reducedCorFMD=new(Corfmd[nFMD]) AliReducedCorFMD();
//    std::cout<<mult<<"  "<<eta<<"  "<<phi<<std::endl;
//	}
//
//	//fReducedEvent->fNCorFmdChannels += 1;
//	nFMD += 1;
////reducedFMD->fCorMultiplicity = mult;
////reducedFMD->fCorEta          = eta;
////reducedFMD->fCorPhi          = phi;
//
//
//}
//
//
//}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTree::FillTrackInfo() 
{
  //
  // fill reduced track information
  //
  AliVEvent* event = InputEvent();
  Bool_t isESD = (event->IsA()==AliESDEvent::Class());
  Bool_t isAOD = (event->IsA()==AliAODEvent::Class());

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  AliPIDResponse* pidResponse = inputHandler->GetPIDResponse();
  
  // find all the tracks which belong to a V0 stored in the reduced event
  UShort_t trackIdsV0[4][20000]={{0}};
  UShort_t trackIdsPureV0[4][20000]={{0}};
  Int_t nV0LegsTagged[4] = {0}; Int_t nPureV0LegsTagged[4] = {0};
  Bool_t leg1Found[4]; Bool_t leg2Found[4];
  for(Int_t iv0=0;iv0<fReducedEvent->fNV0candidates[1];++iv0) {
    AliReducedPair* pair = fReducedEvent->GetV0Pair(iv0);
    if(!pair) continue;
    Int_t pairId = 0; Bool_t isPureV0 = kFALSE;
    if(pair->fCandidateId==AliReducedPair::kGammaConv) {
      pairId=0;
      if(pair->IsPureV0Gamma()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPair::kK0sToPiPi) {
      pairId=1;
      if(pair->IsPureV0K0s()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPair::kLambda0ToPPi) {
      pairId=2;
      if(pair->IsPureV0Lambda()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPair::kALambda0ToPPi) {
      pairId=3;
      if(pair->IsPureV0ALambda()) isPureV0 = kTRUE;
    }
    
    leg1Found[pairId] = kFALSE; leg2Found[pairId] = kFALSE;
    for(Int_t it=0;it<nV0LegsTagged[pairId];++it) {
      if(trackIdsV0[pairId][it]==pair->fLegIds[0]) leg1Found[pairId]=kTRUE;
      if(trackIdsV0[pairId][it]==pair->fLegIds[1]) leg2Found[pairId]=kTRUE;
    }
    // if the legs of this V0 were not already stored then add them now to the list
    if(!leg1Found[pairId]) {trackIdsV0[pairId][nV0LegsTagged[pairId]] = pair->fLegIds[0]; ++nV0LegsTagged[pairId];}
    if(!leg2Found[pairId]) {trackIdsV0[pairId][nV0LegsTagged[pairId]] = pair->fLegIds[1]; ++nV0LegsTagged[pairId];}
    
    if(isPureV0) {
      leg1Found[pairId] = kFALSE; leg2Found[pairId] = kFALSE;
      for(Int_t it=0;it<nPureV0LegsTagged[pairId];++it) {
        if(trackIdsPureV0[pairId][it]==pair->fLegIds[0]) leg1Found[pairId]=kTRUE;
        if(trackIdsPureV0[pairId][it]==pair->fLegIds[1]) leg2Found[pairId]=kTRUE;
      }
      // if the legs of this pure V0 were not already stored then add them now to the list
      if(!leg1Found[pairId]) {trackIdsPureV0[pairId][nPureV0LegsTagged[pairId]] = pair->fLegIds[0]; ++nPureV0LegsTagged[pairId];}
      if(!leg2Found[pairId]) {trackIdsPureV0[pairId][nPureV0LegsTagged[pairId]] = pair->fLegIds[1]; ++nPureV0LegsTagged[pairId];}
    }
  }
  
  // find all the tracks which belong to a stored dielectron pair
  UShort_t trackIdsDiele[20000]={0};
  Int_t nDieleLegsTagged = 0;
  for(Int_t idie=0;idie<fReducedEvent->NDielectrons();++idie) {
    AliReducedPair* pair = fReducedEvent->GetDielectronPair(idie);
    leg1Found[0]=kFALSE; leg2Found[0]=kFALSE;
    for(Int_t it=0; it<nDieleLegsTagged; ++it) {
      if(trackIdsDiele[it]==pair->fLegIds[0]) leg1Found[0]=kTRUE;
      if(trackIdsDiele[it]==pair->fLegIds[1]) leg2Found[0]=kTRUE;
    }
    // if the legs of this dielectron were not already stored then add them now to the list
    if(!leg1Found[0]) {trackIdsDiele[nDieleLegsTagged] = pair->fLegIds[0]; ++nDieleLegsTagged;}
    if(!leg2Found[0]) {trackIdsDiele[nDieleLegsTagged] = pair->fLegIds[1]; ++nDieleLegsTagged;}
  }
    
  AliESDtrack* esdTrack=0;
  AliAODTrack* aodTrack=0;
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t trackId = 0; 
  Bool_t usedForV0[4] = {kFALSE}; 
  Bool_t usedForPureV0[4] = {kFALSE};
  Bool_t usedForV0Or = kFALSE;
  Bool_t usedForDielectron = kFALSE;
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
          //cout << "track " << trackId << " used for V0 type " << i << endl;
          break;
        }
      }
      usedForV0Or = usedForV0Or || usedForV0[i];
      usedForPureV0[i] = kFALSE;
      for(Int_t ii=0; ii<nPureV0LegsTagged[i]; ++ii) {
        if(UShort_t(trackId)==trackIdsPureV0[i][ii]) {
          usedForPureV0[i] = kTRUE;
          //cout << "track " << trackId << " used for pure V0 type " << i << endl;
          break;
        }
      }
    }
    // check whether this track belongs to a dielectron stored in the reduced event
    usedForDielectron = kFALSE;
    for(Int_t ii=0; ii<nDieleLegsTagged; ++ii) {
      if(UShort_t(trackId)==trackIdsDiele[ii]) {
        usedForDielectron = kTRUE;
        break;
      }
    }
    
    ULong_t status = (isESD ? esdTrack->GetStatus() : aodTrack->GetStatus());
    //cout << "TRACK" << endl;
    for(Int_t ibit=0; ibit<32; ++ibit) {
      if(status & (ULong_t(1)<<ibit)) {
        //cout << "bit " << ibit << endl;
        fReducedEvent->fNtracksPerTrackingFlag[ibit] += 1;
      }
    }
    
    //apply track cuts
    if(!usedForV0Or && !usedForDielectron && fTrackFilter && !fTrackFilter->IsSelected(particle)) continue;
    //cout << "storing track " << trackId << endl;
    
    TClonesArray& tracks = *(fReducedEvent->fTracks);
    AliReducedTrack *reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliReducedTrack();
        
    Double_t values[AliDielectronVarManager::kNMaxValues];
    AliDielectronVarManager::Fill(particle, values);
    reducedParticle->fStatus        = status;//(ULong_t)values[AliDielectronVarManager::kTrackStatus];
    reducedParticle->fGlobalPhi     = values[AliDielectronVarManager::kPhi];
    reducedParticle->fGlobalPt      = values[AliDielectronVarManager::kPt]*values[AliDielectronVarManager::kCharge];
    reducedParticle->fGlobalEta     = values[AliDielectronVarManager::kEta];    
    reducedParticle->fMomentumInner = values[AliDielectronVarManager::kPIn];
    reducedParticle->fDCA[0]        = values[AliDielectronVarManager::kImpactParXY];
    reducedParticle->fDCA[1]        = values[AliDielectronVarManager::kImpactParZ];
    reducedParticle->fTrackLength   = values[AliDielectronVarManager::kTrackLength];
    
    reducedParticle->fITSclusterMap = (UChar_t)values[AliDielectronVarManager::kITSclusterMap];
    reducedParticle->fITSsignal     = values[AliDielectronVarManager::kITSsignal];
    reducedParticle->fITSnSig[0]    = values[AliDielectronVarManager::kITSnSigmaEle];
    reducedParticle->fITSnSig[1]    = values[AliDielectronVarManager::kITSnSigmaPio];
    reducedParticle->fITSnSig[2]    = values[AliDielectronVarManager::kITSnSigmaKao];
    reducedParticle->fITSnSig[3]    = values[AliDielectronVarManager::kITSnSigmaPro];
    reducedParticle->fITSchi2       = values[AliDielectronVarManager::kITSchi2Cl];
    
    reducedParticle->fTPCNcls      = (UChar_t)values[AliDielectronVarManager::kNclsTPC];
    reducedParticle->fTPCNclsF     = (UChar_t)values[AliDielectronVarManager::kNFclsTPC];
    reducedParticle->fTPCNclsIter1 = (UChar_t)values[AliDielectronVarManager::kNclsTPCiter1];
    reducedParticle->fTPCsignal    = values[AliDielectronVarManager::kTPCsignal];
    reducedParticle->fTPCsignalN   = values[AliDielectronVarManager::kTPCsignalN];
    reducedParticle->fTPCnSig[0]   = values[AliDielectronVarManager::kTPCnSigmaEle];
    reducedParticle->fTPCnSig[1]   = values[AliDielectronVarManager::kTPCnSigmaPio];
    reducedParticle->fTPCnSig[2]   = values[AliDielectronVarManager::kTPCnSigmaKao];
    reducedParticle->fTPCnSig[3]   = values[AliDielectronVarManager::kTPCnSigmaPro];
    reducedParticle->fTPCClusterMap = EncodeTPCClusterMap(particle, isAOD);
    reducedParticle->fTPCchi2       = values[AliDielectronVarManager::kTPCchi2Cl];
    
    reducedParticle->fTOFbeta      = values[AliDielectronVarManager::kTOFbeta];
    reducedParticle->fTOFnSig[0]   = values[AliDielectronVarManager::kTOFnSigmaEle];
    reducedParticle->fTOFnSig[1]   = values[AliDielectronVarManager::kTOFnSigmaPio];
    reducedParticle->fTOFnSig[2]   = values[AliDielectronVarManager::kTOFnSigmaKao];
    reducedParticle->fTOFnSig[3]   = values[AliDielectronVarManager::kTOFnSigmaPro];
    
    Double_t trdProbab[AliPID::kSPECIES]={0.0};
        
    if(fFlowTrackFilter) {
      // switch on the first bit if this particle should be used for the event plane
      if(fFlowTrackFilter->IsSelected(particle)) reducedParticle->fFlags |= (UShort_t(1)<<0);
    }
    for(Int_t iV0type=0;iV0type<4;++iV0type) {
      if(usedForV0[iV0type]) reducedParticle->fFlags |= (UShort_t(1)<<(iV0type+1));
      if(usedForPureV0[iV0type]) reducedParticle->fFlags |= (UShort_t(1)<<(iV0type+8));
    }
    
    if(isESD){
      //AliESDtrack *track=static_cast<AliESDtrack*>(particle);
      reducedParticle->fTrackId          = (UShort_t)esdTrack->GetID();
      reducedParticle->fTPCCrossedRows   = (UChar_t)esdTrack->GetTPCCrossedRows();
      //reducedParticle->fTPCClusterMap    = EncodeTPCClusterMap(esdTrack);
      const AliExternalTrackParam* tpcInner = esdTrack->GetInnerParam();
      reducedParticle->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
      reducedParticle->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
      reducedParticle->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
      
      reducedParticle->fTOFdeltaBC    = esdTrack->GetTOFDeltaBC();
      
      reducedParticle->fTRDntracklets[0] = esdTrack->GetTRDntracklets();
      reducedParticle->fTRDntracklets[1] = esdTrack->GetTRDntrackletsPID();
      pidResponse->ComputeTRDProbability(esdTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ1D);
      reducedParticle->fTRDpid[0]    = trdProbab[AliPID::kElectron];
      reducedParticle->fTRDpid[1]    = trdProbab[AliPID::kPion];
      pidResponse->ComputeTRDProbability(esdTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ2D);
      reducedParticle->fTRDpidLQ2D[0]    = trdProbab[AliPID::kElectron];
      reducedParticle->fTRDpidLQ2D[1]    = trdProbab[AliPID::kPion];
            
      for(Int_t idx=0; idx<3; ++idx) if(esdTrack->GetKinkIndex(idx)>0) reducedParticle->fFlags |= (1<<(5+idx));
      if(esdTrack->IsEMCAL()) reducedParticle->fCaloClusterId = esdTrack->GetEMCALcluster();
      if(esdTrack->IsPHOS()) reducedParticle->fCaloClusterId = esdTrack->GetPHOScluster();
    }
    if(isAOD) {
      //AliAODTrack *track=static_cast<AliAODTrack*>(particle);
      const AliExternalTrackParam* tpcInner = aodTrack->GetInnerParam();
      reducedParticle->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
      reducedParticle->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
      reducedParticle->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
      
      reducedParticle->fTrackId = aodTrack->GetID(); 
      reducedParticle->fITSsignal = aodTrack->GetITSsignal();
      if(pidResponse) {
        reducedParticle->fITSnSig[0]    = pidResponse->NumberOfSigmasITS(aodTrack,AliPID::kElectron);
        reducedParticle->fITSnSig[1]    = pidResponse->NumberOfSigmasITS(aodTrack,AliPID::kPion);
        reducedParticle->fITSnSig[2]    = pidResponse->NumberOfSigmasITS(aodTrack,AliPID::kKaon);
        reducedParticle->fITSnSig[3]    = pidResponse->NumberOfSigmasITS(aodTrack,AliPID::kProton);
      }
      reducedParticle->fTRDntracklets[0] = aodTrack->GetTRDntrackletsPID();
      reducedParticle->fTRDntracklets[1] = aodTrack->GetTRDntrackletsPID();
      pidResponse->ComputeTRDProbability(aodTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ1D);
      reducedParticle->fTRDpid[0]    = trdProbab[AliPID::kElectron];
      reducedParticle->fTRDpid[1]    = trdProbab[AliPID::kPion];
      pidResponse->ComputeTRDProbability(aodTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ2D);
      reducedParticle->fTRDpidLQ2D[0]    = trdProbab[AliPID::kElectron];
      reducedParticle->fTRDpidLQ2D[1]    = trdProbab[AliPID::kPion];
      
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
      reducedParticle->fLegIds[0]    = (UShort_t)(static_cast<AliVTrack*>(pair->GetFirstDaughterP()))->GetID();
      reducedParticle->fLegIds[1]    = (UShort_t)(static_cast<AliVTrack*>(pair->GetSecondDaughterP()))->GetID();
      reducedParticle->fMass[0]      = values[AliDielectronVarManager::kM];
      reducedParticle->fMass[1]      = -999.;
      reducedParticle->fMass[2]      = -999.;
      reducedParticle->fMass[3]      = -999.;
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
  
  if(!(fFillK0s || fFillLambda || fFillALambda || fFillGammaConversions)) return;
    
  Double_t valuesPos[AliDielectronVarManager::kNMaxValues];
  Double_t valuesNeg[AliDielectronVarManager::kNMaxValues];
  
  if(fV0OpenCuts) {
    fV0OpenCuts->SetEvent(esd);
    fV0OpenCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  if(fV0StrongCuts) {
    fV0StrongCuts->SetEvent(esd);
    fV0StrongCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  
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
    
    pdgV0=0; pdgP=0; pdgN=0;
    Bool_t goodK0s = kTRUE; Bool_t goodLambda = kTRUE; Bool_t goodALambda = kTRUE; Bool_t goodGamma = kTRUE;
    if(fV0OpenCuts) {
      goodK0s = kFALSE; goodLambda = kFALSE; goodALambda = kFALSE; goodGamma = kFALSE;
      Bool_t processV0 = fV0OpenCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
      if(processV0 && TMath::Abs(pdgV0)==310 && TMath::Abs(pdgP)==211 && TMath::Abs(pdgN)==211) {
        goodK0s = kTRUE;
        if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(legPos) || !fK0sPionCuts->IsSelected(legNeg))) goodK0s = kFALSE;
      }
      if(processV0 && pdgV0==3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212)) {
        goodLambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legPos)) goodLambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legNeg)) goodLambda = kFALSE;
      }
      if(processV0 && pdgV0==-3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212)) {
        goodALambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(legNeg)) goodALambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(legPos)) goodALambda = kFALSE;
      }
      if(processV0 && TMath::Abs(pdgV0)==22 && TMath::Abs(pdgP)==11 && TMath::Abs(pdgN)==11) {
        goodGamma = kTRUE;
        if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(legPos) || !fGammaElectronCuts->IsSelected(legNeg))) goodGamma = kFALSE;
      }
      //cout << "open cuts  pdgV0/pdgP/pdgN/processV0 : " << pdgV0 << "/" << pdgP << "/" << pdgN << "/" << processV0 << endl;     
      //cout << "good K0s/Lambda/ALambda/Gamma : " << goodK0s << "/" << goodLambda << "/" << goodALambda << "/" << goodGamma << endl;
    }
    
    Bool_t veryGoodK0s = kFALSE; Bool_t veryGoodLambda = kFALSE; Bool_t veryGoodALambda = kFALSE; Bool_t veryGoodGamma = kFALSE;
    if(fV0StrongCuts && (goodK0s || goodLambda || goodALambda || goodGamma)) {
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0StrongCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
      if(processV0 && goodK0s && TMath::Abs(pdgV0)==310 && TMath::Abs(pdgP)==211 && TMath::Abs(pdgN)==211)
        veryGoodK0s = kTRUE;
      if(processV0 && goodLambda && pdgV0==3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212))
        veryGoodLambda = kTRUE;
      if(processV0 && goodALambda && pdgV0==-3122 && (TMath::Abs(pdgP)==211 || TMath::Abs(pdgP)==2212) && (TMath::Abs(pdgN)==211 || TMath::Abs(pdgN)==2212))
        veryGoodALambda = kTRUE;
      if(processV0 && goodGamma && TMath::Abs(pdgV0)==22 && TMath::Abs(pdgP)==11 && TMath::Abs(pdgN)==11)
        veryGoodGamma = kTRUE;
      //cout << "strong cuts  pdgV0/pdgP/pdgN/processV0 : " << pdgV0 << "/" << pdgP << "/" << pdgN << "/" << processV0 << endl;     
      //cout << "very good K0s/Lambda/ALambda/Gamma : " << veryGoodK0s << "/" << veryGoodLambda << "/" << veryGoodALambda << "/" << veryGoodGamma << endl;
    }
              
    if(!((goodK0s && fFillK0s) || 
         (goodLambda && fFillLambda) || 
         (goodALambda && fFillALambda) || 
         (goodGamma && fFillGammaConversions))) continue;
    
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
      if(veryGoodK0s) goodK0sPair->fMCid |= (UInt_t(1)<<1);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodK0s=kFALSE;}
    if(fFillLambda && goodLambda && lambdaReducedPair->fMass[0]>fLambdaMassRange[0] && lambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodLambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*lambdaReducedPair);
      goodLambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodLambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodLambda) goodLambdaPair->fMCid |= (UInt_t(1)<<2);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodLambda=kFALSE;}
    if(fFillALambda && goodALambda && alambdaReducedPair->fMass[0]>fLambdaMassRange[0] && alambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodALambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*alambdaReducedPair);
      goodALambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodALambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodALambda) goodALambdaPair->fMCid |= (UInt_t(1)<<3);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodALambda = kFALSE;}
    //cout << "gamma mass: " << gammaReducedPair->fMass[0] << endl;
    if(fFillGammaConversions && goodGamma && gammaReducedPair->fMass[0]>fGammaMassRange[0] && gammaReducedPair->fMass[0]<fGammaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPair *goodGammaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPair(*gammaReducedPair);
      goodGammaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodGammaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodGammaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodGammaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodGamma) goodGammaPair->fMCid |= (UInt_t(1)<<4);
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
    reducedPair->fChisquare    = v0->GetChi2V0();
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
    reducedPair->fChisquare = v0Refit.GetChi2();                              
  }
  return reducedPair;
}


//_________________________________________________________________________________
UChar_t AliAnalysisTaskReducedTree::EncodeTPCClusterMap(AliVParticle* track, Bool_t isAOD) {
  //
  // Encode the TPC cluster map into an UChar_t
  // Divide the 159 bits from the bit map into 8 groups of adiacent clusters
  // For each group enable its corresponding bit if in that group there are more clusters compared to
  // a threshold.
  //
  AliESDtrack* esdTrack=0x0;
  AliAODTrack* aodTrack=0x0;
  if(isAOD)
    aodTrack=static_cast<AliAODTrack*>(track);
  else
    esdTrack=static_cast<AliESDtrack*>(track);
  
  const UChar_t threshold=5;
  TBits tpcClusterMap = (isAOD ? aodTrack->GetTPCClusterMap() : esdTrack->GetTPCClusterMap());
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
Int_t AliAnalysisTaskReducedTree::GetSPDTrackletMultiplicity(AliVEvent* event, Float_t lowEta, Float_t highEta) {
  //
  // Count the number of SPD tracklets in a given eta range
  //
  if (!event) return -1;
  
  Int_t nTracklets = 0;
  Int_t nAcc = 0;
  
  if(event->IsA() == AliAODEvent::Class()) {
    AliAODTracklets *tracklets = ((AliAODEvent*)event)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for(Int_t nn=0; nn<nTracklets; ++nn) {
      Double_t theta = tracklets->GetTheta(nn);
      Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
      if(eta < lowEta) continue;
      if(eta > highEta) continue;
      ++nAcc;
    }
  } else if(event->IsA() == AliESDEvent::Class()) {
    nTracklets = ((AliESDEvent*)event)->GetMultiplicity()->GetNumberOfTracklets();
    for(Int_t nn=0; nn<nTracklets; ++nn) {
      Double_t eta = ((AliESDEvent*)event)->GetMultiplicity()->GetEta(nn);
      if(eta < lowEta) continue;
      if(eta > highEta) continue; 
      ++nAcc;
    }
  } else return -1;
  
  return nAcc;
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
