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
#include <TBits.h>

#include <AliAnalysisTaskSE.h>
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
#include <AliAODForwardMult.h>
#include <AliForwardUtil.h>
#include <AliTriggerAnalysis.h>
#include <AliAODTrack.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include <AliVZDC.h>
#include <AliESDv0.h>
#include <AliESDv0Cuts.h>
#include <AliESDv0KineCuts.h>
#include <AliESDFMD.h>
#include <AliVCluster.h>
#include <AliAODTracklets.h>
#include <AliMultiplicity.h>
#include <AliAODTracklets.h>
#include <AliPIDResponse.h>
//#include <AliFlowBayesianPID.h>
#include <AliMCParticle.h>
#include <AliAODMCParticle.h>
#include "AliAnalysisUtils.h"
#include <AliMultSelection.h>
#include <AliMultEstimator.h>
#include <AliCentrality.h>
#include "AliDielectronVarManager.h"
//#include "AliFlowTrackCuts.h"
#include "AliReducedEventInfo.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedFMDInfo.h"
#include "AliReducedEventPlaneInfo.h"
#include "AliAnalysisTaskReducedTreeMaker.h"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskReducedTreeMaker)

//_________________________________________________________________________________
AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker() :
  AliAnalysisTaskSE(),
  fAnalysisUtils(0x0),
  fUseAnalysisUtils(kFALSE),
  fMinVtxContributors(0),
  fMaxVtxZ(100.),
  fCutOnSPDVtxZ(kFALSE),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fTreeWritingOption(kBaseEventsWithBaseTracks),
  fWriteTree(kTRUE),
  fWriteEventsWithNoSelectedTracks(kTRUE),
  fFillTrackInfo(kTRUE),
  fFillV0Info(kTRUE),
  fFillGammaConversions(kTRUE),
  fFillK0s(kTRUE),
  fFillLambda(kTRUE),
  fFillALambda(kTRUE),
  fFillCaloClusterInfo(kTRUE),
  fFillFMDInfo(kFALSE),
  //fFillBayesianPIDInfo(kFALSE),
  fFillEventPlaneInfo(kFALSE),
  fFillMCInfo(kFALSE),
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
  fFMDhist(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fGammaMassRange(),
  fActiveBranches(""),
  fInactiveBranches(""),
  //fAliFlowTrackCuts(0x0),
  //fBayesianResponse(0x0),
  fTreeFile(0x0),
  fTree(0x0),
  fReducedEvent(0x0),
  fUsedVars(0x0),
  fNevents(0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskReducedTreeMaker::AliAnalysisTaskReducedTreeMaker(const char *name, Bool_t writeTree /*=kTRUE*/) :
  AliAnalysisTaskSE(name),
  fAnalysisUtils(0x0),
  fUseAnalysisUtils(kFALSE),
  fMinVtxContributors(0),
  fMaxVtxZ(100.),
  fCutOnSPDVtxZ(kFALSE),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kAny),
  fRejectPileup(kFALSE),
  fTreeWritingOption(kBaseEventsWithBaseTracks),
  fWriteTree(writeTree),
  fWriteEventsWithNoSelectedTracks(kTRUE),
  fFillTrackInfo(kTRUE),
  fFillV0Info(kTRUE),
  fFillGammaConversions(kTRUE),
  fFillK0s(kTRUE),
  fFillLambda(kTRUE),
  fFillALambda(kTRUE),
  fFillCaloClusterInfo(kTRUE),
  fFillFMDInfo(kFALSE),
  //fFillBayesianPIDInfo(kFALSE),
  fFillEventPlaneInfo(kFALSE),
  fFillMCInfo(kFALSE),
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
  fFMDhist(0x0),
  fK0sMassRange(),
  fLambdaMassRange(),
  fGammaMassRange(),
  fActiveBranches(""),
  fInactiveBranches(""),
  //fAliFlowTrackCuts(0x0),
  //fBayesianResponse(0x0),
  fTreeFile(0x0),
  fTree(0x0),
  fReducedEvent(0x0),
  fUsedVars(0x0),
  fNevents(0)
{
  //
  // Constructor
  //
  fK0sMassRange[0] = 0.4; fK0sMassRange[1] = 0.6;
  fLambdaMassRange[0] = 1.08; fLambdaMassRange[1] = 1.15;
  fGammaMassRange[0] = 0.0; fGammaMassRange[1] = 0.1;

  //fAliFlowTrackCuts = new AliFlowTrackCuts();
 
  DefineInput(0,TChain::Class());
  //DefineInput(2,AliAODForwardMult::Class());
  DefineOutput(1, AliReducedBaseEvent::Class());   // reduced information tree
  if(writeTree)
    DefineOutput(2, TTree::Class());   // reduced information tree
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //
  if(fUseAnalysisUtils) fAnalysisUtils = new AliAnalysisUtils();
  if (fTree) return; //already initialised
  
  if(fWriteTree) {
    OpenFile(2);
    fTree = new TTree("DstTree","Reduced ESD information");
  }
  
  switch(fTreeWritingOption) {
     case kBaseEventsWithBaseTracks:
        fReducedEvent = new AliReducedBaseEvent("DstEvent", AliReducedBaseEvent::kUseBaseTracks);
        break;
     case kBaseEventsWithFullTracks:
        fReducedEvent = new AliReducedBaseEvent("DstEvent", AliReducedBaseEvent::kUseReducedTracks);
        break;
     case kFullEventsWithBaseTracks:
        fReducedEvent = new AliReducedEventInfo("DstEvent", AliReducedBaseEvent::kUseBaseTracks);   
        break;
     case kFullEventsWithFullTracks:
        fReducedEvent = new AliReducedEventInfo("DstEvent", AliReducedBaseEvent::kUseReducedTracks);   
        break;
     default:
        break;
  };
 
  if(fWriteTree)
    fTree->Branch("Event",&fReducedEvent,16000,99);

  // if user set active branches
  TObjArray* aractive=fActiveBranches.Tokenize(";");
  if(aractive->GetEntries()>0) {fTree->SetBranchStatus("*", 0);}
  for(Int_t i=0; i<aractive->GetEntries(); i++){
    fTree->SetBranchStatus(aractive->At(i)->GetName(), 1);
  }
  
  // if user set inactive branches
  TObjArray* arinactive=fInactiveBranches.Tokenize(";");
  for(Int_t i=0; i<arinactive->GetEntries(); i++){
    fTree->SetBranchStatus(arinactive->At(i)->GetName(), 0);
  }
 
  // if MC info is not requested, then set the respective branches off
  if(!fFillMCInfo) {
    fTree->SetBranchStatus("fTracks.fMC*", 0); 
  }
  if(!fFillEventPlaneInfo) {
    fTree->SetBranchStatus("fEventPlane.*", 0);   
  }
 
  /*if(fFillBayesianPIDInfo) {
    fBayesianResponse = new AliFlowBayesianPID();
    fBayesianResponse->SetNewTrackParam();
  }*/
  
  // enable all variables in the VarManager
  fUsedVars = new TBits(AliDielectronVarManager::kNMaxValues);
  for(Int_t i=0;i<AliDielectronVarManager::kNMaxValues;++i) fUsedVars->SetBitNumber(i,kTRUE);
  
  
  PostData(1, fReducedEvent);
  if(fWriteTree)
    PostData(2, fTree);
  //if(fFillFriendInfo) PostData(3, fFriendTree);
  //PostData(2, fFriendTree);
  //PostData(1, fTree);
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::UserExec(Option_t *option)
{
  //
  // Main loop. Called for every event
  //  
  option = option;
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  fNevents++;
  
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

  fReducedEvent->ClearEvent();
  
  if(isSelected==0) {
    cout << "AliAnalysisTaskReducedTreeMaker::UserExec() event is not selected" << endl;
    //PostData(1, fReducedEvent);
    return;
  }

  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
  }
  
  if(fFillMCInfo) {
    Bool_t hasMC=AliDielectronMC::Instance()->HasMC();
    if(hasMC) AliDielectronMC::Instance()->ConnectMCEvent();
  }

  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );
  
  //Fill event wise information
  FillEventInfo();
  
  // NOTE: It is important that FillV0PairInfo() is called before FillTrackInfo()
  if(fFillMCInfo) FillMCTruthInfo();
  if(fFillV0Info && isESD) FillV0PairInfo();
  if(fFillTrackInfo) FillTrackInfo();
 
  if(fWriteTree) {
    if(fWriteEventsWithNoSelectedTracks) fTree->Fill();
    if(!fWriteEventsWithNoSelectedTracks && fReducedEvent->fNtracks[1]>0) fTree->Fill();
  }
        
  // if there are candidate pairs, add the information to the reduced tree
  //if(fFillFriendInfo) PostData(3, fFriendTree);
  PostData(1, fReducedEvent);
  //PostData(2, fFriendTree);
  if(fWriteTree)
    PostData(2, fTree);
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillEventInfo() 
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
  
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::SetEvent(event);

  if(fUseAnalysisUtils) {
    if(fAnalysisUtils->IsVertexSelected2013pA(event))  // 2013 p-Pb event selection    
      fReducedEvent->fEventTag |= (ULong64_t(1)<<0);
    fAnalysisUtils->SetMinPlpContribMV(5); fAnalysisUtils->SetMaxPlpChi2MV(5.);
    fAnalysisUtils->SetCheckPlpFromDifferentBCMV(kTRUE);
    fAnalysisUtils->SetMinWDistMV(15.);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup
      fReducedEvent->fEventTag |= (ULong64_t(1)<<1);
    fAnalysisUtils->SetCheckPlpFromDifferentBCMV(kFALSE);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup, with no BC check
      fReducedEvent->fEventTag |= (ULong64_t(1)<<2);
    fAnalysisUtils->SetCheckPlpFromDifferentBCMV(kTRUE);
    fAnalysisUtils->SetMinWDistMV(10.);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup, with min weighted distance 10
      fReducedEvent->fEventTag |= (ULong64_t(1)<<3);
    fAnalysisUtils->SetMinWDistMV(5.);
    if(fAnalysisUtils->IsPileUpMV(event))              // multi-vertexer pileup, with min weighted distance 5
      fReducedEvent->fEventTag |= (ULong64_t(1)<<4);
  }
  
  if(event->IsPileupFromSPD(3,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<5);
  if(event->IsPileupFromSPD(4,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<6);
  if(event->IsPileupFromSPD(5,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<7);
  if(event->IsPileupFromSPD(6,0.6,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<8);
  if(event->IsPileupFromSPD(3,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<9);
  if(event->IsPileupFromSPD(4,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<10);
  if(event->IsPileupFromSPD(5,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<11);
  if(event->IsPileupFromSPD(6,0.8,3.,2.,5.)) fReducedEvent->fEventTag |= (ULong64_t(1)<<12);
  
  fReducedEvent->fRunNo       = event->GetRunNumber();
  AliVVertex* eventVtx = 0x0;
  if(isESD) eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTracks());
  if(isAOD) eventVtx = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertex());
  if(eventVtx) {
     fReducedEvent->fVtx[0] = (isESD ? ((AliESDVertex*)eventVtx)->GetX() : ((AliAODVertex*)eventVtx)->GetX());
     fReducedEvent->fVtx[1] = (isESD ? ((AliESDVertex*)eventVtx)->GetY() : ((AliAODVertex*)eventVtx)->GetY());
     fReducedEvent->fVtx[2] = (isESD ? ((AliESDVertex*)eventVtx)->GetZ() : ((AliAODVertex*)eventVtx)->GetZ());
     fReducedEvent->fNVtxContributors = eventVtx->GetNContributors();
  }
  
  AliCentrality *centrality = 0x0;
  AliMultSelection* multSelection = 0x0;
  if(event->GetRunNumber()<200000) {
     centrality = event->GetCentrality();
     if(centrality) {
       fReducedEvent->fCentrality[0] = centrality->GetCentralityPercentile("V0M");
       fReducedEvent->fCentrality[1] = centrality->GetCentralityPercentile("CL1");
       fReducedEvent->fCentrality[2] = centrality->GetCentralityPercentile("TRK");
       fReducedEvent->fCentrality[3] = centrality->GetCentralityPercentile("ZEMvsZDC");
       fReducedEvent->fCentrality[4] = centrality->GetCentralityPercentile("V0A");
       fReducedEvent->fCentrality[5] = centrality->GetCentralityPercentile("V0C");
       fReducedEvent->fCentrality[6] = centrality->GetCentralityPercentile("ZNA");
       fReducedEvent->fCentQuality   = centrality->GetQuality();
     }
  }
  else {
     multSelection = (AliMultSelection*)event->FindListObject("MultSelection");
     if(multSelection) {
        fReducedEvent->fCentrality[0] = multSelection->GetMultiplicityPercentile("V0M");
        fReducedEvent->fCentrality[1] = multSelection->GetMultiplicityPercentile("CL1");
        fReducedEvent->fCentrality[2] = multSelection->GetMultiplicityPercentile("TRK");
        fReducedEvent->fCentrality[3] = multSelection->GetMultiplicityPercentile("ZEMvsZDC");
        fReducedEvent->fCentrality[4] = multSelection->GetMultiplicityPercentile("V0A");
        fReducedEvent->fCentrality[5] = multSelection->GetMultiplicityPercentile("V0C");
        fReducedEvent->fCentrality[6] = multSelection->GetMultiplicityPercentile("ZNA");
        fReducedEvent->fCentQuality   = multSelection->GetEvSelCode();
    }
  }
  fReducedEvent->fNtracks[0] = event->GetNumberOfTracks();

  // In case we want to write just basic event information, we stop here
  if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kBaseEventsWithFullTracks) 
     return;
  
  AliReducedEventInfo* eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
  if(!eventInfo) return;
  
  if(multSelection) {
     eventInfo->fMultiplicityEstimatorPercentiles[0] = multSelection->GetMultiplicityPercentile("OnlineV0M");
     eventInfo->fMultiplicityEstimatorPercentiles[1] = multSelection->GetMultiplicityPercentile("OnlineV0A");
     eventInfo->fMultiplicityEstimatorPercentiles[2] = multSelection->GetMultiplicityPercentile("OnlineV0C");
     eventInfo->fMultiplicityEstimatorPercentiles[3] = multSelection->GetMultiplicityPercentile("ADM");
     eventInfo->fMultiplicityEstimatorPercentiles[4] = multSelection->GetMultiplicityPercentile("ADA");
     eventInfo->fMultiplicityEstimatorPercentiles[5] = multSelection->GetMultiplicityPercentile("ADC");
     eventInfo->fMultiplicityEstimatorPercentiles[6] = multSelection->GetMultiplicityPercentile("SPDClusters");
     eventInfo->fMultiplicityEstimatorPercentiles[7] = multSelection->GetMultiplicityPercentile("SPDTracklets");
     eventInfo->fMultiplicityEstimatorPercentiles[8] = multSelection->GetMultiplicityPercentile("RefMult05");
     eventInfo->fMultiplicityEstimatorPercentiles[9] = multSelection->GetMultiplicityPercentile("RefMult08");
     AliMultEstimator* estimator = 0x0;
     estimator = multSelection->GetEstimator("OnlineV0M"); if(estimator) eventInfo->fMultiplicityEstimators[0] = estimator->GetValue();
     estimator = multSelection->GetEstimator("OnlineV0A"); if(estimator) eventInfo->fMultiplicityEstimators[1] = estimator->GetValue();
     estimator = multSelection->GetEstimator("OnlineV0C"); if(estimator) eventInfo->fMultiplicityEstimators[2] = estimator->GetValue();
     estimator = multSelection->GetEstimator("ADM"); if(estimator) eventInfo->fMultiplicityEstimators[3] = estimator->GetValue();
     estimator = multSelection->GetEstimator("ADA"); if(estimator) eventInfo->fMultiplicityEstimators[4] = estimator->GetValue();
     estimator = multSelection->GetEstimator("ADC"); if(estimator) eventInfo->fMultiplicityEstimators[5] = estimator->GetValue();
     estimator = multSelection->GetEstimator("SPDClusters"); if(estimator) eventInfo->fMultiplicityEstimators[6] = estimator->GetValue();
     estimator = multSelection->GetEstimator("SPDTracklets"); if(estimator) eventInfo->fMultiplicityEstimators[7] = estimator->GetValue();
     estimator = multSelection->GetEstimator("RefMult05"); if(estimator) eventInfo->fMultiplicityEstimators[8] = estimator->GetValue();
     estimator = multSelection->GetEstimator("RefMult08"); if(estimator) eventInfo->fMultiplicityEstimators[9] = estimator->GetValue();     
  }
  
  AliVVertex* eventVtxSPD = 0x0;
  if(isESD) eventVtxSPD = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexSPD());
  if(isAOD) eventVtxSPD = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertexSPD());
  if(eventVtxSPD) {
     eventInfo->fVtxSPD[0] = (isESD ? ((AliESDVertex*)eventVtxSPD)->GetX() : ((AliAODVertex*)eventVtxSPD)->GetX());
     eventInfo->fVtxSPD[1] = (isESD ? ((AliESDVertex*)eventVtxSPD)->GetY() : ((AliAODVertex*)eventVtxSPD)->GetY());
     eventInfo->fVtxSPD[2] = (isESD ? ((AliESDVertex*)eventVtxSPD)->GetZ() : ((AliAODVertex*)eventVtxSPD)->GetZ());
     eventInfo->fNVtxSPDContributors = eventVtxSPD->GetNContributors();
  }
  
  // ------------------------------------------------------------------------------------------------------------------
  // Improved cut on the distance between SPD and track vertices 
  // See Francesco Prino's slides during Physics Forum from 5 october 2016, slide 33
  //  based on input from Ruben Shahoyan and Alex Dobrin
  //------------------------------------------------------------------------------------------------------------------
  Bool_t vertexDistanceSelected = kTRUE;
  if(!eventVtx) vertexDistanceSelected = kFALSE;
  if(!eventVtxSPD) vertexDistanceSelected = kFALSE;
  if(vertexDistanceSelected) {
     if(eventVtx->GetNContributors()<2 || eventVtxSPD->GetNContributors()<1) vertexDistanceSelected = kFALSE;   
  }
  if(vertexDistanceSelected) {
     Double_t covTracks[6], covSPD[6];
     eventVtx->GetCovarianceMatrix(covTracks);
     eventVtxSPD->GetCovarianceMatrix(covSPD);
     Double_t dz = eventVtx->GetZ() - eventVtxSPD->GetZ();
     Double_t errTot = TMath::Sqrt(covTracks[5]+covSPD[5]);
     Double_t errTrk = TMath::Sqrt(covTracks[5]);
     Double_t nsigTot = TMath::Abs(dz)/errTot;
     Double_t nsigTrk = TMath::Abs(dz)/errTrk;
     if(TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrk>20) vertexDistanceSelected = kFALSE;
  }
  if(vertexDistanceSelected) fReducedEvent->fEventTag |= (ULong64_t(1)<<13);
  //-------------------------------------------------------------------------------------------------------------------------
  
  eventInfo->fBC          = event->GetBunchCrossNumber();
  eventInfo->fEventType   = event->GetEventType();
  //eventInfo->fTriggerMask = event->GetTriggerMask();
  eventInfo->fTriggerMask = inputHandler->IsEventSelected();
  eventInfo->fIsPhysicsSelection = (isSelected!=0 ? kTRUE : kFALSE);
  eventInfo->fIsSPDPileup = event->IsPileupFromSPD(3,0.8,3.,2.,5.);
  eventInfo->fIsSPDPileupMultBins = event->IsPileupFromSPDInMultBins();
  
  if(isESD) {
    eventInfo->fEventNumberInFile = esdEvent->GetEventNumberInFile();
    eventInfo->fL0TriggerInputs = esdEvent->GetHeader()->GetL0TriggerInputs();
    eventInfo->fL1TriggerInputs = esdEvent->GetHeader()->GetL1TriggerInputs();
    eventInfo->fL2TriggerInputs = esdEvent->GetHeader()->GetL2TriggerInputs();
    eventInfo->fIRIntClosestIntMap[0] = esdEvent->GetHeader()->GetIRInt1ClosestInteractionMap();
    eventInfo->fIRIntClosestIntMap[1] = esdEvent->GetHeader()->GetIRInt2ClosestInteractionMap();
    eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTPC());
    if(eventVtx) {
      eventInfo->fVtxTPC[0] = ((AliESDVertex*)eventVtx)->GetX();
      eventInfo->fVtxTPC[1] = ((AliESDVertex*)eventVtx)->GetY();
      eventInfo->fVtxTPC[2] = ((AliESDVertex*)eventVtx)->GetZ();
      eventInfo->fNVtxTPCContributors = eventVtx->GetNContributors();
    }
    eventInfo->fTimeStamp     = esdEvent->GetTimeStamp();
    eventInfo->fNpileupSPD    = esdEvent->GetNumberOfPileupVerticesSPD();
    eventInfo->fNpileupTracks = esdEvent->GetNumberOfPileupVerticesTracks();
    eventInfo->fNPMDtracks    = esdEvent->GetNumberOfPmdTracks();
    eventInfo->fNTRDtracks    = esdEvent->GetNumberOfTrdTracks();
    eventInfo->fNTRDtracklets = esdEvent->GetNumberOfTrdTracklets();
    eventInfo->fNTPCclusters  = esdEvent->GetNumberOfTPCClusters();
    
    for(Int_t ilayer=0; ilayer<2; ++ilayer)
      eventInfo->fSPDFiredChips[ilayer] = esdEvent->GetMultiplicity()->GetNumberOfFiredChips(ilayer);
    for(Int_t ilayer=0; ilayer<6; ++ilayer)
      eventInfo->fITSClusters[ilayer] = esdEvent->GetMultiplicity()->GetNumberOfITSClusters(ilayer);
    eventInfo->fSPDnSingle = esdEvent->GetMultiplicity()->GetNumberOfSingleClusters();
    
    AliESDZDC* zdc = esdEvent->GetESDZDC();
    if(zdc) {
      eventInfo->fZDCnTotalEnergy[0] = zdc->GetZN2TowerEnergy()[0];
      eventInfo->fZDCnTotalEnergy[1] = zdc->GetZN1TowerEnergy()[0];
      eventInfo->fZDCpTotalEnergy[0] = zdc->GetZP2TowerEnergy()[0];
      eventInfo->fZDCpTotalEnergy[1] = zdc->GetZP1TowerEnergy()[0];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZN1TowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZN2TowerEnergy()[i-5];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZP1TowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZP2TowerEnergy()[i-5];
      
    }
  }
  if(isAOD) {
    eventInfo->fIRIntClosestIntMap[0] = aodEvent->GetHeader()->GetIRInt1ClosestInteractionMap();
    eventInfo->fIRIntClosestIntMap[1] = aodEvent->GetHeader()->GetIRInt2ClosestInteractionMap();
    eventInfo->fEventNumberInFile = aodEvent->GetEventNumberInFile();
    eventInfo->fL0TriggerInputs = aodEvent->GetHeader()->GetL0TriggerInputs();
    eventInfo->fL1TriggerInputs = aodEvent->GetHeader()->GetL1TriggerInputs();
    eventInfo->fL2TriggerInputs = aodEvent->GetHeader()->GetL2TriggerInputs();
    eventInfo->fTimeStamp     = aodEvent->GetTimeStamp();
    eventInfo->fNpileupSPD    = aodEvent->GetNumberOfPileupVerticesSPD();
    eventInfo->fNpileupTracks = aodEvent->GetNumberOfPileupVerticesTracks();
    eventInfo->fNPMDtracks    = aodEvent->GetNPmdClusters();
    eventInfo->fNTRDtracks    = aodEvent->GetNumberOfTrdTracks();
    eventInfo->fNTRDtracklets = 0;
    eventInfo->fNTPCclusters  = aodEvent->GetNumberOfTPCClusters();
    
    eventVtx = const_cast<AliAODVertex*>(aodEvent->GetPrimaryVertexTPC());
    if(eventVtx) {
       eventInfo->fVtxTPC[0] = ((AliAODVertex*)eventVtx)->GetX();
       eventInfo->fVtxTPC[1] = ((AliAODVertex*)eventVtx)->GetY();
       eventInfo->fVtxTPC[2] = ((AliAODVertex*)eventVtx)->GetZ();
       eventInfo->fNVtxTPCContributors = eventVtx->GetNContributors();
    }
    
    for(Int_t ilayer=0; ilayer<2; ++ilayer)
      eventInfo->fSPDFiredChips[ilayer] = aodEvent->GetMultiplicity()->GetNumberOfFiredChips(ilayer);
    for(Int_t ilayer=0; ilayer<6; ++ilayer)
       eventInfo->fITSClusters[ilayer] = aodEvent->GetMultiplicity()->GetNumberOfITSClusters(ilayer);
    
    AliAODZDC* zdc = aodEvent->GetZDCData();
    if(zdc) {
       eventInfo->fZDCnTotalEnergy[0] = zdc->GetZNATowerEnergy()[0];
       eventInfo->fZDCnTotalEnergy[1] = zdc->GetZNCTowerEnergy()[0];
       eventInfo->fZDCpTotalEnergy[0] = zdc->GetZPATowerEnergy()[0];
       eventInfo->fZDCpTotalEnergy[1] = zdc->GetZPCTowerEnergy()[0];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZNATowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCnEnergy[i]   = zdc->GetZNCTowerEnergy()[i-5];
      for(Int_t i=0; i<5; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZPATowerEnergy()[i];
      for(Int_t i=5; i<10; ++i)  eventInfo->fZDCpEnergy[i]   = zdc->GetZPCTowerEnergy()[i-5];
    }
  }
  
  // Fill TZERO information
  if(isESD) {
    const AliESDTZERO* tzero = esdEvent->GetESDTZERO();
    if(tzero) {
      eventInfo->fT0start = tzero->GetT0();
      eventInfo->fT0zVertex = tzero->GetT0zVertex();
      for(Int_t i = 0;i<24;i++)
        eventInfo->fT0amplitude[i] = tzero->GetT0amplitude()[i];
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOF[i] = tzero->GetT0TOF()[i];
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOFbest[i] = tzero->GetT0TOFbest()[i];
      eventInfo->fT0pileup = tzero->GetPileupFlag();
      eventInfo->fT0sattelite = tzero->GetSatellite();
    }
  }
  if(isAOD) {
    AliAODTZERO* tzero = aodEvent->GetTZEROData();
    if(tzero) {
      eventInfo->fT0start = -999.;   // not available
      eventInfo->fT0zVertex = tzero->GetT0zVertex();
      for(Int_t i = 0;i<26;i++)
        eventInfo->fT0amplitude[i] = tzero->GetAmp(i);
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOF[i] = tzero->GetT0TOF()[i];
      for(Int_t i = 0;i<3;i++)
        eventInfo->fT0TOFbest[i] = tzero->GetT0TOFbest()[i];
      eventInfo->fT0pileup = tzero->GetPileupFlag();
      eventInfo->fT0sattelite = tzero->GetSatellite();
    }
  }

  // lines from PWG/FLOW/Tasks/AliFlowTrackCuts.cxx
  /*if(isESD && fFillBayesianPIDInfo){
    //fAliFlowTrackCuts->GetBayesianResponse()->SetDetResponse(esdEvent, eventInfo->fCentrality[1],AliESDpid::kTOF_T0,kFALSE); // centrality = PbPb centrality class (0-100%) or -1 for pp collisions
    fBayesianResponse->SetDetResponse(esdEvent, eventInfo->fCentrality[1],AliESDpid::kTOF_T0,kFALSE); // centrality = PbPb centrality class (0-100%) or -1 for pp collisions
    //fAliFlowTrackCuts->GetESDpid().SetTOFResponse(esdEvent,AliESDpid::kTOF_T0);
  }*/
  //fAliFlowTrackCuts->GetBayesianResponse()->ResetDetOR(1);

  eventInfo->fSPDntracklets = GetSPDTrackletMultiplicity(event, -1.0, 1.0);
  for(Int_t ieta=0; ieta<32; ++ieta)
    eventInfo->fSPDntrackletsEta[ieta] = GetSPDTrackletMultiplicity(event, -1.6+0.1*ieta, -1.6+0.1*(ieta+1));
  
  AliVVZERO* vzero = event->GetVZEROData();
  for(Int_t i=0;i<64;++i) 
    eventInfo->fVZEROMult[i] = vzero->GetMultiplicity(i);  
  Float_t multVZERO = 0.0;
  for(Int_t i=0;i<32;++i) multVZERO +=  vzero->GetMultiplicity(i);
  eventInfo->fVZEROTotalMult[1] = multVZERO;
  multVZERO = 0.0;
  for(Int_t i=32;i<64;++i) multVZERO +=  vzero->GetMultiplicity(i);
  eventInfo->fVZEROTotalMult[0] = multVZERO;
  
  if(fFillEventPlaneInfo) {
    AliReducedEventPlaneInfo* ep=new AliReducedEventPlaneInfo();     
    ep->fQvector[AliReducedEventPlaneInfo::kTPC][1][0] = AliDielectronVarManager::GetValue(AliDielectronVarManager::kQnTPCxH2);
    ep->fQvector[AliReducedEventPlaneInfo::kTPC][1][1] = AliDielectronVarManager::GetValue(AliDielectronVarManager::kQnTPCyH2);
    ep->fEventPlaneStatus[AliReducedEventPlaneInfo::kTPC][1] = AliReducedEventPlaneInfo::kRecentered;
    ep->fQvector[AliReducedEventPlaneInfo::kVZEROA][1][0] = AliDielectronVarManager::GetValue(AliDielectronVarManager::kQnV0AxH2);
    ep->fQvector[AliReducedEventPlaneInfo::kVZEROA][1][1] = AliDielectronVarManager::GetValue(AliDielectronVarManager::kQnV0AyH2);
    ep->fEventPlaneStatus[AliReducedEventPlaneInfo::kVZEROA][1] = AliReducedEventPlaneInfo::kRecentered;
    ep->fQvector[AliReducedEventPlaneInfo::kVZEROC][1][0] = AliDielectronVarManager::GetValue(AliDielectronVarManager::kQnV0CxH2);
    ep->fQvector[AliReducedEventPlaneInfo::kVZEROC][1][1] = AliDielectronVarManager::GetValue(AliDielectronVarManager::kQnV0CyH2);
    ep->fEventPlaneStatus[AliReducedEventPlaneInfo::kVZEROC][1] = AliReducedEventPlaneInfo::kRecentered;
    eventInfo->SetEventPlane(ep);
  }
  
  // EMCAL/PHOS clusters
  if(fFillCaloClusterInfo) FillCaloClusters();
  
  // FMD information
  if(fFillFMDInfo&&isESD) FillFMDInfo();
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillCaloClusters() {
  //
  // Fill info about the calorimeter clusters
  //
  AliVEvent* event = InputEvent();
  Int_t nclusters = event->GetNumberOfCaloClusters();

  AliReducedEventInfo* eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
  if(!eventInfo) return;
  
  eventInfo->fNCaloClusters = 0;
  for(Int_t iclus=0; iclus<nclusters; ++iclus) {
    AliVCluster* cluster = event->GetCaloCluster(iclus);
    
    TClonesArray& clusters = *(eventInfo->fCaloClusters);
    AliReducedCaloClusterInfo *reducedCluster=new(clusters[eventInfo->fNCaloClusters]) AliReducedCaloClusterInfo();
    
    reducedCluster->fType    = (cluster->IsEMCAL() ? AliReducedCaloClusterInfo::kEMCAL : AliReducedCaloClusterInfo::kPHOS);
    reducedCluster->fEnergy  = cluster->E();
    reducedCluster->fTrackDx = cluster->GetTrackDx();
    reducedCluster->fTrackDz = cluster->GetTrackDz();
    reducedCluster->fM20     = cluster->GetM20();
    reducedCluster->fM02     = cluster->GetM02();
    reducedCluster->fDispersion = cluster->GetDispersion();
    cluster->GetPosition(reducedCluster->fPosition);
    reducedCluster->fTOF = cluster->GetTOF();
    reducedCluster->fNCells = cluster->GetNCells();
    eventInfo->fNCaloClusters += 1;
  }  // end loop over clusters
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillFMDInfo()
{
  AliAODEvent* aodEvent = AliForwardUtil::GetAODEvent(this);
  if (!aodEvent) {cout<<"didn't get AOD"<<endl; return;}

  //TObject* obj = aodEvent->FindListObject("Forward");  
  //if (!obj) return;

  TH2D* histos[5];
  histos[0] = static_cast<TH2D*>(aodEvent->FindListObject("FMD1I_cache"));  
  histos[1] = static_cast<TH2D*>(aodEvent->FindListObject("FMD2I_cache"));  
  histos[2] = static_cast<TH2D*>(aodEvent->FindListObject("FMD2O_cache"));  
  histos[3] = static_cast<TH2D*>(aodEvent->FindListObject("FMD3I_cache"));  
  histos[4] = static_cast<TH2D*>(aodEvent->FindListObject("FMD3O_cache"));  

  //AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(obj);
  //const TH2D& d2Ndetadphi = aodForward->GetHistogram();

  Float_t m;

  AliReducedEventInfo* eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
  if(!eventInfo) return;
  
  TClonesArray& fmd = *(eventInfo->GetFMD());

  // Loop over eta 
  Int_t nFMD=-1;
  for (Int_t ih = 0; ih < 5; ih++) {
    if(!histos[ih]) continue;
    for (Int_t iEta = 1; iEta <= histos[ih]->GetNbinsX(); iEta++) {

      //Int_t valid = histos[ih]->GetBinContent(iEta, 0);
      //etabin=axeta->FindBin(histos[ih]->GetXaxis()->GetBinCenter(iEta));
      //if (!valid) continue; // No data expected for this eta 
      // Loop over phi 
      for (Int_t iPhi = 1; iPhi <= histos[ih]->GetNbinsY(); iPhi++) {
      m     =  histos[ih]->GetBinContent(iEta, iPhi);
      if(m<1E-6) continue;
      //phibin=axphi->FindBin(histos[ih]->GetYaxis()->GetBinCenter(iPhi));
      nFMD++;
      AliReducedFMDInfo   *reducedFMD=new(fmd[nFMD]) AliReducedFMDInfo();
      reducedFMD->fMultiplicity     =  m;
      reducedFMD->fId               =  iEta*histos[ih]->GetNbinsY()+iPhi;
      if(ih==2||ih==4) reducedFMD->fId*=-1;

      //cout<<ih<<"  "<<iEta<<"  "<<iPhi<<"  "<<reducedFMD->PhiBin()<<"  "<<histos[ih]->GetXaxis()->GetBinCenter(iEta)<<"  "<<histos[ih]->GetYaxis()->GetBinCenter(iPhi)<<"  "<<reducedFMD->Eta()<<"  "<<reducedFMD->Phi()<<"  "<<m<<endl;
      //cout<<ih<<"  "<<iEta<<"  "<<iPhi<<"  "<<etabin<<"  "<<phibin<<"  "<<etabin*nPhi+phibin<<"  "<<reducedFMD->EtaBin(etabin*nPhi+phibin)<<"  "<<reducedFMD->PhiBin(etabin*nPhi+phibin)<<"  "<<endl;
      
        //xc[ih]+=m*TMath::Cos(2.*reducedFMD->Phi());
        //yc[ih]+=m*TMath::Sin(2.*reducedFMD->Phi());

      //cout<<"MINE  "<<iEta<<"  "<<iPhi<<"  "<<d2Ndetadphi.GetXaxis()->GetBinCenter(iEta)<<"  "<<d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi)<<"  "<<reducedFMD->Multiplicity()<<endl;;

      //cout<<iEta<<"  "<<iPhi<<"  "<<d2Ndetadphi.GetXaxis()->GetBinCenter(iEta)<<"  "<<d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi)<<"  "<<reducedFMD->Eta()<<"  "<<reducedFMD->Phi()<<endl;

      }
    }
  }

  //for (Int_t ih = 0; ih < 5; ih++) {
  //cout<<"MINE "<<ih<<"  "<<xc[ih]<<"  "<<yc[ih]<<endl;
  //}

  //AliAODForwardEP fAODEP = AliAODForwardEP();
  //AliFMDEventPlaneFinder  fEventPlaneFinder = AliFMDEventPlaneFinder() ;

  //fEventPlaneFinder.FindEventPlane(event, fAODEP, 


}

//________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeMaker::Rapidity(Double_t r, Double_t z){
  //
  // calculate eta based on radius from beampipe r and distance from interaction point z
  //
  Double_t x = r/z;
  if(z<0) x = x*-1;

  Double_t eta = -1.*TMath::Log((TMath::Sqrt(x*x+1)-1)/x);

  if(z<0) eta = eta*-1;
  return eta;
}

//________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeMaker::Radius(Double_t eta, Double_t z){
  //
  // calculate radius from beampipe based on distance from interaction point z and eta
  //
  Double_t r = 2*TMath::Power(TMath::E(), eta)*z/(TMath::Power(TMath::E(), 2.0*eta)-1);
  return r;
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillMCTruthInfo() 
{
   //
   // fill MC truth info
   //
   Bool_t hasMC = AliDielectronMC::Instance()->HasMC();
   if(!hasMC) return;
   
   AliDielectronMC* mcHandler = AliDielectronMC::Instance();
   
   Int_t nPrimary = mcHandler->GetNPrimaryFromStack();
   
   //cout << "Event+++++++++++++++++++++++++" << endl;
   
   for(Int_t i=0; i<nPrimary; ++i) {
      AliVParticle* particle = mcHandler->GetMCTrackFromMCEvent(i);
      
      // write J/psi's and electrons from J/psi decays
      // TODO: Create a dynamical way to define which particles from the MC stack will be written
      if(!particle) continue;
      Bool_t acceptParticle = kFALSE;
      if(particle->PdgCode()==443) acceptParticle = kTRUE;
      if(TMath::Abs(particle->PdgCode())==411) acceptParticle = kTRUE;
      if(TMath::Abs(particle->PdgCode())==421) acceptParticle = kTRUE;
      if(TMath::Abs(particle->PdgCode())==431) acceptParticle = kTRUE;
      if(TMath::Abs(particle->PdgCode())==4122) acceptParticle = kTRUE;
      AliVParticle* mother = mcHandler->GetMCTrackFromMCEvent(particle->GetMother());
      if(mother && mother->PdgCode()==443) acceptParticle = kTRUE;
      if(mother && TMath::Abs(mother->PdgCode())==411) acceptParticle = kTRUE;
      if(mother && TMath::Abs(mother->PdgCode())==421) acceptParticle = kTRUE;
      if(mother && TMath::Abs(mother->PdgCode())==431) acceptParticle = kTRUE;
      if(mother && TMath::Abs(mother->PdgCode())==4122) acceptParticle = kTRUE;      
      if(!acceptParticle) continue;
      
      TClonesArray& tracks = *(fReducedEvent->fTracks);
      AliReducedBaseTrack* reducedParticle=NULL;
      if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithBaseTracks)
         reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliReducedBaseTrack();
      if(fTreeWritingOption==kBaseEventsWithFullTracks || fTreeWritingOption==kFullEventsWithFullTracks)
         reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliReducedTrackInfo();
      
      reducedParticle->PxPyPz(particle->Px(), particle->Py(), particle->Pz());
      reducedParticle->fQualityFlags |= (ULong_t(1)<<63);               // this means that this is a pure MC track
      
      Int_t nDaughters = (particle->PdgCode()==443 ? particle->GetLastDaughter() - particle->GetFirstDaughter() + 1 : 0);
      if(nDaughters==2) reducedParticle->fQualityFlags |= (ULong_t(1)<<62);    // J/psi -> e+e-
      if(nDaughters>2) reducedParticle->fQualityFlags |= (ULong_t(1)<<61);       // J/psi -> e+e- + X 
      
      AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(reducedParticle);
      if(!trackInfo) continue;
      
      trackInfo->fMCLabels[0] = particle->GetLabel();
      trackInfo->fMCPdg[0] = particle->PdgCode();
      trackInfo->fMCMom[0] = particle->Px();
      trackInfo->fMCMom[1] = particle->Py();
      trackInfo->fMCMom[2] = particle->Pz();
      
      if(mother) {
        trackInfo->fMCLabels[1] = mother->GetLabel();
        trackInfo->fMCPdg[1] = mother->PdgCode();
        if(particle->PdgCode()==443)
          reducedParticle->fQualityFlags |= (ULong_t(1)<<60);    // secondary J/psi
        
        AliVParticle* grandmother = mcHandler->GetMCTrackFromMCEvent(mother->GetMother());
        if(grandmother) {
           trackInfo->fMCLabels[2] = grandmother->GetLabel();
           trackInfo->fMCPdg[2] = grandmother->PdgCode();
           
           AliVParticle* grandgrandmother = mcHandler->GetMCTrackFromMCEvent(grandmother->GetMother());
           if(grandgrandmother) {
              trackInfo->fMCLabels[3] = grandgrandmother->GetLabel();
              trackInfo->fMCPdg[3] = grandgrandmother->PdgCode();
           }
        }
      }
      
      /*cout << "particle label/pdg/mlabel/mpdg/px/py/pz/ndaughters/first/last :: " << trackInfo->fMCLabels[0] << "/" << trackInfo->fMCPdg[0] << "/"
        << trackInfo->fMCLabels[1] << "/" << trackInfo->fMCPdg[1] << "/" << reducedParticle->Px() << "/"
        << reducedParticle->Py() << "/" << reducedParticle->Pz() << "/" << nDaughters << "/" << particle->GetFirstDaughter() << "/"
        << particle->GetLastDaughter() << endl; */
        
      fReducedEvent->fNtracks[1] += 1;  
   }
}


//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillTrackInfo() 
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
  
  Bool_t hasMC = AliDielectronMC::Instance()->HasMC();
  
  // find all the tracks which belong to a V0 stored in the reduced event
  UShort_t trackIdsV0[4][20000]={{0}};
  UShort_t trackIdsPureV0[4][20000]={{0}};
  Int_t nV0LegsTagged[4] = {0}; Int_t nPureV0LegsTagged[4] = {0};
  Bool_t leg1Found[4]; Bool_t leg2Found[4];
  for(Int_t iv0=0;iv0<fReducedEvent->fNV0candidates[1];++iv0) {
    AliReducedPairInfo* pair = fReducedEvent->GetV0Pair(iv0);
    if(!pair) continue;
    Int_t pairId = 0; Bool_t isPureV0 = kFALSE;
    if(pair->fCandidateId==AliReducedPairInfo::kGammaConv) {
      pairId=0;
      if(pair->IsPureV0Gamma()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPairInfo::kK0sToPiPi) {
      pairId=1;
      if(pair->IsPureV0K0s()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPairInfo::kLambda0ToPPi) {
      pairId=2;
      if(pair->IsPureV0Lambda()) isPureV0 = kTRUE;
    }
    if(pair->fCandidateId==AliReducedPairInfo::kALambda0ToPPi) {
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
      
  Int_t pidtypes[4] = {AliPID::kElectron,AliPID::kPion,AliPID::kKaon,AliPID::kProton};
  AliESDtrack* esdTrack=0;
  AliAODTrack* aodTrack=0;
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t trackId = 0; 
  Bool_t usedForV0[4] = {kFALSE}; 
  Bool_t usedForPureV0[4] = {kFALSE};
  Bool_t usedForV0Or = kFALSE;
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
          break;
        }
      }
    }
        
    ULong_t status = (isESD ? esdTrack->GetStatus() : aodTrack->GetStatus());
    
    AliReducedEventInfo* eventInfo = NULL; 
        
    if(fTreeWritingOption==kFullEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithFullTracks) {
      eventInfo = dynamic_cast<AliReducedEventInfo*>(fReducedEvent);
      for(Int_t ibit=0; ibit<32; ++ibit) {
         if(status & (ULong_t(1)<<ibit)) {
            eventInfo->fNtracksPerTrackingFlag[ibit] += 1;
         }
      }
    }
    
    //apply track cuts
    if(!usedForV0Or && fTrackFilter && !fTrackFilter->IsSelected(particle)) continue;
    
    TClonesArray& tracks = *(fReducedEvent->fTracks);
    AliReducedBaseTrack* reducedParticle=NULL;
    if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithBaseTracks)
      reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliReducedBaseTrack();
    if(fTreeWritingOption==kBaseEventsWithFullTracks || fTreeWritingOption==kFullEventsWithFullTracks)
       reducedParticle=new(tracks[fReducedEvent->fNtracks[1]]) AliReducedTrackInfo();
        
    Double_t values[AliDielectronVarManager::kNMaxValues];
    // set the fill map (all 1's) for the AliDielectronVarManager
    AliDielectronVarManager::SetFillMap(fUsedVars);
    AliDielectronVarManager::Fill(particle, values);
    
    reducedParticle->PtPhiEta(values[AliDielectronVarManager::kPt],values[AliDielectronVarManager::kPhi],values[AliDielectronVarManager::kEta]);
    reducedParticle->fCharge        = values[AliDielectronVarManager::kCharge];
    
    if(fFlowTrackFilter) {
       // switch on the first bit if this particle should be used for the event plane
       if(fFlowTrackFilter->IsSelected(particle)) reducedParticle->fQualityFlags |= (ULong_t(1)<<0);
    }
    for(Int_t iV0type=0;iV0type<4;++iV0type) {
       if(usedForV0[iV0type]) reducedParticle->fQualityFlags |= (ULong_t(1)<<(iV0type+1));
       if(usedForPureV0[iV0type]) reducedParticle->fQualityFlags |= (ULong_t(1)<<(iV0type+8));
    }
        
    if(isESD) {
       for(Int_t idx=0; idx<3; ++idx) if(esdTrack->GetKinkIndex(idx)>0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(5+idx));
       for(Int_t idx=0; idx<3; ++idx) if(esdTrack->GetKinkIndex(idx)<0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(12+idx));
       
       //check is track passes bayesian combined TOF+TPC pid cut
       //Bool_t goodtrack = (esdTrack->GetStatus() & AliESDtrack::kTOFout) &&
       //                   (esdTrack->GetStatus() & AliESDtrack::kTIME) &&
       //                   (esdTrack->GetTOFsignal() > 12000) &&
       //                   (esdTrack->GetTOFsignal() < 100000) &&
       //                   (esdTrack->GetIntegratedLength() > 365);
       //Float_t mismProb = fBayesianResponse->GetTOFMismProb(); // mismatch Bayesian probabilities
       /*if(fFillBayesianPIDInfo) {
         fBayesianResponse->ComputeProb(esdTrack,fReducedEvent->fCentrality[1]); // centrality is needed for mismatch fraction
         Int_t kTPC = fBayesianResponse->GetCurrentMask(0); // is TPC on
         if( kTPC){
            //fAliFlowTrackCuts->GetBayesianResponse()->ComputeProb(esdTrack,fReducedEvent->fCentrality[1]); // centrality is needed for mismatch fraction
            //Int_t kTPC = fAliFlowTrackCuts->GetBayesianResponse()->GetCurrentMask(0); // is TPC on
            //Float_t *probabilities = fAliFlowTrackCuts->GetBayesianResponse()->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mismatch parameterization
            Float_t *probabilities = fBayesianResponse->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mismatch parameterization
            if(probabilities[AliPID::kElectron]>0.5) reducedParticle->fQualityFlags |= (ULong_t(1)<<15);
            if(probabilities[AliPID::kPion    ]>0.5) reducedParticle->fQualityFlags |= (ULong_t(1)<<16);
            if(probabilities[AliPID::kKaon    ]>0.5) reducedParticle->fQualityFlags |= (ULong_t(1)<<17);
            if(probabilities[AliPID::kProton  ]>0.5) reducedParticle->fQualityFlags |= (ULong_t(1)<<18);
          
            for(Int_t ipid=0; ipid<4; ipid++){
               if(probabilities[pidtypes[ipid]]>0.7) reducedParticle->fQualityFlags |= (ULong_t(1)<<19);
               if(probabilities[pidtypes[ipid]]>0.8) reducedParticle->fQualityFlags |= (ULong_t(1)<<20);
               if(probabilities[pidtypes[ipid]]>0.9) reducedParticle->fQualityFlags |= (ULong_t(1)<<21);
            }
          
            //reducedParticle->fBayes[0]   = probabilities[0];
            //reducedParticle->fBayes[1]   = probabilities[2];
            //reducedParticle->fBayes[2]   = probabilities[3];
            //reducedParticle->fBayes[3]   = probabilities[4];
       }
     } */
   }
   if(isAOD) {
      for(Int_t idx=0; idx<3; ++idx) if(aodTrack->GetKinkIndex(idx)>0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(5+idx));
      for(Int_t idx=0; idx<3; ++idx) if(aodTrack->GetKinkIndex(idx)<0) reducedParticle->fQualityFlags |= (ULong_t(1)<<(12+idx));
      for(Int_t idx=0; idx<11; ++idx) if(aodTrack->TestFilterBit(BIT(idx))) reducedParticle->SetQualityFlag(15+idx);
   }
   
   // If we want to write only AliReducedBaseTrack objects, then we stop here
   if(fTreeWritingOption==kBaseEventsWithBaseTracks || fTreeWritingOption==kFullEventsWithBaseTracks) {
     fReducedEvent->fNtracks[1] += 1;
     continue;
   }    
    
    AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(reducedParticle);
    if(!trackInfo) continue;
    
    trackInfo->fStatus        = status;
    trackInfo->fMomentumInner = values[AliDielectronVarManager::kPIn];
    trackInfo->fDCA[0]        = values[AliDielectronVarManager::kImpactParXY];
    trackInfo->fDCA[1]        = values[AliDielectronVarManager::kImpactParZ];
    trackInfo->fTrackLength   = values[AliDielectronVarManager::kTrackLength];
    
    trackInfo->fITSclusterMap = (UChar_t)values[AliDielectronVarManager::kITSclusterMap];
    trackInfo->fITSsignal     = values[AliDielectronVarManager::kITSsignal];
    trackInfo->fITSnSig[0]    = values[AliDielectronVarManager::kITSnSigmaEle];
    trackInfo->fITSnSig[1]    = values[AliDielectronVarManager::kITSnSigmaPio];
    trackInfo->fITSnSig[2]    = values[AliDielectronVarManager::kITSnSigmaKao];
    trackInfo->fITSnSig[3]    = values[AliDielectronVarManager::kITSnSigmaPro];
    trackInfo->fITSchi2       = values[AliDielectronVarManager::kITSchi2Cl];
    
    trackInfo->fTPCNcls      = (UChar_t)values[AliDielectronVarManager::kNclsTPC];
    trackInfo->fTPCNclsF     = (UChar_t)values[AliDielectronVarManager::kNFclsTPC];
    trackInfo->fTPCNclsShared = (UChar_t)values[AliDielectronVarManager::kNclsSTPC];
    trackInfo->fTPCCrossedRows = values[AliDielectronVarManager::kNFclsTPCr];
    trackInfo->fTPCsignal    = values[AliDielectronVarManager::kTPCsignal];
    trackInfo->fTPCsignalN   = values[AliDielectronVarManager::kTPCsignalN];
    trackInfo->fTPCnSig[0]   = values[AliDielectronVarManager::kTPCnSigmaEle];
    trackInfo->fTPCnSig[1]   = values[AliDielectronVarManager::kTPCnSigmaPio];
    trackInfo->fTPCnSig[2]   = values[AliDielectronVarManager::kTPCnSigmaKao];
    trackInfo->fTPCnSig[3]   = values[AliDielectronVarManager::kTPCnSigmaPro];
    trackInfo->fTPCClusterMap = EncodeTPCClusterMap(particle, isAOD);
    trackInfo->fTPCchi2       = values[AliDielectronVarManager::kTPCchi2Cl];
    trackInfo->fTPCActiveLength = values[AliDielectronVarManager::kTPCActiveLength];
    trackInfo->fTPCGeomLength = values[AliDielectronVarManager::kTPCGeomLength];
        
    trackInfo->fTOFbeta      = values[AliDielectronVarManager::kTOFbeta];
    trackInfo->fTOFtime      = values[AliDielectronVarManager::kTOFsignal]-pidResponse->GetTOFResponse().GetTimeZero();
    trackInfo->fTOFmismatchProbab = values[AliDielectronVarManager::kTOFmismProb];
    trackInfo->fTOFnSig[0]   = values[AliDielectronVarManager::kTOFnSigmaEle];
    trackInfo->fTOFnSig[1]   = values[AliDielectronVarManager::kTOFnSigmaPio];
    trackInfo->fTOFnSig[2]   = values[AliDielectronVarManager::kTOFnSigmaKao];
    trackInfo->fTOFnSig[3]   = values[AliDielectronVarManager::kTOFnSigmaPro];
    
    Double_t trdProbab[AliPID::kSPECIES]={0.0};
    if(isESD) {
       trackInfo->fMassForTracking = esdTrack->GetMassForTracking();
       
       AliESDEvent* esdEvent = static_cast<AliESDEvent*>(InputEvent());
       AliESDVertex* eventVtx = const_cast<AliESDVertex*>(esdEvent->GetPrimaryVertexTracks());
       if(fReducedEvent->fRunNo>245000. && fReducedEvent->fRunNo<247000.)
         trackInfo->fChi2TPCConstrainedVsGlobal = esdTrack->GetChi2TPCConstrainedVsGlobal(eventVtx);
       
      trackInfo->fTrackId          = (UShort_t)esdTrack->GetID();
      const AliExternalTrackParam* tpcInner = esdTrack->GetTPCInnerParam();

      //trackInfo->fITSSharedClusterMap = esdTrack->GetITSSharedClusterMap();
      for(Int_t i=0; i<6; ++i) {
         if(esdTrack->HasSharedPointOnITSLayer(i)) trackInfo->fITSSharedClusterMap |= (1<<i);
      }
      
      Float_t xyDCA,zDCA;
      Double_t helixinfo[6];
      if(tpcInner){
        trackInfo->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
        trackInfo->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
        trackInfo->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
        esdTrack->GetImpactParametersTPC(xyDCA,zDCA);
        trackInfo->fTPCDCA[0]     = xyDCA;
        trackInfo->fTPCDCA[1]     = zDCA;
        
        // helix information (Alex Chauvin)
        tpcInner->GetHelixParameters(helixinfo,InputEvent()->GetMagneticField());
        if(helixinfo[2] < 0) helixinfo[2] = helixinfo[2] + 2*TMath::Pi();
        helixinfo[2] -= TMath::Pi()/2.;
        trackInfo->fHelixCenter[0]= helixinfo[5]+(TMath::Cos(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixCenter[1]= helixinfo[0]+(TMath::Sin(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixRadius   = TMath::Abs(1./helixinfo[4]);
      }
      
      trackInfo->fTOFdeltaBC    = esdTrack->GetTOFDeltaBC();
      trackInfo->fTOFdx         = esdTrack->GetTOFsignalDx();
      trackInfo->fTOFdz         = esdTrack->GetTOFsignalDz();
      trackInfo->fTOFchi2       = esdTrack->GetTOFchi2();
      
      trackInfo->fTRDntracklets[0] = esdTrack->GetTRDntracklets();
      trackInfo->fTRDntracklets[1] = esdTrack->GetTRDntrackletsPID();
      pidResponse->ComputeTRDProbability(esdTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ1D);
      trackInfo->fTRDpid[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpid[1]    = trdProbab[AliPID::kPion];
      pidResponse->ComputeTRDProbability(esdTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ2D);
      trackInfo->fTRDpidLQ2D[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpidLQ2D[1]    = trdProbab[AliPID::kPion];
                
      if(esdTrack->IsEMCAL()) trackInfo->fCaloClusterId = esdTrack->GetEMCALcluster();
      if(esdTrack->IsPHOS()) trackInfo->fCaloClusterId = esdTrack->GetPHOScluster();
      
      if(fFillMCInfo && hasMC) {
         AliMCParticle* truthParticle = AliDielectronMC::Instance()->GetMCTrack(esdTrack);
         if(truthParticle) {
           trackInfo->fMCMom[0] = truthParticle->Px();
           trackInfo->fMCMom[1] = truthParticle->Py();
           trackInfo->fMCMom[2] = truthParticle->Pz();
           trackInfo->fMCFreezeout[0] = truthParticle->Xv();
           trackInfo->fMCFreezeout[1] = truthParticle->Yv();
           trackInfo->fMCFreezeout[2] = truthParticle->Zv();
           trackInfo->fMCLabels[0] = esdTrack->GetLabel();
           trackInfo->fMCPdg[0] = truthParticle->PdgCode();
           trackInfo->fMCGeneratorIndex = truthParticle->GetGeneratorIndex();
           if(truthParticle->PdgCode()!=-9999 && esdTrack->GetLabel()!=-9999) trackInfo->fQualityFlags |= (ULong_t(1)<<22);   // means the track has MC truth info
           
           AliMCParticle* motherTruth = AliDielectronMC::Instance()->GetMCTrackMother(truthParticle);
           if(motherTruth) {
             trackInfo->fMCLabels[1] = truthParticle->GetMother();
             trackInfo->fMCPdg[1] = motherTruth->PdgCode();
          }
          
          AliMCParticle* grandmotherTruth = NULL;
          if(motherTruth) grandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(motherTruth);
          if(grandmotherTruth) {
             trackInfo->fMCLabels[2] = motherTruth->GetMother();
             trackInfo->fMCPdg[2] = grandmotherTruth->PdgCode();
          }
           
           AliMCParticle* grandgrandmotherTruth = NULL;
           if(grandmotherTruth) grandgrandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(grandmotherTruth);
           if(grandgrandmotherTruth) {
              trackInfo->fMCLabels[3] = grandmotherTruth->GetMother();
              trackInfo->fMCPdg[3] = grandgrandmotherTruth->PdgCode();
           }
         }
      }
    }  // end if(isESD)
    if(isAOD) {
      trackInfo->fMassForTracking = aodTrack->GetMassForTracking();
      trackInfo->fChi2TPCConstrainedVsGlobal = aodTrack->GetChi2TPCConstrainedVsGlobal(); 
      
      //trackInfo->fITSSharedClusterMap = aodTrack->GetITSSharedClusterMap();
      for(Int_t i=0; i<6; ++i) {
         if(aodTrack->HasSharedPointOnITSLayer(i)) trackInfo->fITSSharedClusterMap |= (1<<i);
      }
      
      const AliExternalTrackParam* tpcInner = aodTrack->GetInnerParam();
      Float_t xyDCA,zDCA;
      Double_t helixinfo[6];
      if(tpcInner){
        trackInfo->fTPCPhi        = (tpcInner ? tpcInner->Phi() : 0.0);
        trackInfo->fTPCPt         = (tpcInner ? tpcInner->Pt() : 0.0);
        trackInfo->fTPCEta        = (tpcInner ? tpcInner->Eta() : 0.0);
      
        aodTrack->GetImpactParametersTPC(xyDCA,zDCA);
        trackInfo->fTPCDCA[0]     = xyDCA;
        trackInfo->fTPCDCA[1]     = zDCA;
      
        // helix information (Alex Chauvin)
        tpcInner->GetHelixParameters(helixinfo,InputEvent()->GetMagneticField());
        if(helixinfo[2] < 0) helixinfo[2] = helixinfo[2] + 2*TMath::Pi();
        helixinfo[2] -= TMath::Pi()/2.;
        trackInfo->fHelixCenter[0]= helixinfo[5]+(TMath::Cos(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixCenter[1]= helixinfo[0]+(TMath::Sin(helixinfo[2])*TMath::Abs(1./helixinfo[4])*copysignf(1.0, InputEvent()->GetMagneticField()*values[AliDielectronVarManager::kCharge]));
        trackInfo->fHelixRadius   = TMath::Abs(1./helixinfo[4]);
      }
      
      trackInfo->fTOFdz         = aodTrack->GetTOFsignalDz();
      trackInfo->fTOFdeltaBC = eventInfo->fBC - aodTrack->GetTOFBunchCrossing();
      
      trackInfo->fTrackId = aodTrack->GetID(); 
      trackInfo->fTRDntracklets[0] = aodTrack->GetTRDntrackletsPID();
      trackInfo->fTRDntracklets[1] = aodTrack->GetTRDntrackletsPID();
      pidResponse->ComputeTRDProbability(aodTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ1D);
      trackInfo->fTRDpid[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpid[1]    = trdProbab[AliPID::kPion];
      pidResponse->ComputeTRDProbability(aodTrack,AliPID::kSPECIES,trdProbab,AliTRDPIDResponse::kLQ2D);
      trackInfo->fTRDpidLQ2D[0]    = trdProbab[AliPID::kElectron];
      trackInfo->fTRDpidLQ2D[1]    = trdProbab[AliPID::kPion];
      
      if(aodTrack->IsEMCAL()) trackInfo->fCaloClusterId = aodTrack->GetEMCALcluster();
      if(aodTrack->IsPHOS()) trackInfo->fCaloClusterId = aodTrack->GetPHOScluster();
      
      if(fFillMCInfo && hasMC) {
         AliAODMCParticle* truthParticle = AliDielectronMC::Instance()->GetMCTrack(aodTrack);
         if(truthParticle) {
            trackInfo->fMCMom[0] = truthParticle->Px();
            trackInfo->fMCMom[1] = truthParticle->Py();
            trackInfo->fMCMom[2] = truthParticle->Pz();
            trackInfo->fMCFreezeout[0] = truthParticle->Xv();
            trackInfo->fMCFreezeout[1] = truthParticle->Yv();
            trackInfo->fMCFreezeout[2] = truthParticle->Zv();
            trackInfo->fMCLabels[0] = aodTrack->GetLabel();
            trackInfo->fMCPdg[0] = truthParticle->PdgCode();
            trackInfo->fMCGeneratorIndex = truthParticle->GetGeneratorIndex();
            if(truthParticle->PdgCode()!=-9999 && aodTrack->GetLabel()!=-9999) trackInfo->fQualityFlags |= (ULong_t(1)<<22);
            
            AliAODMCParticle* motherTruth = AliDielectronMC::Instance()->GetMCTrackMother(truthParticle);
            if(motherTruth) {
               trackInfo->fMCLabels[1] = truthParticle->GetMother();
               trackInfo->fMCPdg[1] = motherTruth->PdgCode();
            }
            
            AliAODMCParticle* grandmotherTruth = NULL;
            if(motherTruth) grandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(motherTruth);
            if(grandmotherTruth) {
               trackInfo->fMCLabels[2] = motherTruth->GetMother();
               trackInfo->fMCPdg[2] = grandmotherTruth->PdgCode();
            }
            
            AliAODMCParticle* grandgrandmotherTruth = NULL;
            if(grandmotherTruth) grandgrandmotherTruth = AliDielectronMC::Instance()->GetMCTrackMother(grandmotherTruth);
            if(grandgrandmotherTruth) {
               trackInfo->fMCLabels[3] = grandmotherTruth->GetMother();
               trackInfo->fMCPdg[3] = grandgrandmotherTruth->PdgCode();
            }
         }
      }
    }  // end if(isAOD)

    fReducedEvent->fNtracks[1] += 1;
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskReducedTreeMaker::FillV0PairInfo() 
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
    AliReducedPairInfo* k0sReducedPair     = FillV0PairInfo(v0, AliReducedPairInfo::kK0sToPiPi,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPairInfo* lambdaReducedPair  = FillV0PairInfo(v0, AliReducedPairInfo::kLambda0ToPPi,  legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPairInfo* alambdaReducedPair = FillV0PairInfo(v0, AliReducedPairInfo::kALambda0ToPPi, legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    AliReducedPairInfo* gammaReducedPair   = FillV0PairInfo(v0, AliReducedPairInfo::kGammaConv,     legPos, legNeg, &primaryVertexKF, v0ChargesAreCorrect);
    
    if(fFillK0s && goodK0s && k0sReducedPair->fMass[0]>fK0sMassRange[0] && k0sReducedPair->fMass[0]<fK0sMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodK0sPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*k0sReducedPair);
      goodK0sPair->fMass[0] = k0sReducedPair->fMass[0];
      goodK0sPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodK0sPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodK0sPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodK0s) goodK0sPair->fQualityFlags |= (ULong_t(1)<<1);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodK0s=kFALSE;}
    if(fFillLambda && goodLambda && lambdaReducedPair->fMass[0]>fLambdaMassRange[0] && lambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodLambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*lambdaReducedPair);
      goodLambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodLambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodLambdaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodLambda) goodLambdaPair->fQualityFlags |= (ULong_t(1)<<2);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodLambda=kFALSE;}
    if(fFillALambda && goodALambda && alambdaReducedPair->fMass[0]>fLambdaMassRange[0] && alambdaReducedPair->fMass[0]<fLambdaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodALambdaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*alambdaReducedPair);
      goodALambdaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodALambdaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodALambdaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodALambda) goodALambdaPair->fQualityFlags |= (ULong_t(1)<<3);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodALambda = kFALSE;}
    //cout << "gamma mass: " << gammaReducedPair->fMass[0] << endl;
    if(fFillGammaConversions && goodGamma && gammaReducedPair->fMass[0]>fGammaMassRange[0] && gammaReducedPair->fMass[0]<fGammaMassRange[1]) {
      TClonesArray& tracks = *(fReducedEvent->fCandidates);
      AliReducedPairInfo *goodGammaPair = new (tracks[fReducedEvent->fNV0candidates[1]]) AliReducedPairInfo(*gammaReducedPair);
      goodGammaPair->fMass[0] = k0sReducedPair->fMass[0];
      goodGammaPair->fMass[1] = lambdaReducedPair->fMass[0];
      goodGammaPair->fMass[2] = alambdaReducedPair->fMass[0];
      goodGammaPair->fMass[3] = gammaReducedPair->fMass[0];
      if(veryGoodGamma) goodGammaPair->fQualityFlags |= (ULong_t(1)<<4);
      fReducedEvent->fNV0candidates[1] += 1;
    } else {goodGamma=kFALSE;}
    delete k0sReducedPair;
    delete lambdaReducedPair;
    delete alambdaReducedPair;
    delete gammaReducedPair;
  }   // end loop over V0s
}


//_________________________________________________________________________________
AliReducedPairInfo* AliAnalysisTaskReducedTreeMaker::FillV0PairInfo(AliESDv0* v0, Int_t id, 
						    AliESDtrack* legPos, AliESDtrack* legNeg,
						    AliKFVertex* vtxKF, Bool_t chargesAreCorrect) {
  //
  // Create a reduced V0 object and fill it
  //
  AliReducedPairInfo* reducedPair=new AliReducedPairInfo();  
  reducedPair->fCandidateId = id;
  reducedPair->fPairType    = v0->GetOnFlyStatus();    // on the fly status
  reducedPair->fLegIds[0]   = legPos->GetID();
  reducedPair->fLegIds[1]   = legNeg->GetID();
  if(!reducedPair->fPairType) {    // offline
    UInt_t pidPos = AliPID::kPion;
    if(id==AliReducedPairInfo::kLambda0ToPPi) pidPos = AliPID::kProton;
    if(id==AliReducedPairInfo::kGammaConv) pidPos = AliPID::kElectron;
    UInt_t pidNeg = AliPID::kPion;
    if(id==AliReducedPairInfo::kALambda0ToPPi) pidNeg = AliPID::kProton;
    if(id==AliReducedPairInfo::kGammaConv) pidNeg = AliPID::kElectron;
    reducedPair->fMass[0]      = v0->GetEffMass(pidPos, pidNeg);
    reducedPair->fIsCartesian  = kFALSE;
    reducedPair->fP[1]         = v0->Phi();
    if(reducedPair->fP[1]<0.0) reducedPair->fP[1] = 2.0*TMath::Pi() + reducedPair->fP[1];  // converted to [0,2pi]
    reducedPair->fP[0]         = v0->Pt();
    reducedPair->fP[2]         = v0->Eta();
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
    if(id==AliReducedPairInfo::kLambda0ToPPi) pdgPos = 2212;
    if(id==AliReducedPairInfo::kGammaConv) pdgPos = -11;
    Int_t pdgNeg = -211;
    if(id==AliReducedPairInfo::kALambda0ToPPi) pdgNeg = -2212;
    if(id==AliReducedPairInfo::kGammaConv) pdgNeg = 11;
    AliKFParticle negKF(*(negHelix), pdgPos);
    AliKFParticle posKF(*(posHelix), pdgNeg);
    AliKFParticle v0Refit;
    v0Refit += negKF;
    v0Refit += posKF;
    Double_t massFit=0.0, massErrFit=0.0;
    v0Refit.GetMass(massFit,massErrFit);
    reducedPair->fMass[0] = massFit;
    reducedPair->fIsCartesian  = kFALSE;
    reducedPair->fP[1]         = v0Refit.GetPhi();
    if(reducedPair->fP[1]<0.0) reducedPair->fP[1] = 2.0*TMath::Pi() + reducedPair->fP[1];  // converted to [0,2pi]
    reducedPair->fP[0]         = v0Refit.GetPt();
    reducedPair->fP[2]         = v0Refit.GetEta();
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
UChar_t AliAnalysisTaskReducedTreeMaker::EncodeTPCClusterMap(AliVParticle* track, Bool_t isAOD) {
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
Int_t AliAnalysisTaskReducedTreeMaker::GetSPDTrackletMultiplicity(AliVEvent* event, Float_t lowEta, Float_t highEta) {
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
void AliAnalysisTaskReducedTreeMaker::FinishTaskOutput()
{
  //
  // Finish Task 
  //
  PostData(2, fTree);
}
