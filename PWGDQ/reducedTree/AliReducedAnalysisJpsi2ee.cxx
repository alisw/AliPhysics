//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisJpsi2ee.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisJpsi2ee);


//___________________________________________________________________________
AliReducedAnalysisJpsi2ee::AliReducedAnalysisJpsi2ee() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler()),
  fEventCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts()
{
  //
  // default constructor
  //
   
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2ee::AliReducedAnalysisJpsi2ee(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler()),
  fEventCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts()
{
  //
  // named constructor
  //
   fEventCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
   fPreFilterTrackCuts.SetOwner(kTRUE);
   fPairCuts.SetOwner(kTRUE);
   fPreFilterPairCuts.SetOwner(kTRUE);
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2ee::~AliReducedAnalysisJpsi2ee() 
{
  //
  // destructor
  //
   fEventCuts.Clear("C"); fTrackCuts.Clear("C"); fPreFilterTrackCuts.Clear("C"); fPreFilterPairCuts.Clear("C"); fPairCuts.Clear("C");
   if(fHistosManager) delete fHistosManager;
   if(fMixingHandler) delete fMixingHandler;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2ee::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fEventCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if(values) { if(!cut->IsSelected(event, values)) return kFALSE; }
    else { if(!cut->IsSelected(event)) return kFALSE; }
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2ee::IsTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fTrackCuts.GetEntries()==0) return kTRUE;
  track->ResetFlags();
  
  for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track, values)) track->SetFlag(i); }
    else { if(!cut->IsSelected(track)) return kFALSE; }
    
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2ee::IsTrackPrefilterSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
   //
   // apply event cuts
   //
   if(fPreFilterTrackCuts.GetEntries()==0) return kTRUE;
   
   for(Int_t i=0; i<fPreFilterTrackCuts.GetEntries(); ++i) {
      // if there are more cuts specified, we apply an AND on all of them
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPreFilterTrackCuts.At(i);
      if(values) { if(!cut->IsSelected(track, values)) return kFALSE; }
      else { if(!cut->IsSelected(track)) return kFALSE; }
   }
   return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2ee::IsPairSelected(AliReducedBaseTrack* pair, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fPairCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
    if(values) { if(!cut->IsSelected(pair, values)) return kFALSE; }
    else { if(!cut->IsSelected(pair)) return kFALSE; }
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2ee::IsPairPreFilterSelected(AliReducedBaseTrack* pair, Float_t* values/*=0x0*/) {
   //
   // apply event cuts
   //
   if(fPreFilterPairCuts.GetEntries()==0) return kTRUE;
   // loop over all the cuts and make a logical and between all cuts in the list
   for(Int_t i=0; i<fPreFilterPairCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPreFilterPairCuts.At(i);
      if(values) { if(!cut->IsSelected(pair, values)) return kFALSE; }
      else { if(!cut->IsSelected(pair)) return kFALSE; }
   }
   return kTRUE;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::Process() {
  //
  // process the current event
  //  
  if(!fEvent) return;
  
  AliReducedVarManager::SetEvent(fEvent);
  
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_NoCuts", fValues);
  if(fEvent->IsA()==AliReducedEventInfo::Class()) {
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("OnlineTriggers_NoCuts", fValues);
    }
  }
  if(!IsEventSelected(fEvent)) return;
  
  AliReducedEventInfo* eventInfo = NULL;
  if(fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
  
  if(eventInfo) {
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("OnlineTriggers_AfterCuts", fValues);
      for(UShort_t i=0; i<32; ++i) {
        AliReducedVarManager::FillL0TriggerInputs(eventInfo, i, fValues);
        fHistosManager->FillHistClass("OnlineTriggers_vs_L0TrigInputs", fValues);
      }
      for(UShort_t i=0; i<32; ++i) {
        AliReducedVarManager::FillL1TriggerInputs(eventInfo, i, fValues);
        fHistosManager->FillHistClass("OnlineTriggers_vs_L1TrigInputs", fValues);
      }
      for(UShort_t i=0; i<16; ++i) {
        AliReducedVarManager::FillL2TriggerInputs(eventInfo, i, fValues);
        fHistosManager->FillHistClass("OnlineTriggers_vs_L2TrigInputs", fValues);
      }
    }
  }
    
  for(UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
    fHistosManager->FillHistClass("EvtTags", fValues);
  }
  
  if(eventInfo) {
    for(UShort_t ibit=0; ibit<32; ++ibit) {
      AliReducedVarManager::FillL0TriggerInputs(eventInfo, ibit, fValues);
      fHistosManager->FillHistClass("L0TriggerInput", fValues);
    }
    for(UShort_t ibit=0; ibit<32; ++ibit) {
      AliReducedVarManager::FillL1TriggerInputs(eventInfo, ibit, fValues);
      fHistosManager->FillHistClass("L1TriggerInput", fValues);
    }
    for(UShort_t ibit=0; ibit<16; ++ibit) {
      AliReducedVarManager::FillL2TriggerInputs(eventInfo, ibit, fValues);
      fHistosManager->FillHistClass("L2TriggerInput", fValues);
    }
    
    for(Int_t icl=0; icl<eventInfo->GetNCaloClusters(); ++icl) {
      AliReducedVarManager::FillCaloClusterInfo(eventInfo->GetCaloCluster(icl), fValues);
      fHistosManager->FillHistClass("CaloClusters", fValues);
    }
  }
     
  AliReducedBaseTrack* track = 0x0;
  TClonesArray* trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t it=0; it<fEvent->NTracks(); ++it) {
    track = (AliReducedBaseTrack*)nextTrack();
    
    if(!IsTrackSelected(track)) continue;
    AliReducedVarManager::FillTrackInfo(track,fValues);
    fHistosManager->FillHistClass("TrackQA_AllTracks_ITS_TPC_TRD_TOF_EMCAL", fValues);
    
    AliReducedTrackInfo* trackInfo = NULL;
    if(track->IsA()==AliReducedTrackInfo::Class()) trackInfo = (AliReducedTrackInfo*)track;
    
    if(trackInfo) {
      for(UShort_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
        AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
        fHistosManager->FillHistClass("TrackingFlags", fValues);
      }
    }
    for(UShort_t iflag=0; iflag<64; ++iflag) {
      AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues);
      fHistosManager->FillHistClass("TrackQualityFlags", fValues);
    }
    if(trackInfo) {
      for(Int_t iLayer=0; iLayer<6; ++iLayer) {
        AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
        fHistosManager->FillHistClass("ITSclusterMap", fValues);
      }
      for(Int_t iLayer=0; iLayer<8; ++iLayer) {
        AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
        fHistosManager->FillHistClass("TPCclusterMap", fValues);
      }
    }
    if(track->IsGammaLeg()) fHistosManager->FillHistClass("TrackQA_GammaLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    if(track->IsPureGammaLeg()) fHistosManager->FillHistClass("TrackQA_PureGammaLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    if(track->IsK0sLeg()) fHistosManager->FillHistClass("TrackQA_K0sLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    if(track->IsPureK0sLeg()) fHistosManager->FillHistClass("TrackQA_PureK0sLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    if(track->IsLambdaLeg()) {
      if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_LambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
      else fHistosManager->FillHistClass("TrackQA_LambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    }
    if(track->IsPureLambdaLeg()) {
      if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_PureLambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
      else fHistosManager->FillHistClass("TrackQA_PureLambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    }
    if(track->IsALambdaLeg()) {
      if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_ALambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
      else fHistosManager->FillHistClass("TrackQA_ALambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    }
    if(track->IsPureALambdaLeg()) {
      if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_PureALambdaPosLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
      else fHistosManager->FillHistClass("TrackQA_PureALambdaNegLeg_ITS_TPC_TRD_TOF_EMCAL", fValues);
    }  
  }  // end loop over tracks
    
  TClonesArray* pairList = fEvent->GetPairs();
  TIter nextPair(pairList);
  AliReducedPairInfo* pair = 0x0;
  fValues[AliReducedVarManager::kNpairsSelected] = 0.;
  for(Int_t ip=0; ip<pairList->GetEntries(); ++ip) {
    pair = (AliReducedPairInfo*)nextPair();
    
    if(!IsPairSelected(pair)) continue;
    
    TString pairTypeStr = "";
    if(pair->PairType()==0) pairTypeStr = "Offline";
    if(pair->PairType()==1) pairTypeStr = "OnTheFly";
      
    for(UShort_t iflag=0; iflag<32; ++iflag) {
      AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
      fHistosManager->FillHistClass(Form("PairQualityFlags_%s",pairTypeStr.Data()), fValues);
    }
      
    AliReducedVarManager::FillPairInfo(pair, fValues);
    switch (pair->CandidateId()) {
      case AliReducedPairInfo::kGammaConv :
        fHistosManager->FillHistClass(Form("PairQA_%sGamma",pairTypeStr.Data()),fValues);
	if(pair->IsPureV0Gamma()) fHistosManager->FillHistClass(Form("PairQA_%sPureGamma",pairTypeStr.Data()),fValues);
        break;
      case AliReducedPairInfo::kK0sToPiPi :
	fHistosManager->FillHistClass(Form("PairQA_%sK0s",pairTypeStr.Data()),fValues);
	if(pair->IsPureV0K0s()) fHistosManager->FillHistClass(Form("PairQA_%sPureK0s",pairTypeStr.Data()),fValues);
	break;
      case AliReducedPairInfo::kLambda0ToPPi :
	fHistosManager->FillHistClass(Form("PairQA_%sLambda", pairTypeStr.Data()),fValues);
	if(pair->IsPureV0Lambda()) fHistosManager->FillHistClass(Form("PairQA_%sPureLambda",pairTypeStr.Data()),fValues);
	break;
      case AliReducedPairInfo::kALambda0ToPPi :
	fHistosManager->FillHistClass(Form("PairQA_%sALambda", pairTypeStr.Data()),fValues);
	if(pair->IsPureV0ALambda()) fHistosManager->FillHistClass(Form("PairQA_%sPureALambda",pairTypeStr.Data()),fValues);
	break;
    };
    fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
  }  // end loop over pairs
    
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::Finish() {
  //
  // run stuff after the event loop
  //
}
