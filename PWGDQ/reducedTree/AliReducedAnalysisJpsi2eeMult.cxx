//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisJpsi2eeMult.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisJpsi2eeMult);


//___________________________________________________________________________
AliReducedAnalysisJpsi2eeMult::AliReducedAnalysisJpsi2eeMult() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler()),
  fOptionRunMixing(kTRUE),
  fOptionRunPairing(kTRUE),
  fOptionRunOverMC(kFALSE),
  fOptionRunLikeSignPairing(kTRUE),
  fOptionLoopOverTracks(kTRUE),
  fEventCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts(),
  fPosTracks(),
  fNegTracks(),
  fPrefilterPosTracks(),
  fPrefilterNegTracks(),
  fEventCounter(0)
{
  //
  // default constructor
  //
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2eeMult::AliReducedAnalysisJpsi2eeMult(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler()),
  fOptionRunMixing(kTRUE),
  fOptionRunPairing(kTRUE),
  fOptionRunOverMC(kFALSE),
  fOptionRunLikeSignPairing(kTRUE),
  fOptionLoopOverTracks(kTRUE),
  fEventCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts(),
  fPosTracks(),
  fNegTracks(),
  fPrefilterPosTracks(),
  fPrefilterNegTracks(),
  fEventCounter(0)
{
  //
  // named constructor
  //
   fEventCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
   fPreFilterTrackCuts.SetOwner(kTRUE);
   fPairCuts.SetOwner(kTRUE);
   fPreFilterPairCuts.SetOwner(kTRUE);
   fPosTracks.SetOwner(kFALSE);
   fNegTracks.SetOwner(kFALSE);
   fPrefilterPosTracks.SetOwner(kFALSE);
   fPrefilterNegTracks.SetOwner(kFALSE);
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2eeMult::~AliReducedAnalysisJpsi2eeMult() 
{
  //
  // destructor
  //
   fEventCuts.Clear("C"); fTrackCuts.Clear("C"); fPreFilterTrackCuts.Clear("C"); fPreFilterPairCuts.Clear("C"); fPairCuts.Clear("C");
   fPosTracks.Clear("C"); fNegTracks.Clear("C"); fPrefilterPosTracks.Clear("C"); fPrefilterNegTracks.Clear("C");
   if(fHistosManager) delete fHistosManager;
   if(fMixingHandler) delete fMixingHandler;
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::AddTrackCut(AliReducedInfoCut* cut) {
   //
   // Add a new cut
   //
   fTrackCuts.Add(cut); 
   fMixingHandler->SetNParallelCuts(fMixingHandler->GetNParallelCuts()+1);
   TString histClassNames = fMixingHandler->GetHistClassNames();
   histClassNames += Form("PairMEPP_%s;", cut->GetName());
   histClassNames += Form("PairMEPM_%s;", cut->GetName());
   histClassNames += Form("PairMEMM_%s;", cut->GetName());
   fMixingHandler->SetHistClassNames(histClassNames.Data());
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
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
Bool_t AliReducedAnalysisJpsi2eeMult::IsTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fTrackCuts.GetEntries()==0) return kTRUE;
  track->ResetFlags();
  
  for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track, values)) track->SetFlag(i); }
    else { if(cut->IsSelected(track)) track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsTrackPrefilterSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
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
Bool_t AliReducedAnalysisJpsi2eeMult::IsPairSelected(Float_t* values) {
  //
  // apply event cuts
  //
  if(fPairCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
    if(!cut->IsSelected(values)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsPairPreFilterSelected(Float_t* values) {
   //
   // apply event cuts
   //
   if(fPreFilterPairCuts.GetEntries()==0) return kTRUE;
   // loop over all the cuts and make a logical OR between all cuts in the list
   for(Int_t i=0; i<fPreFilterPairCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPreFilterPairCuts.At(i);
      if(cut->IsSelected(values)) return kTRUE;
   }
   return kFALSE;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
   
   fMixingHandler->SetHistogramManager(fHistosManager);
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::Process() {
  //
  // process the current event
  //  
  if(!fEvent) return;
  AliReducedEventInfo* eventInfo = NULL;
  if(fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
  else {
     cout << "ERROR: AliReducedAnalysisJpsi2eeMult::Process() needs AliReducedEventInfo events" << endl;
     return;
  }
  if(fOptionRunOverMC) {
     if(fEventCounter%10000==0) 
        cout << "Event no. " << fEventCounter << endl;
  }
  else {
    if(fEventCounter%100000==0) 
       cout << "Event no. " << fEventCounter << endl;
  }
  fEventCounter++;
  
  AliReducedVarManager::SetEvent(fEvent);
  
  // reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: the run wise data will be updated automatically in the VarManager in case a run change is detected
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i) fValues[i]=-9999.;
  
  // fill event information before event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_BeforeCuts", fValues);
  for(UShort_t ibit=0; ibit<64; ++ibit) {
     AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
     fHistosManager->FillHistClass("EventTag_BeforeCuts", fValues);
  }
  for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_BeforeCuts", fValues);
  }
  
  
  // apply event selection
  if(!IsEventSelected(fEvent)) return;
  
  if(fOptionRunOverMC) FillMCTruthHistograms();
  
  // select tracks
  if(fOptionLoopOverTracks)
    RunTrackSelection();
    
  // Run the prefilter  
  // NOTE: Pair each track from the selected tracks list with all selected tracks in the prefilter track list
  //         If the created pair fails the pair prefilter criteria, then the selected trak is removed from the track list
  //          and further pairing
  //FillTrackHistograms("Track_BeforePrefilter");
  //RunSameEventPairing("PairPrefilterSE");
  if(fOptionLoopOverTracks)
    RunPrefilter();
  
  if(fOptionLoopOverTracks) {
    fValues[AliReducedVarManager::kNtracksPosAnalyzed] = fPosTracks.GetEntries();
    fValues[AliReducedVarManager::kNtracksNegAnalyzed] = fNegTracks.GetEntries();
    fValues[AliReducedVarManager::kNtracksAnalyzed] = fValues[AliReducedVarManager::kNtracksNegAnalyzed]+fValues[AliReducedVarManager::kNtracksPosAnalyzed];
    fValues[AliReducedVarManager::kEvAverageTPCchi2] /= (fPosTracks.GetEntries()+fNegTracks.GetEntries()>0 ? fValues[AliReducedVarManager::kNtracksAnalyzed] : 1.0); 
  }
  
  // Fill track histograms
  if(fOptionLoopOverTracks)
    FillTrackHistograms();
  
  // Feed the selected tracks to the event mixing handler 
  if(!fOptionRunOverMC && fOptionRunMixing && fOptionRunPairing)
    fMixingHandler->FillEvent(&fPosTracks, &fNegTracks, fValues, AliReducedPairInfo::kJpsiToEE);
  
  // Do the same event pairing
  if(fOptionRunPairing)
    RunSameEventPairing();
 
  // fill event info histograms after cuts
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
  for(UShort_t ibit=0; ibit<64; ++ibit) {
     AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
     fHistosManager->FillHistClass("EventTag_AfterCuts", fValues);
  }
  for(UShort_t ibit=0; ibit<64; ++ibit) {
     AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
     fHistosManager->FillHistClass("EventTriggers_AfterCuts", fValues);
  }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillTrackHistograms(TString trackClass /*= "Track"*/) {
   //
   // Fill all track histograms
   //
   for(Int_t i=0;i<36; ++i) fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+i] = 0.;
   AliReducedTrackInfo* track=0;
   TIter nextPosTrack(&fPosTracks);
   for(Int_t i=0;i<fPosTracks.GetEntries();++i) {
      track = (AliReducedTrackInfo*)nextPosTrack();
      //Int_t tpcSector = TMath::FloorNint(18.*track->Phi()/TMath::TwoPi());
      fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
      AliReducedVarManager::FillTrackInfo(track, fValues);
      FillTrackHistograms(track, Form("%s+", trackClass.Data()) );
      FillTrackHistograms(track, Form("%s", trackClass.Data()) );
   }
   TIter nextNegTrack(&fNegTracks);
   for(Int_t i=0;i<fNegTracks.GetEntries();++i) {
      track = (AliReducedTrackInfo*)nextNegTrack();
      //Int_t tpcSector = TMath::FloorNint(18.*track->Phi()/TMath::TwoPi());
      fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
      AliReducedVarManager::FillTrackInfo(track, fValues);
      FillTrackHistograms(track, Form("%s-", trackClass.Data()) );
      FillTrackHistograms(track, Form("%s", trackClass.Data()) );
      //cout << "Neg track " << i << ": "; AliReducedVarManager::PrintBits(track->Status()); cout << endl;
   }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillTrackHistograms(AliReducedTrackInfo* track, TString trackClass /*="Track"*/) {
   //
   // fill track level histograms
   //
   Bool_t isMCTruth = fOptionRunOverMC && IsMCTruth(track);
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(track->TestFlag(icut)) {
         fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         if(isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_MCTruth", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
            AliReducedVarManager::FillTrackingFlag(track, iflag, fValues);
            fHistosManager->FillHistClass(Form("%sStatusFlags_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(isMCTruth) fHistosManager->FillHistClass(Form("%sStatusFlags_%s_MCTruth", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         }
         for(Int_t iLayer=0; iLayer<6; ++iLayer) {
            AliReducedVarManager::FillITSlayerFlag(track, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sITSclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(isMCTruth) fHistosManager->FillHistClass(Form("%sITSclusterMap_%s_MCTruth", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         }
         for(Int_t iLayer=0; iLayer<8; ++iLayer) {
            AliReducedVarManager::FillTPCclusterBitFlag(track, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(isMCTruth) fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s_MCTruth", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         }
      } // end if(track->TestFlag(icut))
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillPairHistograms(ULong_t mask, Int_t pairType, TString pairClass /*="PairSE"*/, Bool_t isMCTruth /* = kFALSE*/) {
   //
   // fill pair level histograms
   // NOTE: pairType can be 0,1 or 2 corresponding to ++, +- or -- pairs
   TString typeStr[3] = {"PP", "PM", "MM"};
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(mask & (ULong_t(1)<<icut)) {
         fHistosManager->FillHistClass(Form("%s%s_%s", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(icut)->GetName()), fValues);
         if(isMCTruth && pairType==1) fHistosManager->FillHistClass(Form("%s%s_%s_MCTruth", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(icut)->GetName()), fValues);
      }
         
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::RunTrackSelection() {
   //
   // select electron candidates and prefilter tracks
   //
   // clear the track arrays
   fPosTracks.Clear("C"); fNegTracks.Clear("C"); fPrefilterPosTracks.Clear("C"); fPrefilterNegTracks.Clear("C");
   fValues[AliReducedVarManager::kEvAverageTPCchi2] = 0.0;
   
   // loop over the track list and evaluate all the track cuts
   AliReducedTrackInfo* track = 0x0;
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   Float_t nsigma = 0.;
   for(Int_t it=0; it<fEvent->NTracks(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(fOptionRunOverMC && track->IsMCTruth()) continue;
      //cout << "track " << it << ": "; AliReducedVarManager::PrintBits(track->Status()); cout << endl;
      AliReducedVarManager::FillTrackInfo(track, fValues);
      fHistosManager->FillHistClass("Track_BeforeCuts", fValues);
      for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
         //cout << "track / tracking flags :: " << track << " / "; AliReducedVarManager::PrintBits(track->Status()); cout << endl;
         AliReducedVarManager::FillTrackingFlag(track, iflag, fValues);
         fHistosManager->FillHistClass("TrackStatusFlags_BeforeCuts", fValues);
      }
      for(Int_t iLayer=0; iLayer<6; ++iLayer) {
         AliReducedVarManager::FillITSlayerFlag(track, iLayer, fValues);
         fHistosManager->FillHistClass("TrackITSclusterMap_BeforeCuts", fValues);
      }
      for(Int_t iLayer=0; iLayer<8; ++iLayer) {
         AliReducedVarManager::FillTPCclusterBitFlag(track, iLayer, fValues);
         fHistosManager->FillHistClass("TrackTPCclusterMap_BeforeCuts", fValues);
      }
      if(IsTrackSelected(track, fValues)) {
         fValues[AliReducedVarManager::kEvAverageTPCchi2] += track->TPCchi2();
         if(track->Charge()>0) fPosTracks.Add(track);
         if(track->Charge()<0) fNegTracks.Add(track);
      }
      if(IsTrackPrefilterSelected(track, fValues)) {
         if(track->Charge()>0) fPrefilterPosTracks.Add(track);
         if(track->Charge()<0) fPrefilterNegTracks.Add(track);
      }
   }   // end loop over tracks
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::RunSameEventPairing(TString pairClass /*="PairSE"*/) {
   //
   // Run the same event pairing
   //
   fValues[AliReducedVarManager::kNpairsSelected] = 0;
   
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   
   AliReducedTrackInfo* pTrack=0;
   AliReducedTrackInfo* pTrack2=0;
   AliReducedTrackInfo* nTrack=0;
   AliReducedTrackInfo* nTrack2=0;
   for(Int_t ip=0; ip<fPosTracks.GetEntries(); ++ip) {
      pTrack = (AliReducedTrackInfo*)nextPosTrack();
      
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedTrackInfo*)nextNegTrack();
         
         // verify that the two current tracks have at least 1 common bit
         if(!(pTrack->GetFlags() & nTrack->GetFlags())) continue;
         AliReducedVarManager::FillPairInfo(pTrack, nTrack, AliReducedPairInfo::kJpsiToEE, fValues);
         if(IsPairSelected(fValues)) {
            FillPairHistograms(pTrack->GetFlags() & nTrack->GetFlags(), 1, pairClass, fOptionRunOverMC && IsMCTruth(pTrack, nTrack));    // 1 is for +- pairs 
            fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
         }
      }  // end loop over negative tracks
      
      if(fOptionRunLikeSignPairing) {
         for(Int_t ip2=ip+1; ip2<fPosTracks.GetEntries(); ++ip2) {
            pTrack2 = (AliReducedTrackInfo*)fPosTracks.At(ip2);
         
            // verify that the two current tracks have at least 1 common bit
            if(!(pTrack->GetFlags() & pTrack2->GetFlags())) continue;
            AliReducedVarManager::FillPairInfo(pTrack, pTrack2, AliReducedPairInfo::kJpsiToEE, fValues);
            if(IsPairSelected(fValues)) {
               FillPairHistograms(pTrack->GetFlags() & pTrack2->GetFlags(), 0, pairClass);       // 0 is for ++ pairs 
               fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
            }
         }  // end loop over positive tracks
      }
   }  // end loop over positive tracks
   
   if(fOptionRunLikeSignPairing) {
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedTrackInfo*)nextNegTrack();
      
         for(Int_t in2=in+1; in2<fNegTracks.GetEntries(); ++in2) {
            nTrack2 = (AliReducedTrackInfo*)fNegTracks.At(in2);
         
            // verify that the two current tracks have at least 1 common bit
            if(!(nTrack->GetFlags() & nTrack2->GetFlags())) continue;
            AliReducedVarManager::FillPairInfo(nTrack, nTrack2, AliReducedPairInfo::kJpsiToEE, fValues);
            if(IsPairSelected(fValues)) {
               FillPairHistograms(nTrack->GetFlags() & nTrack2->GetFlags(), 2, pairClass);      // 2 is for -- pairs
               fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
            }
         }  // end loop over negative tracks
      }  // end loop over negative tracks
   }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::RunPrefilter() {
   //
   // Run the prefilter selection
   // At this point it is assumed that the track lists are filled
   //
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   TIter nextPosPrefilterTrack(&fPrefilterPosTracks);
   TIter nextNegPrefilterTrack(&fPrefilterNegTracks);
   
   // First pair the positive trackes with the prefilter selected tracks
   AliReducedTrackInfo* track=0;
   AliReducedTrackInfo* trackPref=0;
   for(Int_t ip = 0; ip<fPosTracks.GetEntries(); ++ip) {
      track = (AliReducedTrackInfo*)nextPosTrack();
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedTrackInfo*)nextPosPrefilterTrack();
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over positive prefilter tracks
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedTrackInfo*)nextNegPrefilterTrack();
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over negative prefilter tracks
   }  // end loop over the positive tracks

   for(Int_t in = 0; in<fNegTracks.GetEntries(); ++in) {
      track = (AliReducedTrackInfo*)nextNegTrack();
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedTrackInfo*)nextPosPrefilterTrack();
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over positive prefilter tracks
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedTrackInfo*)nextNegPrefilterTrack();
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over negative prefilter tracks
   }  // end loop over the negative tracks

   // remove tracks
   nextPosTrack.Reset();
   for(Int_t ip = fPosTracks.GetEntries()-1 ; ip >= 0; --ip) {
     track = (AliReducedTrackInfo*)nextPosTrack();
     if(!track->GetFlags()) fPosTracks.Remove(track);
   }
  nextNegTrack.Reset();
  for(Int_t ip = fNegTracks.GetEntries()-1 ; ip >= 0; --ip) {
    track = (AliReducedTrackInfo*)nextNegTrack();
    if(!track->GetFlags()) fNegTracks.Remove(track);
  }

}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::Finish() {
  //
  // run stuff after the event loop
  //
   if(fOptionRunMixing && !fOptionRunOverMC)
     fMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsMCTruth(AliReducedTrackInfo* track) {
   //
   // check whether the track is an electron from a J/psi decay
   //
   if(TMath::Abs(track->MCPdg(0)) != 11) return kFALSE;
   if(TMath::Abs(track->MCPdg(1)) != 443) return kFALSE;
   if(track->MCPdg(2) != -9999) return kFALSE;
   return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeMult::IsMCTruth(AliReducedTrackInfo* ptrack, AliReducedTrackInfo* ntrack) {
   //
   // check whether the tracks are electrons from a common J/psi decay
   //
   if(TMath::Abs(ptrack->MCPdg(0)) != 11) return kFALSE;
   if(TMath::Abs(ntrack->MCPdg(0)) != 11) return kFALSE;
   if(TMath::Abs(ptrack->MCPdg(1)) != 443) return kFALSE;
   if(TMath::Abs(ntrack->MCPdg(1)) != 443) return kFALSE;
   if(ptrack->MCPdg(2) != -9999) return kFALSE;
   if(ntrack->MCPdg(2) != -9999) return kFALSE;   
   if(TMath::Abs(ptrack->MCLabel(1)) != TMath::Abs(ntrack->MCLabel(1))) return kFALSE;
   return kTRUE;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FillMCTruthHistograms() {
  //
  // fill histograms with pure signal
  //   
  AliReducedTrackInfo* track = 0x0;
  Int_t leg1Id = -1;
  Int_t leg2Id = -1;
  AliReducedTrackInfo* leg1=0x0;
  AliReducedTrackInfo* leg2=0x0;
  TClonesArray* trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  Float_t nsigma = 0.;
  for(Int_t it=0; it<fEvent->NTracks(); ++it) {
     track = (AliReducedTrackInfo*)nextTrack();
     if(!track->IsMCTruth()) continue;
     
     if(track->MCPdg(0)==443 && TMath::Abs(track->Rapidity(3.1))<0.9) {       // TODO: use the correct PDG mass and dynamic kinematic selection
       leg1Id = -1; leg2Id = -1;
       FindJpsiTruthLegs(track, leg1Id, leg2Id);
       leg1 = (leg1Id>-1 ? (AliReducedTrackInfo*)fEvent->GetTrack(leg1Id) : 0x0);
       leg2 = (leg2Id>-1 ? (AliReducedTrackInfo*)fEvent->GetTrack(leg2Id) : 0x0);
       AliReducedVarManager::FillMCTruthInfo(track, fValues, leg1, leg2);
       fHistosManager->FillHistClass("MCTruth_BeforeSelection", fValues);
       if(!leg1) continue;
       if(!leg2) continue;
       if(TMath::Abs(leg1->EtaMC())>0.9) continue;                       // TODO: use dynamic kinematic cut on legs
       if(TMath::Abs(leg2->EtaMC())>0.9) continue;
       if(leg1->PtMC()<1.0) continue;
       if(leg2->PtMC()<1.0) continue;
       fHistosManager->FillHistClass("MCTruth_AfterSelection", fValues);
     }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeMult::FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1, Int_t& leg2) {
   //
   // find the jpsi legs in the list of pure MC truth particles
   //
   Int_t mLabel = mother->MCLabel(0);
   AliReducedTrackInfo* track=0x0;
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   Int_t legsFound = 0;
   for(Int_t it=0; it<fEvent->NTracks(); ++it) {
      if(legsFound==2) return;
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCTruth()) continue;
      if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==11) {
         legsFound += 1;
         if(legsFound==1) leg1 = it;
         if(legsFound==2) leg2 = it;
         //if(TMath::Abs(track->EtaMC())>0.9) return kFALSE;                       // TODO: use dynamic kinematic cut on legs
         //if(track->PtMC()<1.0) return kFALSE;
      }
   }
   return;
}

