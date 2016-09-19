//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisJpsi2ee.h"

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
AliReducedAnalysisJpsi2ee::AliReducedAnalysisJpsi2ee(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler()),
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
AliReducedAnalysisJpsi2ee::~AliReducedAnalysisJpsi2ee() 
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
void AliReducedAnalysisJpsi2ee::AddTrackCut(AliReducedInfoCut* cut) {
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
    else { if(cut->IsSelected(track)) track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
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
Bool_t AliReducedAnalysisJpsi2ee::IsPairSelected(Float_t* values) {
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
Bool_t AliReducedAnalysisJpsi2ee::IsPairPreFilterSelected(Float_t* values) {
   //
   // apply event cuts
   //
   if(fPreFilterPairCuts.GetEntries()==0) return kTRUE;
   // loop over all the cuts and make a logical and between all cuts in the list
   for(Int_t i=0; i<fPreFilterPairCuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fPreFilterPairCuts.At(i);
      if(!cut->IsSelected(values)) return kFALSE;
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
   
   fMixingHandler->SetHistogramManager(fHistosManager);
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::Process() {
  //
  // process the current event
  //  
  if(!fEvent) return;
  AliReducedEventInfo* eventInfo = NULL;
  if(fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
  else {
     cout << "ERROR: AliReducedAnalysisJpsi2ee::Process() needs AliReducedEventInfo events" << endl;
     return;
  }
  
  AliReducedVarManager::SetEvent(fEvent);
  
  //cout << "############ E v e n t #################" << endl;
  // reset the values array
  for(Int_t i=0; i<AliReducedVarManager::kNVars; ++i) fValues[i]=-9999.;
  
  // fill event information before event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_BeforeCuts", fValues);
  //cout << "### Event CentVZERO: " << fValues[AliReducedVarManager::kCentVZERO] << endl;
  //cout << "### Event VtxZ: " << fValues[AliReducedVarManager::kVtxZ] << endl;
  //cout << "### Event VZERO-EP: " << fValues[AliReducedVarManager::kVZERORP+1] << endl;
  
  // apply event selection
  if(!IsEventSelected(fEvent)) return;
  //cout << "### Event SELECTED " << endl;
  
  // fill event info histograms after cuts
  // TODO: move this to the end of the event loop since some event info can be track counters, which are 
  //   updated during the track loop
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
  
  // select tracks
  RunTrackSelection();
    
  // Run the prefilter  
  // NOTE: Pair each track from the selected tracks list with all selected tracks in the prefilter track list
  //         If the created pair fails the pair prefilter criteria, then the selected trak is removed from the track list
  //          and further pairing
  //FillTrackHistograms("Track_NoPrefilter");
  //RunSameEventPairing("PairPrefilterSE");
  RunPrefilter();
  
  // Fill track histograms
  FillTrackHistograms();
  
  // Feed the selected tracks to the event mixing handler 
  fMixingHandler->FillEvent(&fPosTracks, &fNegTracks, fValues, AliReducedPairInfo::kJpsiToEE);
  
  // Do the same event pairing
  RunSameEventPairing();
 
  //fMixingHandler->PrintMixingLists(10);
  // Done
  if(fEventCounter%100000==0) cout << "Event no. " << fEventCounter << endl;
  fEventCounter++;
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillTrackHistograms(TString trackClass /*= "Track"*/) {
   //
   // Fill all track histograms
   //
   AliReducedTrackInfo* track=0;
   TIter nextPosTrack(&fPosTracks);
   for(Int_t i=0;i<fPosTracks.GetEntries();++i) {
      track = (AliReducedTrackInfo*)nextPosTrack();
      AliReducedVarManager::FillTrackInfo(track, fValues);
      FillTrackHistograms(track, trackClass);
   }
   TIter nextNegTrack(&fNegTracks);
   for(Int_t i=0;i<fNegTracks.GetEntries();++i) {
      track = (AliReducedTrackInfo*)nextNegTrack();
      AliReducedVarManager::FillTrackInfo(track, fValues);
      FillTrackHistograms(track, trackClass);
   }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillTrackHistograms(AliReducedTrackInfo* track, TString trackClass /*="Track"*/) {
   //
   // fill track level histograms
   //
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(track->TestFlag(icut)) fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillPairHistograms(ULong_t mask, Int_t pairType, TString pairClass /*="PairSE"*/) {
   //
   // fill pair level histograms
   // NOTE: pairType can be 0,1 or 2 corresponding to ++, +- or -- pairs
   TString typeStr[3] = {"PP", "PM", "MM"};
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(mask & (ULong_t(1)<<icut)) fHistosManager->FillHistClass(Form("%s%s_%s", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(icut)->GetName()), fValues);
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::RunTrackSelection() {
   //
   // select electron candidates and prefilter tracks
   //
   // clear the track arrays
   fPosTracks.Clear("C"); fNegTracks.Clear("C"); fPrefilterPosTracks.Clear("C"); fPrefilterNegTracks.Clear("C");
   // loop over the track list and evaluate all the track cuts
   AliReducedTrackInfo* track = 0x0;
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   //cout << "### Track selection for " << fEvent->NTracks() << " tracks"  << endl;
   for(Int_t it=0; it<fEvent->NTracks(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      //cout << "$$ track #" << it << ":  TrackId " << track->TrackId() << "; (px,py,pz): (" << track->Px() << "," << track->Py() << "," << track->Pz() << ")" << endl;
      
      AliReducedVarManager::FillTrackInfo(track, fValues);
      if(IsTrackSelected(track, fValues)) {
         //cout << " --> selected: "; AliReducedVarManager::PrintTrackFlags(track); cout << endl;
         if(track->Charge()>0) fPosTracks.Add(track);
         if(track->Charge()<0) fNegTracks.Add(track);
      }
      if(IsTrackPrefilterSelected(track, fValues)) {
         //cout << " --> prefilter selected: "; AliReducedVarManager::PrintTrackFlags(track); cout << endl;
         if(track->Charge()>0) fPrefilterPosTracks.Add(track);
         if(track->Charge()<0) fPrefilterNegTracks.Add(track);
      }
   }   // end loop over tracks
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::RunSameEventPairing(TString pairClass /*="PairSE"*/) {
   //
   // Run the same event pairing
   //
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   
   AliReducedTrackInfo* pTrack=0;
   AliReducedTrackInfo* pTrack2=0;
   AliReducedTrackInfo* nTrack=0;
   AliReducedTrackInfo* nTrack2=0;
   //cout << "### Run same event pairing for pair class " << pairClass.Data() << endl;
   //cout << "### Positive / negative tracks: " << fPosTracks.GetEntries() << " / " << fNegTracks.GetEntries() << endl;
   for(Int_t ip=0; ip<fPosTracks.GetEntries(); ++ip) {
      pTrack = (AliReducedTrackInfo*)nextPosTrack();
      //cout << "$$$$ ip #" << ip << " id: " << pTrack->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(pTrack); cout << endl;
      
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedTrackInfo*)nextNegTrack();
         //cout << "$$$$$$ in #" << in << " id: " << nTrack->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(nTrack); cout << endl;
         
         // verify that the two current tracks have at least 1 common bit
         if(!(pTrack->GetFlags() & nTrack->GetFlags())) continue;
         //cout << "$$$$$$ mask: "; AliReducedVarManager::PrintBits(pTrack->GetFlags() & nTrack->GetFlags(), fTrackCuts.GetEntries()); cout<<endl;
         AliReducedVarManager::FillPairInfo(pTrack, nTrack, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(IsPairSelected(fValues)) {
            //cout << "$$$$$$ pair is selected" << endl;
            FillPairHistograms(pTrack->GetFlags() & nTrack->GetFlags(), 1, pairClass);    // 1 is for +- pairs 
         }
      }  // end loop over negative tracks
      
      for(Int_t ip2=ip+1; ip2<fPosTracks.GetEntries(); ++ip2) {
         pTrack2 = (AliReducedTrackInfo*)fPosTracks.At(ip2);
         //cout << "$$$$$$ ip2 #" << ip2 << " id: " << pTrack2->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(pTrack2); cout << endl;
         
         // verify that the two current tracks have at least 1 common bit
         if(!(pTrack->GetFlags() & pTrack2->GetFlags())) continue;
         //cout << "$$$$$$ mask: "; AliReducedVarManager::PrintBits(pTrack->GetFlags() & pTrack2->GetFlags(), fTrackCuts.GetEntries()); cout<<endl;
         AliReducedVarManager::FillPairInfo(pTrack, pTrack2, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(IsPairSelected(fValues)) {
            //cout << "$$$$$$ pair is selected" << endl;
            FillPairHistograms(pTrack->GetFlags() & pTrack2->GetFlags(), 0, pairClass);       // 0 is for ++ pairs 
         }
      }  // end loop over positive tracks
   }  // end loop over positive tracks
   
   nextNegTrack.Reset();
   for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
      nTrack = (AliReducedTrackInfo*)nextNegTrack();
      //cout << "$$$$$$ in #" << in << " id: " << nTrack->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(nTrack); cout << endl;
      
      for(Int_t in2=in+1; in2<fNegTracks.GetEntries(); ++in2) {
         nTrack2 = (AliReducedTrackInfo*)fNegTracks.At(in2);
         //cout << "$$$$$$ in2 #" << in2 << " id: " << nTrack2->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(nTrack2); cout << endl;
         
         // verify that the two current tracks have at least 1 common bit
         if(!(nTrack->GetFlags() & nTrack2->GetFlags())) continue;
         //cout << "$$$$$$ mask: "; AliReducedVarManager::PrintBits(nTrack->GetFlags() & nTrack2->GetFlags(), fTrackCuts.GetEntries()); cout<<endl;
         AliReducedVarManager::FillPairInfo(nTrack, nTrack2, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(IsPairSelected(fValues)) {
            //cout << "$$$$$$ pair is selected" << endl;
            FillPairHistograms(nTrack->GetFlags() & nTrack2->GetFlags(), 2, pairClass);      // 2 is for -- pairs
         }
      }  // end loop over negative tracks
   }  // end loop over negative tracks
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::RunPrefilter() {
   //
   // Run the prefilter selection
   // At this point it is assumed that the track lists are filled
   //
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   TIter nextPosPrefilterTrack(&fPrefilterPosTracks);
   TIter nextNegPrefilterTrack(&fPrefilterNegTracks);
   
   // First pair the positive trackes with the prefilter selected tracks
   //cout << "##### Running the prefilter " << endl;
   AliReducedTrackInfo* track=0;
   AliReducedTrackInfo* trackPref=0;
   for(Int_t ip = 0; ip<fPosTracks.GetEntries(); ++ip) {
      track = (AliReducedTrackInfo*)nextPosTrack();
      //cout << "$$$$$$ ip #" << ip << " id: " << track->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(track); cout << endl;
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedTrackInfo*)nextPosPrefilterTrack();
         //cout << "$$$$$$ ipp #" << ipp << " id: " << trackPref->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(trackPref); cout << endl;
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            //cout << "track failed the prefilter test, breaking the prefilter tracks loop" << endl;
            break;
         }
         else {
            //cout << "pair passed the prefilter" << endl;
         }
      }  // end loop over positive prefilter tracks
      if(!track->GetFlags()) {
         fPosTracks.Remove(track); 
         //cout << "track removed from the positives list and continuing" << endl;
         continue;
      }
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedTrackInfo*)nextNegPrefilterTrack();
         //cout << "$$$$$$ ipn #" << ipn << " id: " << trackPref->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(trackPref); cout << endl;
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            //cout << "track failed the prefilter test, breaking the prefilter tracks loop" << endl;
            break;
         }
         else {
            //cout << "pair passed the prefilter" << endl;
         }
      }  // end loop over negative prefilter tracks
      
      if(!track->GetFlags()) {
         fPosTracks.Remove(track);
         //cout << "track removed from the positives list" << endl;
      }
   }  // end loop over the positive tracks

   for(Int_t in = 0; in<fNegTracks.GetEntries(); ++in) {
      track = (AliReducedTrackInfo*)nextNegTrack();
      //cout << "$$$$$$ in #" << in << " id: " << track->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(track); cout << endl;
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedTrackInfo*)nextPosPrefilterTrack();
         //cout << "$$$$$$ ipp #" << ipp << " id: " << trackPref->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(trackPref); cout << endl;
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            //cout << "track failed the prefilter test, breaking the prefilter tracks loop" << endl;
            break;
         }
      }  // end loop over positive prefilter tracks
      if(!track->GetFlags()) {
         fNegTracks.Remove(track); 
         //cout << "track removed from the negatives list" << endl;
         continue;
      }
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedTrackInfo*)nextNegPrefilterTrack();
         //cout << "$$$$$$ ipn #" << ipn << " id: " << trackPref->TrackId() << ";"; AliReducedVarManager::PrintTrackFlags(trackPref); cout << endl;
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         //cout << "$$$$$$ mass/type: " << fValues[AliReducedVarManager::kMass] << "/" << fValues[AliReducedVarManager::kPairType] << endl;
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            //cout << "track failed the prefilter test, breaking the prefilter tracks loop" << endl;
            break;
         }
      }  // end loop over negative prefilter tracks
      if(!track->GetFlags()) {
         fNegTracks.Remove(track);
         //cout << "track removed from the negatives list" << endl;
      }
   }  // end loop over the negative tracks
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::Finish() {
  //
  // run stuff after the event loop
  //
   fMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
   //cout << "*********************************************************  Finish is ending" << endl;
}
