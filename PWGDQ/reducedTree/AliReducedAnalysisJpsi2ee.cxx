//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisJpsi2ee.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>
#include <TList.h>
#include <TRandom.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedCaloClusterTrackMatcher.h"
#include "AliHistogramManager.h"
#include "AliReducedTrackCut.h"

ClassImp(AliReducedAnalysisJpsi2ee);


//___________________________________________________________________________
AliReducedAnalysisJpsi2ee::AliReducedAnalysisJpsi2ee() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
  fClusterTrackMatcher(0x0),
  fOptionRunMixing(kTRUE),
  fOptionRunPairing(kTRUE),
  fOptionRunOverMC(kFALSE),
  fOptionRunLikeSignPairing(kTRUE),
  fOptionLoopOverTracks(kTRUE),
  fOptionRunPrefilter(kTRUE),
  fOptionStoreJpsiCandidates(kFALSE),
  fFillCaloClusterHistograms(kFALSE),
  fEventCuts(),
  fClusterCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts(),
  fClusters(),
  fPosTracks(),
  fNegTracks(),
  fPrefilterPosTracks(),
  fPrefilterNegTracks(),
  fJpsiCandidates(),
  fLegCandidatesMCcuts(),
  fLegCandidatesMCcuts_RequestSameMother(),
  fJpsiMotherMCcuts(),
  fJpsiElectronMCcuts(),
  fClusterTrackMatcherHistograms(0x0),
  fClusterTrackMatcherMultipleMatchesBefore(0x0),
  fClusterTrackMatcherMultipleMatchesAfter(0x0),
  fSkipMCEvent(kFALSE),
  fMCJpsiPtWeights(0x0)
{
  //
  // default constructor
  //
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2ee::AliReducedAnalysisJpsi2ee(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
  fClusterTrackMatcher(0x0),
  fOptionRunMixing(kTRUE),
  fOptionRunPairing(kTRUE),
  fOptionRunOverMC(kFALSE),
  fOptionRunLikeSignPairing(kTRUE),
  fOptionLoopOverTracks(kTRUE),
  fOptionRunPrefilter(kTRUE),
  fOptionStoreJpsiCandidates(kFALSE),
  fFillCaloClusterHistograms(kFALSE),
  fEventCuts(),
  fClusterCuts(),
  fTrackCuts(),
  fPreFilterTrackCuts(),
  fPairCuts(),
  fPreFilterPairCuts(),
  fClusters(),
  fPosTracks(),
  fNegTracks(),
  fPrefilterPosTracks(),
  fPrefilterNegTracks(),
  fJpsiCandidates(),
  fLegCandidatesMCcuts(),
  fLegCandidatesMCcuts_RequestSameMother(),
  fJpsiMotherMCcuts(),
  fJpsiElectronMCcuts(),
  fClusterTrackMatcherHistograms(0x0),
  fClusterTrackMatcherMultipleMatchesBefore(0x0),
  fClusterTrackMatcherMultipleMatchesAfter(0x0),
  fSkipMCEvent(kFALSE),
  fMCJpsiPtWeights(0x0)
{
  //
  // named constructor
  //
   fEventCuts.SetOwner(kTRUE);
   fClusterCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
   fPreFilterTrackCuts.SetOwner(kTRUE);
   fPairCuts.SetOwner(kTRUE);
   fPreFilterPairCuts.SetOwner(kTRUE);
   fClusters.SetOwner(kFALSE);
   fPosTracks.SetOwner(kFALSE);
   fNegTracks.SetOwner(kFALSE);
   fPrefilterPosTracks.SetOwner(kFALSE);
   fPrefilterNegTracks.SetOwner(kFALSE);
   fJpsiCandidates.SetOwner(kTRUE);
   fLegCandidatesMCcuts.SetOwner(kTRUE);
   fJpsiMotherMCcuts.SetOwner(kTRUE);
   fJpsiElectronMCcuts.SetOwner(kTRUE);
   for(Int_t i=0;i<32;++i) fLegCandidatesMCcuts_RequestSameMother[i] = kTRUE;
}


//___________________________________________________________________________
AliReducedAnalysisJpsi2ee::~AliReducedAnalysisJpsi2ee() 
{
  //
  // destructor
  //
   fEventCuts.Clear("C");
   fClusterCuts.Clear("C");
   fTrackCuts.Clear("C");
   fPreFilterTrackCuts.Clear("C");
   fPreFilterPairCuts.Clear("C");
   fPairCuts.Clear("C");
   fClusters.Clear("C");
   fPosTracks.Clear("C");
   fNegTracks.Clear("C");
   fPrefilterPosTracks.Clear("C");
   fPrefilterNegTracks.Clear("C");
   fJpsiCandidates.Clear("C");
   if(fHistosManager) delete fHistosManager;
   if(fMixingHandler) delete fMixingHandler;
   if (fClusterTrackMatcher) delete fClusterTrackMatcher;
   if (fClusterTrackMatcherHistograms) delete fClusterTrackMatcherHistograms;
   if (fClusterTrackMatcherMultipleMatchesBefore) delete fClusterTrackMatcherMultipleMatchesBefore;
   if (fClusterTrackMatcherMultipleMatchesAfter) delete fClusterTrackMatcherMultipleMatchesAfter;
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::AddTrackCut(AliReducedInfoCut* cut) {
   //
   // Add a new cut
   //
   fTrackCuts.Add(cut);
   fMixingHandler->SetNParallelCuts(fMixingHandler->GetNParallelCuts()+1);
   TString histClassNames = fMixingHandler->GetHistClassNames();
   if (fPairCuts.GetEntries()>1) {
      histClassNames = "";
      for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
         for (Int_t iTrackCut=0; iTrackCut<fTrackCuts.GetEntries(); iTrackCut++) {
            histClassNames += Form("PairMEPP_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
            histClassNames += Form("PairMEPM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
            histClassNames += Form("PairMEMM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
         }
      }
   } else {
      histClassNames += Form("PairMEPP_%s;", cut->GetName());
      histClassNames += Form("PairMEPM_%s;", cut->GetName());
      histClassNames += Form("PairMEMM_%s;", cut->GetName());
   }
   fMixingHandler->SetHistClassNames(histClassNames.Data());
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::AddPairCut(AliReducedInfoCut* cut) {
   //
   // Add a new cut
   //
   fPairCuts.Add(cut);
   fMixingHandler->SetNParallelPairCuts(fMixingHandler->GetNParallelPairCuts()+1);
   if (fPairCuts.GetEntries()>1) {
      TString histClassNamesNew = "";
      for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
         for (Int_t iTrackCut=0; iTrackCut<fTrackCuts.GetEntries(); iTrackCut++) {
            histClassNamesNew += Form("PairMEPP_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
            histClassNamesNew += Form("PairMEPM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
            histClassNamesNew += Form("PairMEMM_%s_%s;", fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName());
          }
      }
      fMixingHandler->SetHistClassNames(histClassNamesNew.Data());
   }
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
Bool_t AliReducedAnalysisJpsi2ee::IsClusterSelected(AliReducedCaloClusterInfo* cluster, Float_t* values/*=0x0*/) {
  //
  // apply cluster cuts
  //
  if (fClusterCuts.GetEntries()==0) return kTRUE;
  cluster->ResetFlags();
  for (Int_t i=0; i<fClusterCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fClusterCuts.At(i);
    if (values) { if (cut->IsSelected(cluster, values)) cluster->SetFlag(i); }
    else { if (cut->IsSelected(cluster)) cluster->SetFlag(i); }
  }
  return (cluster->GetFlags()>0 ? kTRUE : kFALSE);
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
   // apply track prefilter cuts
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
ULong_t AliReducedAnalysisJpsi2ee::IsPairSelected(Float_t* values) {
  //
  // apply pair cuts
  //
  if(fPairCuts.GetEntries()==0) return 1;

  ULong_t mask = 0;
  for (Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
     AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
     if (cut->IsSelected(values)) mask|=(ULong_t(1)<<i);
  }
  return mask;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2ee::IsPairPreFilterSelected(Float_t* values) {
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
void AliReducedAnalysisJpsi2ee::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   // make sure variables needed to create jpsi candidate objects are filled
   if(fOptionStoreJpsiCandidates) {
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPt);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPhi);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kEta);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kMass);
      AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPairType);
   }
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
   
   fMixingHandler->SetHistogramManager(fHistosManager);

  if (fClusterTrackMatcher) {
    fClusterTrackMatcherMultipleMatchesBefore = new TH1I("multipleCounts_beforeMatching", "mulitple counts of matched cluster IDs beofore matching", 50, 0.5, 50.5);
    fClusterTrackMatcherMultipleMatchesAfter = new TH1I("multipleCounts_afterMatching", "mulitple counts of matched cluster IDs after matching", 50, 0.5, 50.5);
    fClusterTrackMatcherHistograms = new TList();
    fClusterTrackMatcherHistograms->SetOwner();
    fClusterTrackMatcherHistograms->SetName("ClusterTrackMatcherHistograms");
    fClusterTrackMatcherHistograms->Add(fClusterTrackMatcherMultipleMatchesBefore);
    fClusterTrackMatcherHistograms->Add(fClusterTrackMatcherMultipleMatchesAfter);
    fHistosManager->AddToOutputList(fClusterTrackMatcherHistograms);
  }
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
  if(fOptionRunOverMC && fEventCounter%10000==0)  cout << "Event no. " << fEventCounter << endl;
  else if(fEventCounter%10000==0)                 cout << "Event no. " << fEventCounter << endl;
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
  if(!IsEventSelected(fEvent, fValues)) return;
  
  // fill calorimeter info histograms
  if(fFillCaloClusterHistograms) {
    RunClusterSelection();
    FillClusterHistograms();
  }

  // fill MC truth histograms
  if(fOptionRunOverMC) {
     fSkipMCEvent = kFALSE;
     FillMCTruthHistograms();
     if(fSkipMCEvent) return;
  }
  
  // select tracks
  if(fOptionLoopOverTracks)
    RunTrackSelection();
    
  // Run the prefilter  
  // NOTE: Pair each track from the selected tracks list with all selected tracks in the prefilter track list
  //         If the created pair fails the pair prefilter criteria, then the selected trak is removed from the track list
  //          and further pairing
  //FillTrackHistograms("Track_BeforePrefilter");
  //RunSameEventPairing("PairPrefilterSE");
  if(fOptionLoopOverTracks && fOptionRunPrefilter)
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
  if(fOptionRunMixing)
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
  for(UShort_t ich=0; ich<64; ++ich) {
     AliReducedVarManager::FillV0Channel(ich, fValues);
     fHistosManager->FillHistClass("V0Channels", fValues);
  }
  
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillTrackHistograms(TString trackClass /*= "Track"*/) {
   //
   // Fill all track histograms
   //
   for(Int_t i=0;i<36; ++i) fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+i] = 0.;
   AliReducedBaseTrack* track=0;
   if (fClusterTrackMatcher) {
     fClusterTrackMatcher->ClearMatchedClusterIDsBefore();
     fClusterTrackMatcher->ClearMatchedClusterIDsAfter();
   }
   TIter nextPosTrack(&fPosTracks);
   for(Int_t i=0;i<fPosTracks.GetEntries();++i) {
      track = (AliReducedBaseTrack*)nextPosTrack();
      //Int_t tpcSector = TMath::FloorNint(18.*track->Phi()/TMath::TwoPi());
      fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
      // reset track variables
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
      
      AliReducedVarManager::FillTrackInfo(track, fValues);
      if (fClusterCuts.GetEntries())  AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, &fClusters, fClusterTrackMatcher);
      else                            AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, NULL, fClusterTrackMatcher);
      FillTrackHistograms(track, trackClass);
   }
   if (fClusterTrackMatcher) {
     fClusterTrackMatcher->FillMultipleMatchesHistogram(fClusterTrackMatcherMultipleMatchesBefore, fClusterTrackMatcher->GetMatchedClusterIDsBefore());
     fClusterTrackMatcher->FillMultipleMatchesHistogram(fClusterTrackMatcherMultipleMatchesAfter, fClusterTrackMatcher->GetMatchedClusterIDsAfter());
     fClusterTrackMatcher->ClearMatchedClusterIDsBefore();
     fClusterTrackMatcher->ClearMatchedClusterIDsAfter();
   }
   TIter nextNegTrack(&fNegTracks);
   for(Int_t i=0;i<fNegTracks.GetEntries();++i) {
      track = (AliReducedBaseTrack*)nextNegTrack();
      //Int_t tpcSector = TMath::FloorNint(18.*track->Phi()/TMath::TwoPi());
      fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
      // reset track variables
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
      
      AliReducedVarManager::FillTrackInfo(track, fValues);
      if (fClusterCuts.GetEntries())  AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, &fClusters, fClusterTrackMatcher);
      else                            AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, NULL, fClusterTrackMatcher);
      FillTrackHistograms(track, trackClass);
   }
   if (fClusterTrackMatcher) {
     fClusterTrackMatcher->FillMultipleMatchesHistogram(fClusterTrackMatcherMultipleMatchesBefore, fClusterTrackMatcher->GetMatchedClusterIDsBefore());
     fClusterTrackMatcher->FillMultipleMatchesHistogram(fClusterTrackMatcherMultipleMatchesAfter, fClusterTrackMatcher->GetMatchedClusterIDsAfter());
   }
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillTrackHistograms(AliReducedBaseTrack* track, TString trackClass /*="Track"*/) {
   //
   // fill track level histograms
   //
   UInt_t mcDecisionMap = 0;
   if(fOptionRunOverMC) mcDecisionMap = CheckReconstructedLegMCTruth(track);      
   
   for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
      if(track->TestFlag(icut)) {
         fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
         if(mcDecisionMap) {       // Fill histograms for tracks identified as MC truth
            for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
               if(mcDecisionMap & (UInt_t(1)<<iMC))
                  fHistosManager->FillHistClass(Form("%s_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(), 
                                                                  fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
            }
         }
         
         if(track->IsA() != AliReducedTrackInfo::Class()) continue;
         
         AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
         if(!trackInfo) continue;
         
         for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
            AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
            fHistosManager->FillHistClass(Form("%sStatusFlags_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sStatusFlags_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                                     fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
         for(UInt_t iflag=0; iflag<64; ++iflag) {
            AliReducedVarManager::FillTrackQualityFlag(trackInfo, iflag, fValues);
            fHistosManager->FillHistClass(Form("%sQualityFlags_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sQualityFlags_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                        fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
         for(Int_t iLayer=0; iLayer<6; ++iLayer) {
            AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sITSclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sITSclusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                                     fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
            AliReducedVarManager::FillITSsharedLayerFlag(trackInfo, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sITSsharedClusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sITSsharedClusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                        fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
         for(Int_t iLayer=0; iLayer<8; ++iLayer) {
            AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
            fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
            if(mcDecisionMap) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisionMap & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                                     fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
      } // end if(track->TestFlag(icut))
   }  // end loop over cuts
}


//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillPairHistograms(ULong_t trackMask, ULong_t pairMask, Int_t pairType, TString pairClass /*="PairSE"*/, UInt_t mcDecisions /* = 0*/) {
   //
   // fill pair level histograms
   // NOTE: pairType can be 0,1 or 2 corresponding to ++, +- or -- pairs
   TString typeStr[3] = {"PP", "PM", "MM"};
   if (fPairCuts.GetEntries()>1) {
      for(Int_t iTrackCut=0; iTrackCut<fTrackCuts.GetEntries(); ++iTrackCut) {
         for(Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); ++iPairCut) {
            if((trackMask & (ULong_t(1)<<iTrackCut)) && (pairMask & (ULong_t(1)<<iPairCut))) {
               fHistosManager->FillHistClass(Form("%s%s_%s_%s", pairClass.Data(), typeStr[pairType].Data(),
                                                  fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName()), fValues);
               if(mcDecisions && pairType==1) {
                  for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                     if(mcDecisions & (UInt_t(1)<<iMC))
                        fHistosManager->FillHistClass(Form("%s%s_%s_%s_%s", pairClass.Data(), typeStr[pairType].Data(),
                                                           fTrackCuts.At(iTrackCut)->GetName(), fPairCuts.At(iPairCut)->GetName(),
                                                           fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
                  }
               }
            }
         }
      }
   } else {
      for(Int_t iTrackCut=0; iTrackCut<fTrackCuts.GetEntries(); ++iTrackCut) {
         if(trackMask & (ULong_t(1)<<iTrackCut)) {
            fHistosManager->FillHistClass(Form("%s%s_%s", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(iTrackCut)->GetName()), fValues);
            if(mcDecisions && pairType==1) {
               for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                  if(mcDecisions & (UInt_t(1)<<iMC))
                     fHistosManager->FillHistClass(Form("%s%s_%s_%s", pairClass.Data(), typeStr[pairType].Data(), fTrackCuts.At(iTrackCut)->GetName(), fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
               }
            }
         }
      }  // end loop over cuts
   }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillClusterHistograms(TString clusterClass/*="CaloCluster"*/) {
  //
  // fill cluster histograms
  //
  AliReducedCaloClusterInfo* cluster = NULL;
  TIter nextCluster(&fClusters);
  for (Int_t i=0; i<fClusters.GetEntries(); ++i) {
    cluster = (AliReducedCaloClusterInfo*)nextCluster();
    for (Int_t i=AliReducedVarManager::kEMCALclusterEnergy; i<=AliReducedVarManager::kNEMCALvars; ++i) fValues[i] = -9999.;
    AliReducedVarManager::FillCaloClusterInfo(cluster, fValues);
    FillClusterHistograms(cluster, clusterClass);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillClusterHistograms(AliReducedCaloClusterInfo* cluster, TString clusterClass/*="CaloCluster"*/) {
  //
  // fill cluster histograms
  //
  for (Int_t icut=0; icut<fClusterCuts.GetEntries(); ++icut) {
    if (cluster->TestFlag(icut)) fHistosManager->FillHistClass(Form("%s_%s", clusterClass.Data(), fClusterCuts.At(icut)->GetName()), fValues);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::RunClusterSelection() {
  //
  // select cluster
  //
  fClusters.Clear("C");

  if (fEvent->IsA() == AliReducedBaseEvent::Class()) return;
  Int_t nCaloCluster = ((AliReducedEventInfo*)fEvent)->GetNCaloClusters();
  if (!nCaloCluster) return;

  AliReducedCaloClusterInfo* cluster = NULL;
  for (Int_t icl=0; icl<nCaloCluster; ++icl) {
    cluster = ((AliReducedEventInfo*)fEvent)->GetCaloCluster(icl);

    for (Int_t i=AliReducedVarManager::kEMCALclusterEnergy; i<=AliReducedVarManager::kNEMCALvars; ++i) fValues[i] = -9999.;

    AliReducedVarManager::FillCaloClusterInfo(cluster, fValues);
    fHistosManager->FillHistClass("CaloCluster_BeforeCuts", fValues);

    if (IsClusterSelected(cluster, fValues)) fClusters.Add(cluster);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::RunTrackSelection() {
   //
   // select electron candidates and prefilter tracks
   //
   // clear the track arrays
   fPosTracks.Clear("C"); fNegTracks.Clear("C"); fPrefilterPosTracks.Clear("C"); fPrefilterNegTracks.Clear("C");
   fValues[AliReducedVarManager::kEvAverageTPCchi2] = 0.0;
   
   // loop over the track list(s) and evaluate all the track cuts
   LoopOverTracks(1);      // first array
   LoopOverTracks(2);      // second array (if used)
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::LoopOverTracks(Int_t arrayOption /*=1*/) {
   //
   // Loop over a given track array, apply cuts and add selected tracks to arrays
   //
   AliReducedBaseTrack* track = 0x0;
   TClonesArray* trackList = (arrayOption==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if (!trackList) return;

   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedBaseTrack*)nextTrack();
      // do not loop over pure MC truth tracks 
      // NOTE: this can be also handled via AliReducedTrackCut::SetRejectPureMC()
      if(fOptionRunOverMC && track->IsMCTruth()) continue;     
      // reset track variables
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;

      AliReducedVarManager::FillTrackInfo(track, fValues);
      if (fClusterCuts.GetEntries())  AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, &fClusters, fClusterTrackMatcher);
      else                            AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, NULL, fClusterTrackMatcher);
      fHistosManager->FillHistClass("Track_BeforeCuts", fValues);
      
      if(track->IsA() == AliReducedTrackInfo::Class()) {
         AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
         if(trackInfo) {
            for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
               AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
               fHistosManager->FillHistClass("TrackStatusFlags_BeforeCuts", fValues);
            }
            for(UInt_t iflag=0; iflag<64; ++iflag) {
               AliReducedVarManager::FillTrackQualityFlag(trackInfo, iflag, fValues);
               fHistosManager->FillHistClass("TrackQualityFlags_BeforeCuts", fValues);
            }
            for(Int_t iLayer=0; iLayer<6; ++iLayer) {
               AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("TrackITSclusterMap_BeforeCuts", fValues);
               AliReducedVarManager::FillITSsharedLayerFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("TrackITSsharedClusterMap_BeforeCuts", fValues);
            }
            for(Int_t iLayer=0; iLayer<8; ++iLayer) {
               AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("TrackTPCclusterMap_BeforeCuts", fValues);
            }
         }
      }
      
      if(IsTrackSelected(track, fValues)) {
         if(track->Charge()>0) fPosTracks.Add(track);
         if(track->Charge()<0) fNegTracks.Add(track);
         
         if(track->IsA() == AliReducedTrackInfo::Class())
            fValues[AliReducedVarManager::kEvAverageTPCchi2] += ((AliReducedTrackInfo*)track)->TPCchi2();
      }
      if(IsTrackPrefilterSelected(track, fValues)) {
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
   if(fOptionStoreJpsiCandidates) fJpsiCandidates.Clear("C");
   fValues[AliReducedVarManager::kNpairsSelected] = 0;
   
   TIter nextPosTrack(&fPosTracks);
   TIter nextNegTrack(&fNegTracks);
   
   AliReducedBaseTrack* pTrack=0;
   AliReducedBaseTrack* pTrack2=0;
   AliReducedBaseTrack* nTrack=0;
   AliReducedBaseTrack* nTrack2=0;
   for(Int_t ip=0; ip<fPosTracks.GetEntries(); ++ip) {
      pTrack = (AliReducedBaseTrack*)nextPosTrack();
      
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedBaseTrack*)nextNegTrack();
         
         // verify that the two current tracks have at least 1 common bit
         if(!(pTrack->GetFlags() & nTrack->GetFlags())) continue;
         AliReducedVarManager::FillPairInfo(pTrack, nTrack, AliReducedPairInfo::kJpsiToEE, fValues);

         ULong_t pairCutMask = IsPairSelected(fValues);
         if(pairCutMask) {
            FillPairHistograms(pTrack->GetFlags() & nTrack->GetFlags(), pairCutMask, 1, pairClass, (fOptionRunOverMC ? CheckReconstructedLegMCTruth(pTrack, nTrack) : 0));    // 1 is for +- pairs
            fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
            if(fOptionStoreJpsiCandidates) {
               AliReducedPairInfo* pair = new AliReducedPairInfo();
               pair->SetFlags(pTrack->GetFlags() & nTrack->GetFlags());
               pair->SetQualityFlags(pairCutMask);
               pair->PtPhiEta(fValues[AliReducedVarManager::kPt], fValues[AliReducedVarManager::kPhi], fValues[AliReducedVarManager::kEta]);
               pair->SetMass(fValues[AliReducedVarManager::kMass]);
               pair->CandidateId(AliReducedPairInfo::kJpsiToEE);
               pair->PairType(1);
               pair->SetLegIds(pTrack->TrackId(), nTrack->TrackId());
               pair->SetPseudoProper(fValues[AliReducedVarManager::kPseudoProperDecayTime]);
               pair->PairTypeSPD(fValues[AliReducedVarManager::kPairTypeSPD]);
               fJpsiCandidates.Add(pair);
            }
         }
      }  // end loop over negative tracks
      
      if(fOptionRunLikeSignPairing) {
         for(Int_t ip2=ip+1; ip2<fPosTracks.GetEntries(); ++ip2) {
            pTrack2 = (AliReducedBaseTrack*)fPosTracks.At(ip2);
         
            // verify that the two current tracks have at least 1 common bit
            if(!(pTrack->GetFlags() & pTrack2->GetFlags())) continue;
            AliReducedVarManager::FillPairInfo(pTrack, pTrack2, AliReducedPairInfo::kJpsiToEE, fValues);

            ULong_t pairCutMask = IsPairSelected(fValues);
            if(pairCutMask) {
               FillPairHistograms(pTrack->GetFlags() & pTrack2->GetFlags(), pairCutMask, 0, pairClass);       // 0 is for ++ pairs
               fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
               if(fOptionStoreJpsiCandidates) {
                  AliReducedPairInfo* pair = new AliReducedPairInfo();
                  pair->SetFlags(pTrack->GetFlags() & pTrack2->GetFlags());
                  pair->SetQualityFlags(pairCutMask);
                  pair->PtPhiEta(fValues[AliReducedVarManager::kPt], fValues[AliReducedVarManager::kPhi], fValues[AliReducedVarManager::kEta]);
                  pair->SetMass(fValues[AliReducedVarManager::kMass]);
                  pair->CandidateId(AliReducedPairInfo::kJpsiToEE);
                  pair->PairType(0);
                  pair->SetLegIds(pTrack->TrackId(), pTrack2->TrackId());
                  pair->SetPseudoProper(fValues[AliReducedVarManager::kPseudoProperDecayTime]);
                  pair->PairTypeSPD(fValues[AliReducedVarManager::kPairTypeSPD]);
                  fJpsiCandidates.Add(pair);
               }
            }
         }  // end loop over positive tracks
      }
   }  // end loop over positive tracks
   
   if(fOptionRunLikeSignPairing) {
      nextNegTrack.Reset();
      for(Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
         nTrack = (AliReducedBaseTrack*)nextNegTrack();
      
         for(Int_t in2=in+1; in2<fNegTracks.GetEntries(); ++in2) {
            nTrack2 = (AliReducedBaseTrack*)fNegTracks.At(in2);
         
            // verify that the two current tracks have at least 1 common bit
            if(!(nTrack->GetFlags() & nTrack2->GetFlags())) continue;
            AliReducedVarManager::FillPairInfo(nTrack, nTrack2, AliReducedPairInfo::kJpsiToEE, fValues);

            ULong_t pairCutMask = IsPairSelected(fValues);
            if(pairCutMask) {
               FillPairHistograms(nTrack->GetFlags() & nTrack2->GetFlags(), pairCutMask, 2, pairClass);      // 2 is for -- pairs
               fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
               if(fOptionStoreJpsiCandidates) {
                  AliReducedPairInfo* pair = new AliReducedPairInfo();
                  pair->SetFlags(nTrack->GetFlags() & nTrack2->GetFlags());
                  pair->SetQualityFlags(pairCutMask);
                  pair->PtPhiEta(fValues[AliReducedVarManager::kPt], fValues[AliReducedVarManager::kPhi], fValues[AliReducedVarManager::kEta]);
                  pair->SetMass(fValues[AliReducedVarManager::kMass]);
                  pair->CandidateId(AliReducedPairInfo::kJpsiToEE);
                  pair->PairType(2);
                  pair->SetLegIds(nTrack->TrackId(), nTrack2->TrackId());
                  pair->SetPseudoProper(fValues[AliReducedVarManager::kPseudoProperDecayTime]);
                  pair->PairTypeSPD(fValues[AliReducedVarManager::kPairTypeSPD]);
                  fJpsiCandidates.Add(pair);
               }
            }
         }  // end loop over negative tracks
      }  // end loop over negative tracks
   }
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
   AliReducedBaseTrack* track=0;
   AliReducedBaseTrack* trackPref=0;
   for(Int_t ip = 0; ip<fPosTracks.GetEntries(); ++ip) {
      track = (AliReducedBaseTrack*)nextPosTrack();
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedBaseTrack*)nextPosPrefilterTrack();
         
         if(track->TrackId()==trackPref->TrackId()) continue;       // avoid self-pairing
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over positive prefilter tracks
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedBaseTrack*)nextNegPrefilterTrack();
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over negative prefilter tracks
   }  // end loop over the positive tracks

   for(Int_t in = 0; in<fNegTracks.GetEntries(); ++in) {
      track = (AliReducedBaseTrack*)nextNegTrack();
      
      nextPosPrefilterTrack.Reset();
      for(Int_t ipp = 0; ipp<fPrefilterPosTracks.GetEntries(); ++ipp) {
         trackPref = (AliReducedBaseTrack*)nextPosPrefilterTrack();
         
         AliReducedVarManager::FillPairInfo(track, trackPref, AliReducedPairInfo::kJpsiToEE, fValues);
         if(!IsPairPreFilterSelected(fValues)) {
            track->ResetFlags(); 
            break;
         }
      }  // end loop over positive prefilter tracks
      
      nextNegPrefilterTrack.Reset();
      for(Int_t ipn = 0; ipn<fPrefilterNegTracks.GetEntries(); ++ipn) {
         trackPref = (AliReducedBaseTrack*)nextNegPrefilterTrack();
         
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
     track = (AliReducedBaseTrack*)nextPosTrack();
     if(!track->GetFlags()) {
        fPosTracks.Remove(track);
     }
   }
  nextNegTrack.Reset();
  for(Int_t ip = fNegTracks.GetEntries()-1 ; ip >= 0; --ip) {
    track = (AliReducedBaseTrack*)nextNegTrack();
    if(!track->GetFlags()) {
       fNegTracks.Remove(track);
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::Finish() {
  //
  // run stuff after the event loop
  //
   if(fOptionRunMixing && !fOptionRunOverMC)
     fMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
}


//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2ee::CheckReconstructedLegMCTruth(AliReducedBaseTrack* track) {
   //
   // Check a reconstructed track against all the specified MC truth cuts
   //
   // TODO: In the fLegCandidatesMCcuts one can also add AliSignalMC objects which can then be tested
   //             using the AliReducedTrackInfo::fMCPdg[]
   //
   if(fLegCandidatesMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fLegCandidatesMCcuts.At(i);
      if(cut->IsSelected(track)) 
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2ee::CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack) {
   //
   // check the pair of tracks to see if they match the defined MC cuts and in addition
   // that they have the same mother
   // NOTE: The condition for the 2 tracks to have the same mother requires information on the MC label,
   //             which is available just in the full track information (AliReducedTrackInfo::fMCLabels[]).
   //           The consequence is that for the jpsi2ee analysis, the reconstructed tracks need to be always written as full tracks 
   //
   
   // check that both tracks are full tracks
   if(ptrack->IsA() != AliReducedTrackInfo::Class()) return 0;
   if(ntrack->IsA() != AliReducedTrackInfo::Class()) return 0;
   
   // check the MC requirements on each of the leg and their logical intersection
   if(fLegCandidatesMCcuts.GetEntries()==0) return 0;
   UInt_t pTrackDecisions = CheckReconstructedLegMCTruth(ptrack);
   if(!pTrackDecisions) return 0;
   UInt_t nTrackDecisions = CheckReconstructedLegMCTruth(ntrack);
   
   UInt_t decisions = 0;
   for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
      Bool_t pDecision = ( pTrackDecisions & (UInt_t(1)<<i) );
      Bool_t nDecision = ( nTrackDecisions & (UInt_t(1)<<i) );
      
      // if needed, check that the tracks have the same mother
      Bool_t sameMotherDecision = kTRUE;
      if(fLegCandidatesMCcuts_RequestSameMother[i])
         sameMotherDecision = (TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) == TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1)));
      if(sameMotherDecision && pDecision && nDecision)
         decisions |= (UInt_t(1)<<i);
   }
   
   return decisions;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FillMCTruthHistograms() {
  //
  // fill histograms with pure signal
  //   
   // loop over the first track array
  LoopOverMCTracks(1);
  // and over the second
  // NOTE: In the current model, handling the MC truth info requires the labels, which are properties of the full track,
  //         so there is no point in looping over the second track array which, if it exists, contains just base tracks
  //LoopOverMCTracks(2);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::LoopOverMCTracks(Int_t trackArray /*=1*/) {
   //
   // loop over the track array and check the pure MC tracks against the defined MC selections
   //   
   AliReducedTrackInfo* mother=0x0;
   AliReducedTrackInfo* daughter1 = 0x0;
   AliReducedTrackInfo* daughter2 = 0x0;
   
   TClonesArray* trackList = (trackArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if(!trackList) return;
   TIter nextTrack(trackList);
   
   // if the pt dependent weights were set, check the weight and reject randomly the event
   if(fMCJpsiPtWeights) {
      for(Int_t it=0; it<trackList->GetEntries(); ++it) {
         mother = (AliReducedTrackInfo*)nextTrack();
         if(!mother->IsMCKineParticle()) continue;
      
         // apply selections on the jpsi mother
         UInt_t motherDecisions = CheckMotherMCTruth(mother,kTRUE);
         if(!motherDecisions) continue;
         
         Double_t pt = mother->Pt();
         if(pt>fMCJpsiPtWeights->GetXaxis()->GetXmax()) 
            pt = fMCJpsiPtWeights->GetXaxis()->GetXmax();
         Double_t weight = fMCJpsiPtWeights->GetBinContent(fMCJpsiPtWeights->FindBin(pt));
         if(weight>1.0) weight = 1.0;
         Double_t rnd = gRandom->Rndm(); 
         if(weight<rnd) {
            fSkipMCEvent = kTRUE;
            return;
         }
      }
   }
   
   
   // Loop through the MC signals and check whether there are required MC signals without a defined mother 
   UInt_t signalsWithoutMother = 0;
   for(Int_t i=0; i<fJpsiMotherMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiMotherMCcuts.At(i);
      if(!((AliReducedTrackCut*)cut)->GetMCFilterMap()) 
         signalsWithoutMother |= (UInt_t(1)<<i);
   }
   
   // Loop through the list of pure MC particles and find the required mothers and their corresponding daughters
   nextTrack.Reset();
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      mother = (AliReducedTrackInfo*)nextTrack();
      if(!mother->IsMCKineParticle()) continue;
      
      // apply selections on the jpsi mother
      UInt_t motherDecisions = CheckMotherMCTruth(mother);
      if(!motherDecisions) continue;
      
      // find the jpsi daughters (needed to compute 2-track properties like the polarization, etc.)
      Int_t daughter1Label = 0; Int_t daughter2Label = 0;
      FindJpsiTruthLegs(mother, daughter1Label, daughter2Label);   
      daughter1 = FindMCtruthTrackByLabel(daughter1Label);
      daughter2 = FindMCtruthTrackByLabel(daughter2Label);
      
      // reset track variables and fill info
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
      AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);
      
      // loop over jpsi mother selections and fill histograms before the kine cuts on electrons
      for(Int_t iCut = 0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
         if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
         fHistosManager->FillHistClass(Form("%s_PureMCTruth_BeforeSelection", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);         
      }
      
      if(!daughter1) continue;
      if(!daughter2) continue;
      
      // apply selections on pure MC daughter electrons (kine cuts)
      UInt_t daughter1Decisions = CheckDaughterMCTruth(daughter1);
      if(!daughter1Decisions) continue;
      UInt_t daughtersDecisions = daughter1Decisions & CheckDaughterMCTruth(daughter2);
      if(!daughtersDecisions) continue;
      
      for(Int_t iCut = 0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
         if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
         if(!(daughtersDecisions & (UInt_t(1)<<iCut)))  continue;
         fHistosManager->FillHistClass(Form("%s_PureMCTruth_AfterSelection", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);         
      }
   }  // end loop over tracks
   
   
   // Loop through the list of pure MC particles and find those signals which do not require a mother (e.g. pairs of electrons without a mother)
   // NOTE: At the moment, just signals with exactly 2 prongs per event are allowed
   if(signalsWithoutMother) {
      nextTrack.Reset();
      daughter1 = 0x0;
      daughter2 = 0x0;
      UInt_t daughter1Decisions = 0;
      UInt_t daughter2Decisions = 0;
      AliReducedTrackInfo* tempTrack = 0x0;
      for(Int_t it=0; it<trackList->GetEntries(); ++it) {
         tempTrack = (AliReducedTrackInfo*)nextTrack();
         if(!tempTrack->IsMCKineParticle()) continue;
         
         UInt_t daughterDecisions = CheckDaughterMCTruth(tempTrack);
         daughterDecisions &= signalsWithoutMother;
         if(!daughterDecisions) continue;
         
         if(!daughter1) {
            daughter1 = tempTrack;
            daughter1Decisions = daughterDecisions;
            continue;
         }
         if(daughter1 && !daughter2) {
            daughter2 = tempTrack;
            daughter2Decisions = daughterDecisions;
         }
         
         if(daughter1 && daughter2) {
            daughterDecisions = daughter1Decisions & daughter2Decisions;
            for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
            AliReducedVarManager::FillMCTruthInfo(daughter1, daughter2, fValues);
            
            // loop over cuts and fill histograms
            for(Int_t iCut = 0; iCut<fJpsiElectronMCcuts.GetEntries(); ++iCut) {
               if(!(daughterDecisions & (UInt_t(1)<<iCut)))  continue;
               // the same information is filled in both Before and After histogram lists
               fHistosManager->FillHistClass(Form("%s_PureMCTruth_BeforeSelection", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);         
               fHistosManager->FillHistClass(Form("%s_PureMCTruth_AfterSelection", fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);         
            }
         }
      }   // end loop over tracks
   }   // end if(signalsWithoutMother)
   
   return;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2ee::CheckMotherMCTruth(AliReducedTrackInfo* mother, Bool_t checkReweight) {
   //
   // Check the mother pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fJpsiMotherMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fJpsiMotherMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiMotherMCcuts.At(i);
      // If no MC bit was requested for the mother, skip this mother signal
      // The MC cut will be applied just at the daughter level (this is likely an MC signal without a mother, e.g. electrons from gamma-gamma continuum from Starlight)
       
      //check if reweight is needed for this MC signal
       if (checkReweight && (((AliReducedTrackCut*)cut)->GetApplyReweightMCpt())==kFALSE) continue;
      if(!((AliReducedTrackCut*)cut)->GetMCFilterMap()) continue;
      if(cut->IsSelected(mother)) 
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
UInt_t AliReducedAnalysisJpsi2ee::CheckDaughterMCTruth(AliReducedTrackInfo* daughter) {
   //
   // Check the daughter pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fJpsiElectronMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fJpsiElectronMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiElectronMCcuts.At(i);
      if(cut->IsSelected(daughter))
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
AliReducedTrackInfo* AliReducedAnalysisJpsi2ee::FindMCtruthTrackByLabel(Int_t label) {
   //
   // search the track list for pure MC track with label and return the track pointer
   //
   AliReducedTrackInfo* track=0x0;
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCKineParticle()) continue;
      if(track->MCLabel(0)==label) return track;
   }
   return 0x0;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2ee::FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label) {
   //
   // find the jpsi legs in the list of pure MC truth particles
   //
   Int_t mLabel = mother->MCLabel(0);
   AliReducedTrackInfo* track=0x0;
   
   Int_t legsFound = 0;
   
   // loop over the first track array
   TClonesArray* trackList = fEvent->GetTracks();
   TIter nextTrack(trackList);
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      if(legsFound==2) return;
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCKineParticle()) continue;
      if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==11) {
         legsFound += 1;
         if(legsFound==1) leg1Label = track->MCLabel(0);
         if(legsFound==2) leg2Label = track->MCLabel(0);
      }
   }
   return;
}
