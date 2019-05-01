/**************************************************************************
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
// Analysis task for single tracks                                       //
//                                                                       //
// creation date: 10/10/2018                                             //
// author: Lucas Altenkamper, lucas.altenkamper@cern.ch                  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliReducedAnalysisSingleTrack.h"

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
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisSingleTrack);

//___________________________________________________________________________
AliReducedAnalysisSingleTrack::AliReducedAnalysisSingleTrack() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fClusterTrackMatcher(0x0),
  fOptionRunOverMC(kTRUE),
  fOptionRunOverCaloCluster(kFALSE),
  fEventCuts(),
  fTrackCuts(),
  fClusterCuts(),
  fMCSignalCuts(),
  fTracks(),
  fClusters(),
  fClusterTrackMatcherHistograms(0x0),
  fClusterTrackMatcherMultipleMatchesBefore(0x0),
  fClusterTrackMatcherMultipleMatchesAfter(0x0)
{
  //
  // default constructor
  //
}

//___________________________________________________________________________
AliReducedAnalysisSingleTrack::AliReducedAnalysisSingleTrack(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fClusterTrackMatcher(0x0),
  fOptionRunOverMC(kTRUE),
  fOptionRunOverCaloCluster(kFALSE),
  fEventCuts(),
  fTrackCuts(),
  fClusterCuts(),
  fMCSignalCuts(),
  fTracks(),
  fClusters(),
  fClusterTrackMatcherHistograms(0x0),
  fClusterTrackMatcherMultipleMatchesBefore(0x0),
  fClusterTrackMatcherMultipleMatchesAfter(0x0)
{
  //
  // named constructor
  //
  fEventCuts.SetOwner(kTRUE);
  fTrackCuts.SetOwner(kTRUE);
  fClusterCuts.SetOwner(kTRUE);
  fMCSignalCuts.SetOwner(kTRUE);
  fTracks.SetOwner(kFALSE);
  fClusters.SetOwner(kFALSE);
}

//___________________________________________________________________________
AliReducedAnalysisSingleTrack::~AliReducedAnalysisSingleTrack()
{
  //
  // destructor
  //
  fEventCuts.Clear("C");
  fTrackCuts.Clear("C");
  fClusterCuts.Clear("C");
  fMCSignalCuts.Clear("C");
  fTracks.Clear("C");
  fClusters.Clear("C");
  if (fHistosManager) delete fHistosManager;
  if (fClusterTrackMatcher) delete fClusterTrackMatcher;
  if (fClusterTrackMatcherHistograms) delete fClusterTrackMatcherHistograms;
  if (fClusterTrackMatcherMultipleMatchesBefore) delete fClusterTrackMatcherMultipleMatchesBefore;
  if (fClusterTrackMatcherMultipleMatchesAfter) delete fClusterTrackMatcherMultipleMatchesAfter;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisSingleTrack::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if (fEventCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for (Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if (values) { if (!cut->IsSelected(event, values)) return kFALSE; }
    else { if (!cut->IsSelected(event)) return kFALSE; }
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisSingleTrack::IsTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
  //
  // apply track cuts
  //
  if (fTrackCuts.GetEntries()==0) return kTRUE;
  track->ResetFlags();
  for (Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if (values) { if (cut->IsSelected(track, values)) track->SetFlag(i); }
    else { if (cut->IsSelected(track)) track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisSingleTrack::IsClusterSelected(AliReducedCaloClusterInfo* cluster, Float_t* values/*=0x0*/) {
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
UInt_t AliReducedAnalysisSingleTrack::CheckTrackMCTruth(AliReducedBaseTrack* track) {
  //
  // check a track against all the specified MC truth cuts
  //
  if (fMCSignalCuts.GetEntries()==0) return 0;
  UInt_t decisionMap = 0;
  for (Int_t i=0; i<fMCSignalCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fMCSignalCuts.At(i);
    if (cut->IsSelected(track))
      decisionMap |= (UInt_t(1)<<i);
  }
  return decisionMap;
}

//___________________________________________________________________________
void AliReducedAnalysisSingleTrack::FillMCTruthHistograms() {
  //
  // fill histograms with pure MC signal according to defined MC selections
  //
  TClonesArray* trackList = fEvent->GetTracks();
  if (!trackList) return;
  
  TIter nextTrack(trackList);
  AliReducedTrackInfo* track = 0x0;
  for (Int_t it=0; it<trackList->GetEntries(); ++it) {
    track = (AliReducedTrackInfo*)nextTrack();
    if (!track->IsMCKineParticle()) continue;
    
    // apply MC selections on the track
    UInt_t mcDecisionMap = CheckTrackMCTruth(track);
    if (!mcDecisionMap) continue;
    
    // reset track variables and fill info
    for (Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i] = -9999.;
    AliReducedVarManager::FillMCTruthInfo(track, fValues);
    
    // loop over track selections and fill histograms
    for (Int_t iCut = 0; iCut<fMCSignalCuts.GetEntries(); ++iCut) {
      if (!(mcDecisionMap & (UInt_t(1)<<iCut)))  continue;
      fHistosManager->FillHistClass(Form("%s_PureMCTruth", fMCSignalCuts.At(iCut)->GetName()), fValues);
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisSingleTrack::RunTrackSelection() {
  //
  // select tracks
  //
  fTracks.Clear("C");
  fValues[AliReducedVarManager::kEvAverageTPCchi2] = 0.0;
  
  TClonesArray* trackList = fEvent->GetTracks();
  if (!trackList) return;
  if (!trackList->GetEntries()) return;
  
  TIter nextTrack(trackList);
  AliReducedBaseTrack* track = 0x0;
  for (Int_t it=0; it<trackList->GetEntries(); ++it) {
    track = (AliReducedBaseTrack*)nextTrack();
    
    // do not loop over pure MC truth tracks
    // NOTE: taken from AliReducedAnalysisJpsi2ee::LoopOverTracks(), required here?
    if (fOptionRunOverMC && track->IsMCTruth()) continue;
    
    // reset track variables
    for (Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i] = -9999.;
    
    AliReducedVarManager::FillTrackInfo(track, fValues);
    if (fClusterCuts.GetEntries())  AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, &fClusters, fClusterTrackMatcher);
    else                            AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues, NULL, fClusterTrackMatcher);
    fHistosManager->FillHistClass("Track_BeforeCuts", fValues);
    
    if (track->IsA() == AliReducedTrackInfo::Class()) {
      AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
      if (trackInfo) {
        for (UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
          AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
          fHistosManager->FillHistClass("TrackStatusFlags_BeforeCuts", fValues);
        }
        for (Int_t iLayer=0; iLayer<6; ++iLayer) {
          AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
          fHistosManager->FillHistClass("TrackITSclusterMap_BeforeCuts", fValues);
          AliReducedVarManager::FillITSsharedLayerFlag(trackInfo, iLayer, fValues);
          fHistosManager->FillHistClass("TrackITSsharedClusterMap_BeforeCuts", fValues);
        }
        for (Int_t iLayer=0; iLayer<8; ++iLayer) {
          AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
          fHistosManager->FillHistClass("TrackTPCclusterMap_BeforeCuts", fValues);
        }
      }
    }
    
    if (IsTrackSelected(track, fValues)) {
      fTracks.Add(track);
      if (track->IsA()==AliReducedTrackInfo::Class()) fValues[AliReducedVarManager::kEvAverageTPCchi2] += ((AliReducedTrackInfo*)track)->TPCchi2();
    }
  } // end loop over tracks
}

//___________________________________________________________________________
void AliReducedAnalysisSingleTrack::FillTrackHistograms(TString trackClass/*="Track"*/) {
  //
  // fill track histograms
  //
  if (fClusterTrackMatcher) {
    fClusterTrackMatcher->ClearMatchedClusterIDsBefore();
    fClusterTrackMatcher->ClearMatchedClusterIDsAfter();
  }

  for (Int_t i=0;i<36; ++i) fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+i] = 0.;
  AliReducedBaseTrack* track = 0;
  TIter nextTrack(&fTracks);
  for (Int_t i=0;i<fTracks.GetEntries();++i) {
    track = (AliReducedBaseTrack*)nextTrack();
    fValues[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
    
    // reset track variables
    for (Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i] = -9999.;
    
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
void AliReducedAnalysisSingleTrack::FillTrackHistograms(AliReducedBaseTrack* track, TString trackClass/*="Track"*/) {
  //
  // fill track level histograms
  //
  UInt_t mcDecisionMap = 0;
  if (fOptionRunOverMC) mcDecisionMap = CheckTrackMCTruth(track);
  
  for (Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
    if (track->TestFlag(icut)) {
      fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
      if (mcDecisionMap) {
        for (Int_t iMC=0; iMC<=fMCSignalCuts.GetEntries(); ++iMC) {
          if (mcDecisionMap & (UInt_t(1)<<iMC))
            fHistosManager->FillHistClass(Form("%s_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                               fMCSignalCuts.At(iMC)->GetName()), fValues);
        }
      }
      
      if (track->IsA() != AliReducedTrackInfo::Class()) continue;
      
      AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
      if (!trackInfo) continue;
      
      for (UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
        AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
        fHistosManager->FillHistClass(Form("%sStatusFlags_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
        if (mcDecisionMap) {
          for (Int_t iMC=0; iMC<=fMCSignalCuts.GetEntries(); ++iMC) {
            if (mcDecisionMap & (UInt_t(1)<<iMC))
              fHistosManager->FillHistClass(Form("%sStatusFlags_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                 fMCSignalCuts.At(iMC)->GetName()), fValues);
          }
        }
      }
      for (Int_t iLayer=0; iLayer<6; ++iLayer) {
        AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
        fHistosManager->FillHistClass(Form("%sITSclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
        if (mcDecisionMap) {
          for (Int_t iMC=0; iMC<=fMCSignalCuts.GetEntries(); ++iMC) {
            if (mcDecisionMap & (UInt_t(1)<<iMC))
              fHistosManager->FillHistClass(Form("%sITSclusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                 fMCSignalCuts.At(iMC)->GetName()), fValues);
          }
        }
        AliReducedVarManager::FillITSsharedLayerFlag(trackInfo, iLayer, fValues);
        fHistosManager->FillHistClass(Form("%sITSsharedClusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
        if (mcDecisionMap) {
          for (Int_t iMC=0; iMC<=fMCSignalCuts.GetEntries(); ++iMC) {
            if (mcDecisionMap & (UInt_t(1)<<iMC))
              fHistosManager->FillHistClass(Form("%sITSsharedClusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                 fMCSignalCuts.At(iMC)->GetName()), fValues);
          }
        }
      }
      for (Int_t iLayer=0; iLayer<8; ++iLayer) {
        AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
        fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName()), fValues);
        if (mcDecisionMap) {
          for (Int_t iMC=0; iMC<=fMCSignalCuts.GetEntries(); ++iMC) {
            if (mcDecisionMap & (UInt_t(1)<<iMC))
              fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s_%s", trackClass.Data(), fTrackCuts.At(icut)->GetName(),
                                                 fMCSignalCuts.At(iMC)->GetName()), fValues);
          }
        }
      }
    } // end if (track->TestFlag(icut))
  } // end loop over cuts
}

//___________________________________________________________________________
void AliReducedAnalysisSingleTrack::RunClusterSelection() {
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
void AliReducedAnalysisSingleTrack::FillClusterHistograms(TString clusterClass/*="CaloCluster"*/) {
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
void AliReducedAnalysisSingleTrack::FillClusterHistograms(AliReducedCaloClusterInfo* cluster, TString clusterClass/*="CaloCluster"*/) {
  //
  // fill cluster histograms
  //
  for (Int_t icut=0; icut<fClusterCuts.GetEntries(); ++icut) {
    if (cluster->TestFlag(icut)) fHistosManager->FillHistClass(Form("%s_%s", clusterClass.Data(), fClusterCuts.At(icut)->GetName()), fValues);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisSingleTrack::Init() {
  //
  // initialize stuff
  //
  AliReducedVarManager::SetDefaultVarNames();
  fHistosManager->SetUseDefaultVariableNames(kTRUE);
  fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames, AliReducedVarManager::fgVariableUnits);

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
void AliReducedAnalysisSingleTrack::Process() {
  //
  // process the current event
  //
  if (!fEvent) return;
  AliReducedEventInfo* eventInfo = NULL;
  if (fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
  else {
    cout << "ERROR: AliReducedAnalysisSingleTrack::Process() needs AliReducedEventInfo events" << endl;
    return;
  }
  if (fOptionRunOverMC && (fEventCounter%10000==0)) cout << "Event no. " << fEventCounter << endl;
  else if (fEventCounter%100000==0)                 cout << "Event no. " << fEventCounter << endl;
  fEventCounter++;
  
  AliReducedVarManager::SetEvent(fEvent);
  
  // reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: the run wise data will be updated automatically in the VarManager in case a run change is detected
  for (Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i) fValues[i]=-9999.;

  // fill event information before event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_BeforeCuts", fValues);
  for (UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
    fHistosManager->FillHistClass("EventTag_BeforeCuts", fValues);
  }
  for (UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
    fHistosManager->FillHistClass("EventTriggers_BeforeCuts", fValues);
  }

  // apply event selection
  if (!IsEventSelected(fEvent)) return;
  
  // fill MC truth histograms
  if (fOptionRunOverMC) FillMCTruthHistograms();

  // select cluster and fill histograms
  if (fOptionRunOverCaloCluster) {
    RunClusterSelection();
    FillClusterHistograms();
  }

  // select tracks
  RunTrackSelection();
  fValues[AliReducedVarManager::kNtracksAnalyzed] = fTracks.GetEntries();
  fValues[AliReducedVarManager::kEvAverageTPCchi2] /= (fTracks.GetEntries()>0 ? fValues[AliReducedVarManager::kNtracksAnalyzed] : 1.0);

  // fill track histograms
  FillTrackHistograms();
  
  // fill event info histograms after cuts
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
  for (UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
    fHistosManager->FillHistClass("EventTag_AfterCuts", fValues);
  }
  for (UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
    fHistosManager->FillHistClass("EventTriggers_AfterCuts", fValues);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisSingleTrack::Finish() {
  //
  // run stuff after the event loop
  //
}
