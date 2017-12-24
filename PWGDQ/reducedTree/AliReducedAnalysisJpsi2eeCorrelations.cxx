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
// Analysis task for J/psi - hadron correlations                         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliReducedAnalysisJpsi2eeCorrelations.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>
#include <TList.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisJpsi2eeCorrelations);

//___________________________________________________________________________
AliReducedAnalysisJpsi2eeCorrelations::AliReducedAnalysisJpsi2eeCorrelations() :
  AliReducedAnalysisJpsi2ee(),
  fElectronBit(32),
  fOptionRunCorrelation(kTRUE),
  fOptionAssociatedTracks(kFALSE),
  fAssociatedTrackCuts(),
  fAssociatedTracks()
{
  //
  // default constructor
  //
}

//___________________________________________________________________________
AliReducedAnalysisJpsi2eeCorrelations::AliReducedAnalysisJpsi2eeCorrelations(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisJpsi2ee(name, title),
  fElectronBit(32),
  fOptionRunCorrelation(kTRUE),
  fOptionAssociatedTracks(kFALSE),
  fAssociatedTrackCuts(),
  fAssociatedTracks()
{
  //
  // named constructor
  //
  fAssociatedTrackCuts.SetOwner(kTRUE);
  fAssociatedTracks.SetOwner(kFALSE);
}

//___________________________________________________________________________
AliReducedAnalysisJpsi2eeCorrelations::~AliReducedAnalysisJpsi2eeCorrelations()
{
  //
  // destructor
  //
  fAssociatedTrackCuts.Clear("C");
  fAssociatedTracks.Clear("C");
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::Init() {
  //
  // initializer
  //
  AliReducedVarManager::SetDefaultVarNames();
  fHistosManager->SetUseDefaultVariableNames(kTRUE);
  fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::Process() {
  //
  // process current event
  //
  if(!fEvent) return;
  AliReducedEventInfo* eventInfo = NULL;
  if (fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
  else {
    cout << "ERROR: AliReducedAnalysisJpsi2eeCorrelations::Process() needs AliReducedEventInfo events" << endl;
    return;
  }

  if(fOptionRunOverMC) {
    if(fEventCounter%10000==0) cout << "Event no. " << fEventCounter << endl;
  } else {
    if(fEventCounter%100000==0) cout << "Event no. " << fEventCounter << endl;
  }
  fEventCounter++;

  AliReducedVarManager::SetEvent(fEvent);

  // reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: the run wise data will be updated automatically in the VarManager in case a run change is detected
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i) {
    fValues[i]  = -9999.;
    fValues2[i] = -9999.;
  }

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

  // fill MC truth histograms
  if(fOptionRunOverMC) FillMCTruthHistograms();

  // loop over tracks
  if(fOptionLoopOverTracks) {
    // track and prefilter selection
    RunTrackSelection();
    RunPrefilter();

    // set values
    fValues[AliReducedVarManager::kNtracksPosAnalyzed] = fPosTracks.GetEntries();
    fValues[AliReducedVarManager::kNtracksNegAnalyzed] = fNegTracks.GetEntries();
    fValues[AliReducedVarManager::kNtracksAnalyzed] = fValues[AliReducedVarManager::kNtracksNegAnalyzed]+fValues[AliReducedVarManager::kNtracksPosAnalyzed];
    fValues[AliReducedVarManager::kEvAverageTPCchi2] /= (fPosTracks.GetEntries()+fNegTracks.GetEntries()>0 ? fValues[AliReducedVarManager::kNtracksAnalyzed] : 1.0);

    // fill track histograms
    FillTrackHistograms();
  }

  // loop over associated tracks
  if (fOptionRunCorrelation) {
    // track selection
    RunAssociatedTrackSelection();

    // fill associated track histograms
    FillAssociatedTrackHistograms();
  }

  // run same event correlation
  if (fOptionRunCorrelation) RunSameEventCorrelation();
  
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
void AliReducedAnalysisJpsi2eeCorrelations::Finish() {
  //
  // after event loop (left over mixing at some point...)
  //
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::AddAssociatedTrackCut(AliReducedInfoCut* cut) {
  //
  // add a new associated track cut
  // NOTE: this probably needs some updates for event mixing at some point (see AliReducedAnalysisJpsi2ee::AddTrackCut)
  //
  fAssociatedTrackCuts.Add(cut);
}
  
//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeCorrelations::IsAssociatedTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
  //
  // apply associated track cuts
  //
  if(fAssociatedTrackCuts.GetEntries()==0) return kTRUE;
  track->ResetFlags();
  for(Int_t i=0; i<fAssociatedTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fAssociatedTrackCuts.At(i);
    if(values){ if(cut->IsSelected(track, values)) track->SetFlag(i); }
    else      { if(cut->IsSelected(track)) track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::RunAssociatedTrackSelection() {
  //
  // select associated tracks
  //

  // clear the track arrays
  fAssociatedTracks.Clear("C");
  fValues2[AliReducedVarManager::kEvAverageTPCchi2] = 0.0;

  // loop over the track list and evaluate all the track cuts
  AliReducedTrackInfo*  track     = 0x0;
  TClonesArray*         trackList = 0x0;
  if (fOptionAssociatedTracks)  trackList = fEvent->GetTracks2();
  else                          trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for (Int_t itr=0; itr<trackList->GetEntries(); ++itr) {
    track = (AliReducedTrackInfo*)nextTrack();
    if(fOptionRunOverMC && track->IsMCTruth()) continue;
    AliReducedVarManager::FillTrackInfo(track, fValues2);
    fHistosManager->FillHistClass("AssociatedTrack_BeforeCuts", fValues2);
    for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
      AliReducedVarManager::FillTrackingFlag(track, iflag, fValues2);
      fHistosManager->FillHistClass("AssociatedTrackStatusFlags_BeforeCuts", fValues2);
    }
    for(Int_t iLayer=0; iLayer<6; ++iLayer) {
      AliReducedVarManager::FillITSlayerFlag(track, iLayer, fValues2);
      fHistosManager->FillHistClass("AssociatedTrackITSclusterMap_BeforeCuts", fValues2);
    }
    for(Int_t iLayer=0; iLayer<8; ++iLayer) {
      AliReducedVarManager::FillTPCclusterBitFlag(track, iLayer, fValues2);
      fHistosManager->FillHistClass("AssociatedTrackTPCclusterMap_BeforeCuts", fValues2);
    }
    if(IsAssociatedTrackSelected(track, fValues2)) {
      fValues2[AliReducedVarManager::kEvAverageTPCchi2] += track->TPCchi2();
      fAssociatedTracks.Add(track);
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::RunSameEventCorrelation(TString pairClass/*="PairSE"*/) {
  //
  // run the same event pairing for candidates (e+e-) and the correlation to associated tracks
  //
  fValues[AliReducedVarManager::kNpairsSelected] = 0;

  TIter nextPosTrack(&fPosTracks);
  TIter nextNegTrack(&fNegTracks);
  TIter nextAssocTrack(&fAssociatedTracks);

  AliReducedTrackInfo* pTrack     = 0x0;
  AliReducedTrackInfo* nTrack     = 0x0;
  AliReducedTrackInfo* assocTrack = 0x0;
  for (Int_t ip=0; ip<fPosTracks.GetEntries(); ++ip) {
    pTrack = (AliReducedTrackInfo*)nextPosTrack();

    nextNegTrack.Reset();
    for (Int_t in=0; in<fNegTracks.GetEntries(); ++in) {
      nTrack = (AliReducedTrackInfo*)nextNegTrack();

      // verify that the two current tracks have at least 1 common bit
      if(!(pTrack->GetFlags() & nTrack->GetFlags())) continue;

      // fill pair info
      AliReducedVarManager::FillPairInfo(pTrack, nTrack, AliReducedPairInfo::kJpsiToEE, fValues);
      if(IsPairSelected(fValues)) {

        FillPairHistograms(pTrack->GetFlags() & nTrack->GetFlags(), 1, pairClass, fOptionRunOverMC && IsMCTruth(pTrack, nTrack)); // 1 is for +- pairs
        fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
        fValues[AliReducedVarManager::kTriggerPt] = fValues[AliReducedVarManager::kPt];

        // loop over associated tracks
        nextAssocTrack.Reset();
        for (Int_t ia=0; ia<fAssociatedTracks.GetEntries(); ++ia) {
          assocTrack = (AliReducedTrackInfo*)nextAssocTrack();

          // verify that the associated track is not a J/psi candidate leg
          if (assocTrack->GetQualityFlags()&(ULong_t(1)<<fElectronBit)) continue;

          // fill correlation histogram
          AliReducedVarManager::FillCorrelationInfo(assocTrack, fValues);
          FillCorrelationHistograms(pTrack->GetFlags()&nTrack->GetFlags(), assocTrack->GetFlags(), "CorrSE", fOptionRunOverMC && IsMCTruth(pTrack, nTrack) && assocTrack->IsMCTruth());
        }
      }
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::FillAssociatedTrackHistograms(TString trackClass/*="AssociatedTrack"*/) {
  //
  // fill associated track histograms
  //
  for(Int_t i=0;i<36; ++i) fValues2[AliReducedVarManager::kNtracksAnalyzedInPhiBins+i] = 0.;
  AliReducedTrackInfo* track=0;
  TIter nextAssocTrack(&fAssociatedTracks);
  for(Int_t i=0;i<fAssociatedTracks.GetEntries();++i) {
    track = (AliReducedTrackInfo*)nextAssocTrack();
    fValues2[AliReducedVarManager::kNtracksAnalyzedInPhiBins+(track->Eta()<0.0 ? 0 : 18) + TMath::FloorNint(18.*track->Phi()/TMath::TwoPi())] += 1;
    AliReducedVarManager::FillTrackInfo(track, fValues2);
    FillAssociatedTrackHistograms(track, trackClass);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::FillAssociatedTrackHistograms(AliReducedTrackInfo* track, TString trackClass/*="AssociatedTrack"*/) {
  //
  // fill associated track histograms
  //
  Bool_t isMCTruth = fOptionRunOverMC && track->IsMCTruth();
  for(Int_t icut=0; icut<fAssociatedTrackCuts.GetEntries(); ++icut) {
    if(track->TestFlag(icut)) {
      fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
      if(isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
      for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
        AliReducedVarManager::FillTrackingFlag(track, iflag, fValues2);
        fHistosManager->FillHistClass(Form("%sStatusFlags_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
        if(isMCTruth) fHistosManager->FillHistClass(Form("%sStatusFlags_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
      }
      for(Int_t iLayer=0; iLayer<6; ++iLayer) {
        AliReducedVarManager::FillITSlayerFlag(track, iLayer, fValues2);
        fHistosManager->FillHistClass(Form("%sITSclusterMap_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
        if(isMCTruth) fHistosManager->FillHistClass(Form("%sITSclusterMap_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
      }
      for(Int_t iLayer=0; iLayer<8; ++iLayer) {
        AliReducedVarManager::FillTPCclusterBitFlag(track, iLayer, fValues2);
        fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
        if(isMCTruth) fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues2);
      }
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::FillCorrelationHistograms(ULong_t maskTrack, ULong_t maskAssocTrack, TString corrClass/*="CorrSE"*/, Bool_t isMCTruth/*=kFALSE*/) {
  //
  // fill correlation histograms
  //
  for (Int_t iCutTrack=0; iCutTrack<fTrackCuts.GetEntries(); ++iCutTrack) {
    for (Int_t iCutAssocTrack=0; iCutAssocTrack<fAssociatedTrackCuts.GetEntries(); ++iCutAssocTrack) {
      if ( (maskTrack & (ULong_t(1)<<iCutTrack)) && (maskAssocTrack & (ULong_t(1)<<iCutAssocTrack)) ) {
        fHistosManager->FillHistClass(Form("%s_%s_%s", corrClass.Data(), fTrackCuts.At(iCutTrack)->GetName(),
                                           fAssociatedTrackCuts.At(iCutAssocTrack)->GetName()), fValues);
        if (isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_%s_MCTruth", corrClass.Data(), fTrackCuts.At(iCutTrack)->GetName(),
                                                          fAssociatedTrackCuts.At(iCutAssocTrack)->GetName()), fValues);
      }
    }
  }
}
