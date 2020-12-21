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
  fCorrelationsMixingHandler(new AliMixingHandler("J/psi - hadron correlations","",AliMixingHandler::kMixCorrelation)),
  fOptionUseLikeSignPairs(fOptionRunLikeSignPairing),
  fOptionRunCorrelation(kTRUE),
  fOptionRunCorrelationMixing(kTRUE),
  fAssociatedTrackCuts(),
  fAssociatedTracks(),
  fAssociatedTracksMB()
{
  //
  // default constructor
  //
}

//___________________________________________________________________________
AliReducedAnalysisJpsi2eeCorrelations::AliReducedAnalysisJpsi2eeCorrelations(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisJpsi2ee(name, title),
  fCorrelationsMixingHandler(new AliMixingHandler("J/psi - hadron correlations", "", AliMixingHandler::kMixCorrelation)),
  fOptionUseLikeSignPairs(fOptionRunLikeSignPairing),
  fOptionRunCorrelation(kTRUE),
  fOptionRunCorrelationMixing(kTRUE),
  fMBEventCuts(),
  fAssociatedTrackCuts(),
  fAssociatedTracks(),
  fAssociatedTracksMB()
{
  //
  // named constructor
  //
  fMBEventCuts.SetOwner(kTRUE);
  fAssociatedTrackCuts.SetOwner(kTRUE);
  fAssociatedTracks.SetOwner(kFALSE);
  fAssociatedTracksMB.SetOwner(kFALSE);
}

//___________________________________________________________________________
AliReducedAnalysisJpsi2eeCorrelations::~AliReducedAnalysisJpsi2eeCorrelations()
{
  //
  // destructor
  //
  fMBEventCuts.Clear("C");
  fAssociatedTrackCuts.Clear("C");
  fAssociatedTracks.Clear("C");
  fAssociatedTracksMB.Clear("C");
  if(fCorrelationsMixingHandler) delete fCorrelationsMixingHandler;
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::Init() {
  //
  // initializer
  //
  AliReducedAnalysisJpsi2ee::Init();
  fCorrelationsMixingHandler->SetHistogramManager(fHistosManager);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::Process() {
   //
   // process current event
   //
   // The AliReducedAnalysisJpsi2ee ancestor class must be set via the options to run the same event pairing
   // The jpsi pair candidates will then be stored in the array fJpsiCandidates
   AliReducedAnalysisJpsi2ee::Process();

  // check event cuts
  Bool_t eventSelected = IsEventSelected(fEvent, fValues);

  // do mixing between MB and triggered event
  if (fMBEventCuts.GetEntries() && !fOptionRunOverMC && fOptionRunCorrelationMixing && fOptionRunCorrelation) {

    // take hadrons from MB event
    if (IsMBEventSelected(fEvent, fValues)) {
      RunAssociatedTrackSelection(kFALSE, kTRUE);
      fCorrelationsMixingHandler->FillEvent(NULL, &fAssociatedTracksMB, fValues);
    }

    // take Jpsi canidates from triggered event
    if (eventSelected) {
      if (!fOptionUseLikeSignPairs) {
        TList jpsiCandidatesTemp;
        TIter nextJpsi(&fJpsiCandidates);
        AliReducedPairInfo* jpsi = 0x0;
        for (Int_t i=0; i<fJpsiCandidates.GetEntries(); ++i) {
          jpsi = (AliReducedPairInfo*)nextJpsi();
          if ((Int_t)jpsi->PairType()==1) jpsiCandidatesTemp.Add(jpsi->Clone());
        }
        fCorrelationsMixingHandler->FillEvent(&jpsiCandidatesTemp, NULL, fValues);
      } else {
        fCorrelationsMixingHandler->FillEvent(&fJpsiCandidates, NULL, fValues);
      }
    }
  }

  // apply event selection
  if (!eventSelected) return;

  // loop over associated tracks
  RunAssociatedTrackSelection();
  FillAssociatedTrackHistograms();

  // Feed the selected triggers and associated to the event mixing handler
  if(!fMBEventCuts.GetEntries() && !fOptionRunOverMC && fOptionRunCorrelationMixing && fOptionRunCorrelation) {
    if (!fOptionUseLikeSignPairs) {
      TList jpsiCandidatesTemp;
      TIter nextJpsi(&fJpsiCandidates);
      AliReducedPairInfo* jpsi = 0x0;
      for (Int_t i=0; i<fJpsiCandidates.GetEntries(); ++i) {
        jpsi = (AliReducedPairInfo*)nextJpsi();
        if ((Int_t)jpsi->PairType()==1) jpsiCandidatesTemp.Add(jpsi->Clone());
      }
      fCorrelationsMixingHandler->FillEvent(&jpsiCandidatesTemp, &fAssociatedTracks, fValues);
    } else {
      fCorrelationsMixingHandler->FillEvent(&fJpsiCandidates, &fAssociatedTracks, fValues);
    }
  }

   // run same event correlation
   if(fOptionRunCorrelation) RunSameEventCorrelation();
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::Finish() {
  //
  // after event loop (left over mixing at some point...)
  //
   AliReducedAnalysisJpsi2ee::Finish();
   if(fOptionRunCorrelation && fOptionRunCorrelationMixing && !fOptionRunOverMC)
     fCorrelationsMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::AddAssociatedTrackCut(AliReducedInfoCut* cut) {
  //
  // add a new associated track cut
  //
  fAssociatedTrackCuts.Add(cut);

  // NOTE: number of parallel must be equal to number of hist classes in AliMixingHandler::Init()
  if (fOptionUseLikeSignPairs)  fCorrelationsMixingHandler->SetNParallelCuts(fCorrelationsMixingHandler->GetNParallelCuts()+3);
  else                          fCorrelationsMixingHandler->SetNParallelCuts(fCorrelationsMixingHandler->GetNParallelCuts()+1);
  if (fPairCuts.GetEntries())   fCorrelationsMixingHandler->SetNParallelPairCuts(fPairCuts.GetEntries());
  if (fPairCuts.GetEntries()>1) {
    TString histClassNamesNew  = "";
    for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
      for (Int_t iTrackCut=0; iTrackCut<fAssociatedTrackCuts.GetEntries(); iTrackCut++) {
        if (fOptionUseLikeSignPairs) histClassNamesNew += Form("CorrMEPP_%s_%s_%s;",
                                                               fTrackCuts.At(iTrackCut)->GetName(),
                                                               fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                                               fPairCuts.At(iPairCut)->GetName());
        histClassNamesNew += Form("CorrMEPM_%s_%s_%s;",
                                  fTrackCuts.At(iTrackCut)->GetName(),
                                  fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                  fPairCuts.At(iPairCut)->GetName());
        if (fOptionUseLikeSignPairs) histClassNamesNew += Form("CorrMEMM_%s_%s_%s;",
                                                               fTrackCuts.At(iTrackCut)->GetName(),
                                                               fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                                               fPairCuts.At(iPairCut)->GetName());
      }
    }
    fCorrelationsMixingHandler->SetHistClassNames(histClassNamesNew.Data());
  } else {
    TString histClassNames = fCorrelationsMixingHandler->GetHistClassNames();
    if (fOptionUseLikeSignPairs) histClassNames += Form("CorrMEPP_%s_%s;", fTrackCuts.At(fAssociatedTrackCuts.GetEntries()-1)->GetName(), cut->GetName());
    histClassNames += Form("CorrMEPM_%s_%s;", fTrackCuts.At(fAssociatedTrackCuts.GetEntries()-1)->GetName(), cut->GetName());
    if (fOptionUseLikeSignPairs) histClassNames += Form("CorrMEMM_%s_%s;", fTrackCuts.At(fAssociatedTrackCuts.GetEntries()-1)->GetName(), cut->GetName());
    fCorrelationsMixingHandler->SetHistClassNames(histClassNames.Data());
  }
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeCorrelations::IsMBEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //
  if(fMBEventCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fMBEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fMBEventCuts.At(i);
    if(values) { if(!cut->IsSelected(event, values)) return kFALSE; }
    else { if(!cut->IsSelected(event)) return kFALSE; }
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisJpsi2eeCorrelations::IsAssociatedTrackSelected(AliReducedBaseTrack* track, Float_t* values /*=0x0*/) {
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
void AliReducedAnalysisJpsi2eeCorrelations::RunAssociatedTrackSelection(Bool_t fillHistograms /*=kTRUE*/, Bool_t fillMBTracks /*=kFALSE*/) {
  //
  // select associated tracks
  //
  // clear the track arrays
  if (!fillMBTracks)  fAssociatedTracks.Clear("C");
  else                fAssociatedTracksMB.Clear("C");
  // loop over the track list and evaluate all the track cuts
  LoopOverAssociatedTracks(1, fillHistograms, fillMBTracks);
  LoopOverAssociatedTracks(2, fillHistograms, fillMBTracks);
  // set number of associated tracks
  fValues[AliReducedVarManager::kNtracksAnalyzed] = fAssociatedTracks.GetEntries();
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::LoopOverAssociatedTracks(Int_t arrayOption /*=1*/, Bool_t fillHistograms /*=kTRUE*/, Bool_t fillMBTracks /*=kFALSE*/) {
   //
   // loop over the given track array, select tracks and fill histograms
   //
   AliReducedTrackInfo*  track     = 0x0;
   TClonesArray*         trackList = (arrayOption==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if (!trackList) return;
   TIter nextTrack(trackList);
   for (Int_t itr=0; itr<trackList->GetEntries(); ++itr) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(fOptionRunOverMC && track->IsMCTruth()) continue;
      AliReducedVarManager::FillTrackInfo(track, fValues);
      AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
      if (fillHistograms) fHistosManager->FillHistClass("AssociatedTrack_BeforeCuts", fValues);
      for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingStatus; ++iflag) {
         AliReducedVarManager::FillTrackingFlag(track, iflag, fValues);
         if (fillHistograms) fHistosManager->FillHistClass("AssociatedTrackStatusFlags_BeforeCuts", fValues);
      }
      for(Int_t iLayer=0; iLayer<6; ++iLayer) {
         AliReducedVarManager::FillITSlayerFlag(track, iLayer, fValues);
         if (fillHistograms) fHistosManager->FillHistClass("AssociatedTrackITSclusterMap_BeforeCuts", fValues);
      }
      for(Int_t iLayer=0; iLayer<8; ++iLayer) {
         AliReducedVarManager::FillTPCclusterBitFlag(track, iLayer, fValues);
         if (fillHistograms) fHistosManager->FillHistClass("AssociatedTrackTPCclusterMap_BeforeCuts", fValues);
      }
      if(IsAssociatedTrackSelected(track, fValues)) {
         if (!fillMBTracks) fAssociatedTracks.Add(track);
         else               fAssociatedTracksMB.Add(track);
      }
   }  // end loop over tracks
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::RunSameEventCorrelation(TString pairClass /*="PairSE"*/) {
  //
  // run the same event pairing for candidates (e+e-) and the correlation to associated tracks
  //
  if(fJpsiCandidates.GetEntries()==0) return;
  if(fAssociatedTracks.GetEntries()==0) return;
  TString pairTypeNames[3] = {"PP","PM","MM"};

  TIter nextAssocTrack(&fAssociatedTracks);
  TIter nextJpsi(&fJpsiCandidates);
  
  AliReducedPairInfo* jpsi = 0x0;
  AliReducedBaseTrack* assoc = 0x0;
  for(Int_t it=0;it<fJpsiCandidates.GetEntries(); ++it) {
     jpsi = (AliReducedPairInfo*)nextJpsi();
     
     nextAssocTrack.Reset();
     for(Int_t ia=0;ia<fAssociatedTracks.GetEntries(); ++ia) {
        assoc = (AliReducedBaseTrack*)nextAssocTrack();
        
        // make sure we do not correlate with one of the jpsi legs
        if(assoc->TrackId()==jpsi->LegId(0)) continue;
        if(assoc->TrackId()==jpsi->LegId(1)) continue;
        // NOTE: the model is that there is a set of n-selections for the jpsi and n-selections for the assoc
        //       One needs to have at least one matching bit in order to correlate them
        // TODO: we need to make sure there are equal numbers of bits for both trigger and assoc
        //           Right now, the extra bits from the particle with more bits (cuts defined) are ignored
        if(!(jpsi->GetFlags() & assoc->GetFlags())) continue;
        if (!fOptionUseLikeSignPairs && ((Int_t)jpsi->PairType() == 0 || (Int_t)jpsi->PairType() == 2)) continue;

        AliReducedVarManager::FillCorrelationInfo(jpsi, assoc, fValues);

        // fill correlation histograms
        // TODO: isMCTruth must be handled; can be either a Bool or a bit map
        //             Not sure if we need MC truth for correlation -> if so remove from code
        Bool_t isMCTruth = kFALSE;
        FillCorrelationHistograms(jpsi, assoc, Form("CorrSE%s", pairTypeNames[(Int_t)jpsi->PairType()].Data()), isMCTruth);
     }  // end loop over associated tracks
  }  // end loop over jpsi candidates
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::FillAssociatedTrackHistograms(TString trackClass /*="AssociatedTrack"*/) {
  //
  // fill associated track histograms
  //
  AliReducedTrackInfo* track=0;
  TIter nextAssocTrack(&fAssociatedTracks);
  for(Int_t i=0;i<fAssociatedTracks.GetEntries();++i) {
    track = (AliReducedTrackInfo*)nextAssocTrack();
    AliReducedVarManager::FillTrackInfo(track, fValues);
    AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
    FillAssociatedTrackHistograms(track, trackClass);
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::FillAssociatedTrackHistograms(AliReducedTrackInfo* track, TString trackClass  /*="AssociatedTrack"*/) {
  //
  // fill associated track histograms
  //
  //Bool_t isMCTruth = fOptionRunOverMC && track->IsMCTruth();
  Bool_t isMCTruth = kFALSE;     // TODO: handle the MC info if needed
  for(Int_t icut=0; icut<fAssociatedTrackCuts.GetEntries(); ++icut) {
    if(track->TestFlag(icut)) {
      fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
      //if(isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
      for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
        AliReducedVarManager::FillTrackingFlag(track, iflag, fValues);
        fHistosManager->FillHistClass(Form("%sStatusFlags_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
        if(isMCTruth) fHistosManager->FillHistClass(Form("%sStatusFlags_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
      }
      for(Int_t iLayer=0; iLayer<6; ++iLayer) {
        AliReducedVarManager::FillITSlayerFlag(track, iLayer, fValues);
        fHistosManager->FillHistClass(Form("%sITSclusterMap_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
        if(isMCTruth) fHistosManager->FillHistClass(Form("%sITSclusterMap_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
      }
      for(Int_t iLayer=0; iLayer<8; ++iLayer) {
        AliReducedVarManager::FillTPCclusterBitFlag(track, iLayer, fValues);
        fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
        if(isMCTruth) fHistosManager->FillHistClass(Form("%sTPCclusterMap_%s_MCTruth", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
      }
    }
  }
}

//___________________________________________________________________________
void AliReducedAnalysisJpsi2eeCorrelations::FillCorrelationHistograms(AliReducedPairInfo* jpsi, AliReducedBaseTrack* assoc, TString corrClass /*="CorrSE"*/, Bool_t isMCTruth/*=kFALSE*/) {
  //
  // fill correlation histograms
  //
  ULong_t trackMask = jpsi->GetFlags() & assoc->GetFlags();
  ULong_t pairMask  = jpsi->GetQualityFlags();
  if (fPairCuts.GetEntries()>1) {
    for(Int_t iTrackCut=0; iTrackCut<fAssociatedTrackCuts.GetEntries(); ++iTrackCut) {
      for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
        if((trackMask & (ULong_t(1)<<iTrackCut)) && (pairMask & (ULong_t(1)<<iPairCut))) {
          fHistosManager->FillHistClass(Form("%s_%s_%s_%s", corrClass.Data(),
                                             fTrackCuts.At(iTrackCut)->GetName(),
                                             fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                             fPairCuts.At(iPairCut)->GetName()), fValues);
          if(isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_%s_MCTruth", corrClass.Data(),
                                                           fTrackCuts.At(iTrackCut)->GetName(),
                                                           fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                                           fPairCuts.At(iPairCut)->GetName()), fValues);
        }
      }
    }
  } else {
    for(Int_t iCut=0; iCut<fAssociatedTrackCuts.GetEntries(); ++iCut) {
       if(!(trackMask & (ULong_t(1)<<iCut))) continue;
       fHistosManager->FillHistClass(Form("%s_%s_%s", corrClass.Data(), fTrackCuts.At(iCut)->GetName(), fAssociatedTrackCuts.At(iCut)->GetName()), fValues);
       if(isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_%s_MCTruth", corrClass.Data(), fTrackCuts.At(iCut)->GetName(), fAssociatedTrackCuts.At(iCut)->GetName()), fValues);
    }
  }
}
