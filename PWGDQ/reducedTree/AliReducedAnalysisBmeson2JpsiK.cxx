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

#include "AliReducedAnalysisBmeson2JpsiK.h"

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

ClassImp(AliReducedAnalysisBmeson2JpsiK);

//___________________________________________________________________________
AliReducedAnalysisBmeson2JpsiK::AliReducedAnalysisBmeson2JpsiK() :
  AliReducedAnalysisJpsi2ee(),
  fBmesonMixingHandler(new AliMixingHandler("Bmeson -> J/psi + K","",AliMixingHandler::kMixCorrelation)),
  fOptionUseLikeSignPairs(fOptionRunLikeSignPairing),
  fOptionRunTripletPairing(kTRUE),
  fOptionRunTripletPairingMixing(kFALSE),
  fMBEventCuts(),
  fAssociatedTrackCuts(),
  fAssociatedTracks(),
  fAssociatedTracksMB(),
  fAssociatedMCcuts(),
  fCommonAncestor(),
  fBcandidateMCcuts(),
  fBcandidateDaughter1MCcuts(),
  fBcandidateDaughter2MCcuts(),
  fTripletCut(0x0)
{
  //
  // default constructor
  //
}

//___________________________________________________________________________
AliReducedAnalysisBmeson2JpsiK::AliReducedAnalysisBmeson2JpsiK(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisJpsi2ee(name, title),
  fBmesonMixingHandler(new AliMixingHandler("Bmeson -> J/psi + K","",AliMixingHandler::kMixCorrelation)),
  fOptionUseLikeSignPairs(fOptionRunLikeSignPairing),
  fOptionRunTripletPairing(kTRUE),
  fOptionRunTripletPairingMixing(kFALSE),
  fMBEventCuts(),
  fAssociatedTrackCuts(),
  fAssociatedTracks(),
  fAssociatedTracksMB(),
  fAssociatedMCcuts(),
  fCommonAncestor(),
  fBcandidateMCcuts(),
  fBcandidateDaughter1MCcuts(),
  fBcandidateDaughter2MCcuts(),
  fTripletCut(0x0)
{
  //
  // named constructor
  //
  fMBEventCuts.SetOwner(kTRUE);
  fAssociatedTrackCuts.SetOwner(kTRUE);
  fAssociatedTracks.SetOwner(kFALSE);
  fAssociatedTracksMB.SetOwner(kFALSE);
  fAssociatedMCcuts.SetOwner(kFALSE);
  fBcandidateMCcuts.SetOwner(kTRUE);
  fBcandidateDaughter1MCcuts.SetOwner(kTRUE);
  fBcandidateDaughter2MCcuts.SetOwner(kTRUE);
}

//___________________________________________________________________________
AliReducedAnalysisBmeson2JpsiK::~AliReducedAnalysisBmeson2JpsiK()
{
  //
  // destructor
  //
  fMBEventCuts.Clear("C");
  fAssociatedTrackCuts.Clear("C");
  fAssociatedTracks.Clear("C");
  fAssociatedTracksMB.Clear("C");
  if(fBmesonMixingHandler) delete fBmesonMixingHandler;
  fAssociatedMCcuts.Clear("C");
  if(fTripletCut) 
    delete fTripletCut;
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::Init() {
  //
  // initializer
  //
  AliReducedAnalysisJpsi2ee::Init();
  fBmesonMixingHandler->SetHistogramManager(fHistosManager);
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::Process() {
   //
   // process current event
   //
   // The AliReducedAnalysisJpsi2ee ancestor class must be set via the options to run the same event pairing
   // The jpsi pair candidates will then be stored in the array fJpsiCandidates
   AliReducedAnalysisJpsi2ee::Process();

  // check event cuts
  Bool_t eventSelected = IsEventSelected(fEvent, fValues);

  // apply event selection
  if (!eventSelected) return;

  //

  // mixing
  fBmesonMixingHandler->FillEvent(&fJpsiCandidates, &fAssociatedTracks, fValues); 

  // MC
  if (fOptionRunOverMC) {
    LoopOverMCTruthTracks(1);
  }

  // loop over associated tracks
  RunAssociatedTrackSelection();
  FillAssociatedTrackHistograms();

   // run same event triplet pairing
   if(fOptionRunTripletPairing) RunSameEventTripletPairing();
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::Finish() {
  //
  // after event loop (left over mixing at some point...)
  //
   AliReducedAnalysisJpsi2ee::Finish();
   if(fOptionRunTripletPairing && fOptionRunTripletPairingMixing && !fOptionRunOverMC)
     fBmesonMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::SetTripletCut(AliReducedInfoCut* cut) {
  //
  // Set the cut on the triplets 
  // NOTE: This is just one single cut, irespective of how many jpsi or associated track cut combinations are specified
  //
  fTripletCut = cut;
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::AddAssociatedTrackCut(AliReducedInfoCut* cut) {
  //
  // add a new associated track cut
  //
  fAssociatedTrackCuts.Add(cut);

  // NOTE: number of parallel must be equal to number of hist classes in AliMixingHandler::Init()
  if (fOptionUseLikeSignPairs)  fBmesonMixingHandler->SetNParallelCuts(fBmesonMixingHandler->GetNParallelCuts()+3);
  else                          fBmesonMixingHandler->SetNParallelCuts(fBmesonMixingHandler->GetNParallelCuts()+1);
  if (fPairCuts.GetEntries())   fBmesonMixingHandler->SetNParallelPairCuts(fPairCuts.GetEntries());
  if (fPairCuts.GetEntries()>1) {
    TString histClassNamesNew  = "";
    for (Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); iPairCut++) {
      for (Int_t iTrackCut=0; iTrackCut<fAssociatedTrackCuts.GetEntries(); iTrackCut++) {
        if (fOptionUseLikeSignPairs) histClassNamesNew += Form("TripletMEPP_%s_%s_%s;",
                                                               fTrackCuts.At(iTrackCut)->GetName(),
                                                               fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                                               fPairCuts.At(iPairCut)->GetName());
        histClassNamesNew += Form("TripletMEPM_%s_%s_%s;",
                                  fTrackCuts.At(iTrackCut)->GetName(),
                                  fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                  fPairCuts.At(iPairCut)->GetName());
        if (fOptionUseLikeSignPairs) histClassNamesNew += Form("TripletMEMM_%s_%s_%s;",
                                                               fTrackCuts.At(iTrackCut)->GetName(),
                                                               fAssociatedTrackCuts.At(iTrackCut)->GetName(),
                                                               fPairCuts.At(iPairCut)->GetName());
      }
    }
    fBmesonMixingHandler->SetHistClassNames(histClassNamesNew.Data());
  } else {
    TString histClassNames = fBmesonMixingHandler->GetHistClassNames();
    if (fOptionUseLikeSignPairs) histClassNames += Form("TripletMEPP_%s_%s;", fTrackCuts.At(fAssociatedTrackCuts.GetEntries()-1)->GetName(), cut->GetName());
    histClassNames += Form("TripletMEPM_%s_%s;", fTrackCuts.At(fAssociatedTrackCuts.GetEntries()-1)->GetName(), cut->GetName());
    if (fOptionUseLikeSignPairs) histClassNames += Form("TripletMEMM_%s_%s;", fTrackCuts.At(fAssociatedTrackCuts.GetEntries()-1)->GetName(), cut->GetName());
    fBmesonMixingHandler->SetHistClassNames(histClassNames.Data());
  }
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisBmeson2JpsiK::IsMBEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
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
Bool_t AliReducedAnalysisBmeson2JpsiK::IsTripletSelected(Float_t* values) {
  //
  // apply associated track cuts
  //
  if(!fTripletCut) 
    return kTRUE;
  if (! fTripletCut->IsSelected(values))
    return kFALSE;
  else
    return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisBmeson2JpsiK::IsAssociatedTrackSelected(AliReducedBaseTrack* track, Float_t* values /*=0x0*/) {
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
void AliReducedAnalysisBmeson2JpsiK::RunAssociatedTrackSelection(Bool_t fillHistograms /*=kTRUE*/, Bool_t fillMBTracks /*=kFALSE*/) {
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
void AliReducedAnalysisBmeson2JpsiK::LoopOverAssociatedTracks(Int_t arrayOption /*=1*/, Bool_t fillHistograms /*=kTRUE*/, Bool_t fillMBTracks /*=kFALSE*/) {
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
   }  // end loop over trackss
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::RunSameEventTripletPairing(TString pairClass /*="PairSE"*/) {
  // 
  // run the same event pairing for candidates (e+e-) and the correlation to associated tracks, kaon candidate
  //
  if(fJpsiCandidates.GetEntries()==0) return;
  if(fAssociatedTracks.GetEntries()==0) return;
  TString pairTypeNames[3] = {"PP","PM","MM"};

  TIter nextAssocTrack(&fAssociatedTracks);
  TIter nextJpsi(&fJpsiCandidates);
  
  AliReducedPairInfo* jpsi = 0x0;
  AliReducedTrackInfo* assoc = 0x0;
  for(Int_t it=0;it<fJpsiCandidates.GetEntries(); ++it) {
     jpsi = (AliReducedPairInfo*)nextJpsi();

     // find full track info of jpsi candidates
     AliReducedTrackInfo* track=0x0;
     AliReducedTrackInfo* leg1=0x0;
     AliReducedTrackInfo* leg2=0x0;
     TClonesArray* trackList = fEvent->GetTracks();
     TIter nextTrack(trackList);
     // iterate over all tracks in event in order to find pointer to full track info of electron legs, break loop if pointer of both legs are found
     for(Int_t ib=0; ib<trackList->GetEntries(); ++ib) {
       if(leg1&&leg2) break;
       track = (AliReducedTrackInfo*)nextTrack();
       if(track->TrackId()==jpsi->LegId(0)) leg1=track;
       if(track->TrackId()==jpsi->LegId(1)) leg2=track;
     }

     nextAssocTrack.Reset();
     for(Int_t ia=0;ia<fAssociatedTracks.GetEntries(); ++ia) {
        assoc = (AliReducedTrackInfo*)nextAssocTrack();
        
        // make sure we do not correlate with one of the jpsi legs
        if(assoc->TrackId()==jpsi->LegId(0)) continue;
        if(assoc->TrackId()==jpsi->LegId(1)) continue;
        // NOTE: the model is that there is a set of n-selections for the jpsi and n-selections for the assoc
        //       One needs to have at least one matching bit in order to correlate them
        // TODO: we need to make sure there are equal numbers of bits for both trigger and assoc
        //           Right now, the extra bits from the particle with more bits (cuts defined) are ignored
        if(!(jpsi->GetFlags() & assoc->GetFlags())) continue;
        if (!fOptionUseLikeSignPairs && ((Int_t)jpsi->PairType() == 0 || (Int_t)jpsi->PairType() == 2)) continue;

        AliReducedVarManager::FillBcandidateInfo(jpsi, leg1, leg2, assoc, fValues);

        if (!IsTripletSelected(fValues))
          continue;
        
        // MC matching of associated (kaon candidate) track
        UInt_t decisionMap = 0;
        if(fAssociatedMCcuts.GetEntries()!=0) {
            for(Int_t i=0; i<fAssociatedMCcuts.GetEntries(); ++i) {
              AliReducedInfoCut* cut=(AliReducedInfoCut*)fAssociatedMCcuts.At(i);
              if(cut->IsSelected(assoc))
                decisionMap |= (UInt_t(1)<<i);
            }
         }

        // decide wether kaon and jpsi candidate has common mother. 
        UInt_t tripletMCdecisions = decisionMap & jpsi->GetMCFlags();
        for(Int_t imc = 0; imc<fAssociatedMCcuts.GetEntries();imc++) {
          if(!(tripletMCdecisions & (UInt_t(1)<<imc))) 
            continue;
          // fCommonAncector = index of B candidate in kaon candidates history
          Bool_t sameMotherDecision = (TMath::Abs(leg1->MCLabel(2)) == TMath::Abs(assoc->MCLabel(fCommonAncestor[imc]))); 
          // if the mother is not the same, flip bit back, i.e. empty number
          if(!sameMotherDecision)
            tripletMCdecisions ^= (UInt_t(1)<<imc);
        } 

        FillTripletHistograms(jpsi, assoc, Form("TripletSE%s", pairTypeNames[(Int_t)jpsi->PairType()].Data()), tripletMCdecisions);

     }  // end loop over associated tracks
  }  // end loop over jpsi candidates
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::LoopOverMCTruthTracks(Int_t trackArray /*=1*/) {
   //
   // loop over the track array and check the pure MC tracks against the defined MC selections
   //   
   AliReducedTrackInfo* mother=0x0;
   AliReducedTrackInfo* daughter1 = 0x0;
   AliReducedTrackInfo* daughter2 = 0x0;

   TClonesArray* trackList = (trackArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
   if(!trackList) return;
   TIter nextTrack(trackList);

   // Loop through the list of pure MC particles and find the required mothers and their corresponding daughters
   nextTrack.Reset();
   for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      mother = (AliReducedTrackInfo*)nextTrack();
 
      if(!mother->IsMCKineParticle()) continue;
      
      // apply selections on the jpsi mother
      UInt_t motherDecisions = CheckMotherMCTruth(mother); // returns 0
      
      if(!motherDecisions) continue;
      
      // find the Bmeson daughters (needed to compute 2-track properties like the polarization, etc.)
      Int_t daughter1Label = 0; Int_t daughter2Label = 0;
      FindBmesonTruthLegs(mother, daughter1Label, daughter2Label);
      
      daughter1 = AliReducedAnalysisJpsi2ee::FindMCtruthTrackByLabel(daughter1Label);
      daughter2 = AliReducedAnalysisJpsi2ee::FindMCtruthTrackByLabel(daughter2Label);
      
      // reset track variables and fill info
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i) fValues[i]=-9999.;
      AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);
      
      // loop over bmeson mother selections and fill histograms before the kine cuts on electrons
      for(Int_t iCut = 0; iCut<fBcandidateMCcuts.GetEntries(); ++iCut) {
         if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
         fHistosManager->FillHistClass(Form("%s_PureMCTruth_BeforeSelection", fBcandidateMCcuts.At(iCut)->GetName()), fValues);         
      }
      
      if(!daughter1) continue;
      if(!daughter2) continue;
      
      // apply selections on pure MC daughter (kine cuts
      // TODO: apply several cuts to one kaon and one jpsi
      UInt_t daughter1Decisions=0;
      UInt_t daughter2Decisions=0;
      CheckDaughterMCTruth(daughter1,daughter2,daughter1Decisions,daughter2Decisions);
      if(!daughter1Decisions) continue;
      //UInt_t daughtersDecisions = daughter1Decisions & CheckDaughterMCTruth(daughter2);
      if(!daughter2Decisions) continue;
 
      for(Int_t iCut = 0; iCut<fBcandidateMCcuts.GetEntries(); ++iCut) {
         if(!(motherDecisions & (UInt_t(1)<<iCut)))  continue;
         if(!(daughter1Decisions & (UInt_t(1)<<iCut)))  continue;
         if(!(daughter2Decisions & (UInt_t(1)<<iCut)))  continue;
         fHistosManager->FillHistClass(Form("%s_PureMCTruth_AfterSelection", fBcandidateMCcuts.At(iCut)->GetName()), fValues);         
      }
   }  // end loop over tracks

}

//___________________________________________________________________________
UInt_t AliReducedAnalysisBmeson2JpsiK::CheckMotherMCTruth(AliReducedTrackInfo* mother) {
   //
   // Check the mother pure MC truth against all defined selections and return a bit map with all decisions
   //
   if(fBcandidateMCcuts.GetEntries()==0) return 0;
   
   UInt_t decisionMap = 0;
   for(Int_t i=0; i<fBcandidateMCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fBcandidateMCcuts.At(i);
      // If no MC bit was requested for the mother, skip this mother signal
      // The MC cut will be applied just at the daughter level (this is likely an MC signal without a mother, e.g. electrons from gamma-gamma continuum from Starlight)
       
      //check if reweight is needed for this MC signal
      // if (checkReweight && (((AliReducedTrackCut*)cut)->GetApplyReweightMCpt())==kFALSE) continue;
      //if(!((AliReducedTrackCut*)cut)->GetMCFilterMap()) continue;
      //cout << "enters loop " << endl; // ok! 
      if(cut->IsSelected(mother)) // no output 
         decisionMap |= (UInt_t(1)<<i);
   }
   
   return decisionMap;
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::CheckDaughterMCTruth(AliReducedTrackInfo* daughter1, AliReducedTrackInfo* daughter2, UInt_t& daughter1Decisions, UInt_t& daughter2Decisions) {
   //
   // Check the daughter pure MC truth against all defined selections and return a bit map with all decisions
   //''
   if(fBcandidateDaughter1MCcuts.GetEntries()==0) return;

   daughter1Decisions = 0;
   daughter2Decisions = 0;
   for(Int_t i=0; i<fBcandidateDaughter1MCcuts.GetEntries(); ++i) {
      AliReducedInfoCut* cut = (AliReducedInfoCut*)fBcandidateDaughter1MCcuts.At(i);
      if(cut->IsSelected(daughter1))
         daughter1Decisions |= (UInt_t(1)<<i);
      cut = (AliReducedInfoCut*)fBcandidateDaughter2MCcuts.At(i);
      if(cut->IsSelected(daughter2))
         daughter2Decisions |= (UInt_t(1)<<i);
   }

}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::FindBmesonTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label) {
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

      // find jpsi leg, only first occurence stored 
      if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==443) {
         legsFound += 1;
         if(leg1Label==0) leg1Label = track->MCLabel(0);
      }

      // find kaon leg, only first occurence stored
      if ((track->MCLabel(1)==mLabel || track->MCLabel(2)==mLabel) && TMath::Abs(track->MCPdg(0))==321) {
          legsFound += 1;
          if(leg2Label==0) leg2Label = track->MCLabel(0);
      }
   }
   return;
}

//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::FillAssociatedTrackHistograms(TString trackClass /*="AssociatedTrack"*/) {
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
void AliReducedAnalysisBmeson2JpsiK::FillAssociatedTrackHistograms(AliReducedTrackInfo* track, TString trackClass  /*="AssociatedTrack"*/) {
  //
  // fill associated track histograms
  //
  // && track->IsMCTruth();
  //Bool_t isMCTruth = kFALSE;     // TODO: handle the MC info if needed
  for(Int_t icut=0; icut<fAssociatedTrackCuts.GetEntries(); ++icut) {
    if(track->TestFlag(icut)) {
      fHistosManager->FillHistClass(Form("%s_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName()), fValues);
    }
    if(fOptionRunOverMC){
        for(Int_t i=0; i<fAssociatedMCcuts.GetEntries(); ++i) {
          AliReducedInfoCut* cut=(AliReducedInfoCut*)fAssociatedMCcuts.At(i);
          if(cut->IsSelected(track))
            fHistosManager->FillHistClass(Form("%s_%s_%s", trackClass.Data(), fAssociatedTrackCuts.At(icut)->GetName(), cut->GetName()), fValues);
          }
    }
  }
}



//___________________________________________________________________________
void AliReducedAnalysisBmeson2JpsiK::FillTripletHistograms(AliReducedPairInfo* jpsi, AliReducedBaseTrack* assoc, TString corrClass /*="CorrSE"*/, UInt_t tripletMCdecisions/*=0*/) {
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
          //if(isMCTruth) fHistosManager->FillHistClass(Form("%s_%s_%s_MCTruth", corrClass.Data(),
          //                                                 fTrackCuts.At(iTrackCut)->GetName(),
          //                                                 fAssociatedTrackCuts.At(iTrackCut)->GetName(),
          //                                                 fPairCuts.At(iPairCut)->GetName()), fValues);
        }
      }
    }
  } else {
    for(Int_t iCut=0; iCut<fAssociatedTrackCuts.GetEntries(); ++iCut) {
       if(!(trackMask & (ULong_t(1)<<iCut))) continue;
       fHistosManager->FillHistClass(Form("%s_%s_%s", corrClass.Data(), fTrackCuts.At(iCut)->GetName(), fAssociatedTrackCuts.At(iCut)->GetName()), fValues);
       // if kaon and jpsi is not from the same mother, continue: else, fill MC histograms
       if(!tripletMCdecisions)
         continue;
       for(Int_t imc = 0; imc<fAssociatedTrackCuts.GetEntries();imc++) {
          if(!(tripletMCdecisions & (UInt_t(1)<<imc)))
            continue;
          fHistosManager->FillHistClass(Form("%s_%s_%s_%s_%s", corrClass.Data(), fTrackCuts.At(iCut)->GetName(), fAssociatedTrackCuts.At(iCut)->GetName(), fLegCandidatesMCcuts.At(imc)->GetName(), fAssociatedMCcuts.At(imc)->GetName()), fValues);
      }
    }
  }
}
