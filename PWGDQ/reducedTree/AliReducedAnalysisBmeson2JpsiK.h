//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Analysis task for J/psi - hadron correlations                        //
//                                                                      //
// creation date: 19/08/2021                                            //
// authors:       I. Arsene, ionut.arsene@cern.ch                       //
//                I. Storehaug, ida.storehaug@cern.ch                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ALIREDUCEDANALYSISBMESON2JPSIK_H
#define ALIREDUCEDANALYSISBMESON2JPSIK_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedAnalysisJpsi2ee.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"

//________________________________________________________________
class AliReducedAnalysisBmeson2JpsiK : public AliReducedAnalysisJpsi2ee {

public:
  AliReducedAnalysisBmeson2JpsiK();
  AliReducedAnalysisBmeson2JpsiK(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisBmeson2JpsiK();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();

  // setters
  void AddMBEventCut(AliReducedInfoCut* cut) {fMBEventCuts.Add(cut);}
  void AddAssociatedTrackCut(AliReducedInfoCut* cut);
  void AddJpsiElectronCut(AliReducedInfoCut* cut) {AddTrackCut(cut);}      // synonim function to AddTrackCut
  void SetUseLikeSignPairs(Bool_t option) {fOptionUseLikeSignPairs = option;}
  void SetRunTripletPairing(Bool_t option) {fOptionRunTripletPairing = option;}
  void SetLoopOverTracks(Bool_t option) {
    fOptionLoopOverTracks = option;
    if(!fOptionLoopOverTracks) {fOptionRunPairing = kFALSE; fOptionRunMixing = kFALSE; fOptionUseLikeSignPairs = kFALSE; fOptionRunTripletPairing = kFALSE;}
  }
  void SetRunTripletPairingMixing(Bool_t option) {fOptionRunTripletPairingMixing = option;};

  void AddAssocMCcut(AliReducedInfoCut* cut, Int_t commonAncestor) {
    if(fAssociatedMCcuts.GetEntries()>=32) return;
    fAssociatedMCcuts.Add(cut);
    fCommonAncestor.push_back(commonAncestor);
  }

  void AddBMotherMCCut(AliReducedInfoCut* cutMother, AliReducedInfoCut* cutJpsi, AliReducedInfoCut* cutKaon) {
     if(fBcandidateMCcuts.GetEntries()>=32) return;
     fBcandidateMCcuts.Add(cutMother);
     fBcandidateDaughter1MCcuts.Add(cutJpsi);
     fBcandidateDaughter2MCcuts.Add(cutKaon);
  }

  void SetTripletCut(AliReducedInfoCut* cut); 

  void    LoopOverMCTruthTracks(Int_t trackArray =1);
  UInt_t  CheckMotherMCTruth(AliReducedTrackInfo* mother);
  void    CheckDaughterMCTruth(AliReducedTrackInfo* daughter1, AliReducedTrackInfo* daughter2, UInt_t& daughter1Decisions, UInt_t& daughter2Decisions);
  void    FindBmesonTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label);
  

  // getters
  Int_t GetNAssociatedTrackCuts() const {return fAssociatedTrackCuts.GetEntries();}
  const Char_t* GetAssociatedTrackCutName(Int_t i) const {return (i<fAssociatedTrackCuts.GetEntries() ? fAssociatedTrackCuts.At(i)->GetName() : "");}
  Bool_t GetRunTripletPairing() {return fOptionRunTripletPairing;}
  AliMixingHandler* GetTripletPairingMixingHandler() const {return fBmesonMixingHandler;};
  Bool_t GetUseLikeSignPairs() const {return fOptionUseLikeSignPairs;}
  Bool_t GetRunTripletPairingMixing() const {return fOptionRunTripletPairingMixing;}
  Int_t GetNMCAssocCuts() const {return fAssociatedMCcuts.GetEntries();}
  const Char_t* GetAssocMCcutName(Int_t i) const {return (i<fAssociatedMCcuts.GetEntries() ? fAssociatedMCcuts.At(i)->GetName() : "");}
  Int_t GetNBcandidateMCcuts() const {return fBcandidateMCcuts.GetEntries();}
  const Char_t* GetBcandidateMCcutsName(Int_t i) const {return (i<fBcandidateMCcuts.GetEntries() ? fBcandidateMCcuts.At(i)->GetName() : "");}
  
 
protected:
  AliMixingHandler* fBmesonMixingHandler;

  Bool_t fOptionUseLikeSignPairs;     // true: use like sign pairs for correlation
  Bool_t fOptionRunTripletPairing;       // true: run correlation of candidates and associated tracks, false: apply only associated track cuts
  Bool_t fOptionRunTripletPairingMixing; // true: run the event mixing for the correlation function

  TList fMBEventCuts;           // array of MB event cuts (for mixing with EMCal trigger)
  TList fAssociatedTrackCuts;   // array of associated track cuts
  TList fAssociatedTracks;      // list of selected associated tracks
  TList fAssociatedTracksMB;    // list of selected associated tracks in MB events (for mixing with EMCal trigger)
  TList fAssociatedMCcuts;      // list of MC truth selesctions for associated track
  std::vector<Int_t> fCommonAncestor;  // list of common ancestors

  TList fBcandidateMCcuts;
  TList fBcandidateDaughter1MCcuts;
  TList fBcandidateDaughter2MCcuts;

  AliReducedInfoCut* fTripletCut;    // selection cut on the triplets

  Bool_t IsMBEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t IsAssociatedTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsTripletSelected(Float_t* values=0x0);

  void RunAssociatedTrackSelection(Bool_t fillHistograms=kTRUE, Bool_t fillMBTracks=kFALSE);
  void LoopOverAssociatedTracks(Int_t arrayOption=1, Bool_t fillHistograms=kTRUE, Bool_t fillMBTracks=kFALSE);
  void RunSameEventTripletPairing(TString pairClass="PairSE");
  void FillAssociatedTrackHistograms(TString trackClass = "AssociatedTrack");
  void FillAssociatedTrackHistograms(AliReducedTrackInfo* track, TString trackClass = "AssociatedTrack");
  void FillTripletHistograms(AliReducedPairInfo* jpsi, AliReducedBaseTrack* assoc, TString corrClass="CorrSE", UInt_t tripletMCdecisions=0);

  ClassDef(AliReducedAnalysisBmeson2JpsiK, 1);
};

#endif
