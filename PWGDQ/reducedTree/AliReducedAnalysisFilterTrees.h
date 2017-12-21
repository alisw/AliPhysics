//
// Creation date: 2017/08/09
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISFILTERTREES_H
#define ALIREDUCEDANALYSISFILTERTREES_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

//________________________________________________________________
class AliReducedAnalysisFilterTrees : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisFilterTrees();
  AliReducedAnalysisFilterTrees(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisFilterTrees();
  
  virtual void Init();
  virtual void Process();
  virtual void Finish();
  
  // setters
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrackCut(AliReducedInfoCut* cut) {fTrackCuts.Add(cut);}
  void AddPairCut(AliReducedInfoCut* cut) {fPairCuts.Add(cut);}
  void SetWriteFilteredTracks(Bool_t option=kTRUE) {fWriteFilteredTracks=option;}
  void SetWriteFilteredPairs(Bool_t option=kTRUE) {fWriteFilteredPairs=option;}
  
  void SetBuildCandidatePairs(AliReducedPairInfo::CandidateType type) {fBuildCandidatePairs=kTRUE; fCandidateType=type;}
  void SetBuildCandidateLikePairs(Bool_t option=kTRUE) {fBuildCandidateLikePairs=option;}
  void AddCandidateLeg1Cut(AliReducedInfoCut* cut) {fLeg1Cuts.Add(cut);}
  void AddCandidateLeg2Cut(AliReducedInfoCut* cut) {fLeg2Cuts.Add(cut);}
  void AddCandidatePairCut(AliReducedInfoCut* cut) {fCandidatePairCuts.Add(cut);}
  void SetRunCandidatePrefilter(Bool_t option=kTRUE) {fRunCandidatePrefilter=option;}
  void SetRunCandidatePrefilterOnSameCharge(Bool_t option=kTRUE) {fRunCandidatePrefilterOnSameCharge=option;}
  void AddCandidateLeg1PrefilterCut(AliReducedInfoCut* cut) {fLeg1PrefilterCuts.Add(cut);}
  void AddCandidateLeg2PrefilterCut(AliReducedInfoCut* cut) {fLeg2PrefilterCuts.Add(cut);}
  void AddCandidateLeg1PairPrefilterCut(AliReducedInfoCut* cut) {fLeg1PairPrefilterCuts.Add(cut);}
  void AddCandidateLeg2PairPrefilterCut(AliReducedInfoCut* cut) {fLeg2PairPrefilterCuts.Add(cut);}
  
  // getters
  virtual AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  Bool_t GetWriteFilteredTracks() const {return fWriteFilteredTracks;}
  Int_t GetNTrackCuts() const {return fTrackCuts.GetEntries();}
  const Char_t* GetTrackCutName(Int_t i) const {return (i<fTrackCuts.GetEntries() ? fTrackCuts.At(i)->GetName() : "");} 
  Bool_t GetWriteFilteredPairs() const {return fWriteFilteredPairs;}
  Int_t GetNPairCuts() const {return fPairCuts.GetEntries();}
  const Char_t* GetPairCutName(Int_t i) const {return (i<fPairCuts.GetEntries() ? fPairCuts.At(i)->GetName() : "");} 
  Bool_t GetBuildCandidatePairs() const {return fBuildCandidatePairs;}
  Bool_t GetBuildCandidateLikePairs() const {return fBuildCandidateLikePairs;}
  Int_t GetCandidateType() const {return fCandidateType;}
  Bool_t GetRunCandidatePrefilter() const {return fRunCandidatePrefilter;}
  Int_t GetNCandidateLegCuts() const {return fLeg1Cuts.GetEntries();}
  const Char_t* GetCandidateLegCutName(Int_t i, Int_t leg);
  Bool_t IsAsymmetricDecayChannel();
  
protected:
   AliHistogramManager* fHistosManager;   // Histogram manager
   
   TList fEventCuts;               // array of event cuts used for filtering
   TList fTrackCuts;               // array of track cuts used for filtering
   Bool_t fWriteFilteredTracks;   // filter the track list
   TList fPairCuts;                  // array of pair cuts used for filtering
   Bool_t fWriteFilteredPairs;   // filter the pair list
   
   Bool_t fBuildCandidatePairs;   // if true, build additional candidate pairs from selected tracks 
   Bool_t fBuildCandidateLikePairs;  // if true, build also like pairs (e.g. like-sign for symmetric decay channels)
   Int_t fCandidateType;             // candidate type, see AliReducedPairInfo::CandidateType 
   TList fLeg1Cuts;                      // list of track cuts for LEG1  (these cuts will be used also for LEG2 if the decay channel is symmetric)
   TList fLeg2Cuts;                      // list of tracks cuts for LEG2 (NOTE: fLeg1Cuts and fLeg2Cuts must contain the same number of cuts)
   TList fCandidatePairCuts;          // list of cuts for pair candidates
   Bool_t fRunCandidatePrefilter;   // if true, run a prefilter on the selected legs
   Bool_t fRunCandidatePrefilterOnSameCharge;   // default FALSE (unlike charged pairs only);
                                                                // if true, run the prefilter on same charge pairs also;
   TList fLeg1PrefilterCuts;            // cuts for tracks used in the prefilter for LEG1  
   TList fLeg2PrefilterCuts;            // cuts for tracks used in the prefilter for LEG2
   TList fLeg1PairPrefilterCuts;      // cuts on the prefilter pairs for LEG1
   TList fLeg2PairPrefilterCuts;      // cuts on the prefilter pairs for LEG2
   TList fLeg1Tracks;                    // list of selected LEG1 tracks in the current event
   TList fLeg2Tracks;                    // list of selected LEG2 tracks in the current event
   TList fLeg1PrefilteredTracks;    // list of prefilter selected LEG1 tracks in the current event
   TList fLeg2PrefilteredTracks;    // list of prefilter selected LEG2 tracks in the current event
   
   Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
   Bool_t IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
   Bool_t IsPairSelected(AliReducedPairInfo* pair, Float_t* values=0x0);
   void CreateFilteredEvent();
   Bool_t TrackIsCandidateLeg(AliReducedBaseTrack* track);
   void WriteFilteredPairs();
   void WriteFilteredTracks();
   Bool_t IsCandidateLegSelected(AliReducedBaseTrack* track, Float_t* values=0x0, Int_t whichLeg=1); 
   Bool_t IsCandidatePairSelected(Float_t* values);
   Bool_t IsCandidateLegPrefilterSelected(AliReducedBaseTrack* track, Float_t* values=0x0, Int_t whichLeg=1);
   Bool_t IsCandidateLegPairPrefilterSelected(Float_t* values, Int_t whichLeg=1);
   void BuildCandidatePairs();
   void RunCandidateLegsSelection();
   void RunCandidateLegsPrefilter(Int_t leg);
   void RunSameEventPairing();
   void SetupPair(AliReducedPairInfo* pair, Float_t* values);
   ULong_t CheckTrackCompatibility(AliReducedBaseTrack* leg1, AliReducedBaseTrack* leg2, Bool_t isAsymmetricDecayChannel);
   void FillCandidateLegHistograms(TString histClass, AliReducedBaseTrack* track, Float_t* values, Int_t leg, Bool_t isAsymmetricDecayChannel);
   void FillCandidatePairHistograms(TString histClass, AliReducedPairInfo* pair, Float_t* values, Bool_t isAsymmetricDecayChannel);
   
  ClassDef(AliReducedAnalysisFilterTrees,1);
};

#endif
