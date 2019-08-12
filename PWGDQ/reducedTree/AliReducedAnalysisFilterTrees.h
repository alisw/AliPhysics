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
  void SetMCJpsiPtWeights(TH1F* weights) {fMCJpsiPtWeights = weights;}
  
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
  
  void SetRunOverMC(Bool_t option) {fOptionRunOverMC = option;};
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
   Bool_t GetRunOverMC() const {return fOptionRunOverMC;};
  //  Int_t GetNLegCandidateMCcuts() const {return fLegCandidatesMCcuts.GetEntries();}
  //const Char_t* GetLegCandidateMCcutName(Int_t i) const {return (i<fLegCandidatesMCcuts.GetEntries() ? fLegCandidatesMCcuts.At(i)->GetName() : "");}
  const Char_t* GetLegCandidateMCcutName() const {return (fLegCandidatesMCcuts ? fLegCandidatesMCcuts->GetName() : "");}
  Int_t GetNJpsiMotherMCCuts() const {return fJpsiMotherMCcuts.GetEntries();}
  const Char_t* GetJpsiMotherMCcutName(Int_t i) const {return (i<fJpsiMotherMCcuts.GetEntries() ? fJpsiMotherMCcuts.At(i)->GetName() : "");}

 
 /* void AddLegCandidateMCcut(AliReducedInfoCut* cut) {
     if(fLegCandidatesMCcuts.GetEntries()>=32) return;
     fLegCandidatesMCcuts.Add(cut);
  }*/
 void SetLegCandidateMCcut(AliReducedInfoCut* cut) {
     //if(fLegCandidatesMCcuts.GetEntries()>=32) return;
     fLegCandidatesMCcuts  = cut;
  }

  void AddJpsiMotherMCCut(AliReducedInfoCut* cutMother, AliReducedInfoCut* cutElectron) {
     if(fJpsiMotherMCcuts.GetEntries()>=32) return;
     fJpsiMotherMCcuts.Add(cutMother);
     fJpsiElectronMCcuts.Add(cutElectron);
  }

  void FillMCTruthHistograms();
 
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

   Bool_t fOptionRunOverMC;  // true: trees contain MC info -> fill histos to compute efficiencies, false: run normally as on data
   // selection based on the MC truth information of the reconstructed leg candidates
   // NOTE:    The list is a list of AliReducedInfoCut objects which can be used to 
   //              apply cuts on the MC flags of the tracks.
   // NOTE: The names of the cuts are used in the naming of the histogram classes
   AliReducedInfoCut *fLegCandidatesMCcuts; 
 
   // selection cuts for the pure MC truth (select the J/psi from stack)
   // the list should contains cuts which can be applied to a pure MC truth particle (no reconstructed information)
   //  e.g. cuts on the MC flags and on kinematics
   //  For each selection, a separate histogram directory will be created
   TList fJpsiMotherMCcuts;
   TH1F*  fMCJpsiPtWeights;            //! weights vs pt to reject events depending on the jpsi true pt (needed to re-weights jpsi Pt distribution)
   Bool_t fSkipMCEvent; // if true MC event is skipped
   // Selection on the MC truth of the electrons from the jpsi decay
   //  Tipically, here one can specify the kinematic selection on the electrons from jpsi decay
   //       so dividing the jpsi yield at this step by the yield of jpsi selected by the fJpsiMotherMCcuts, one can obtain the
   //       acceptance efficiency.
   //  NOTE: The number of selections on the jpsi electron needs to be the same and in sync with the number of fJpsiMotherMCcuts cuts
   TList fJpsiElectronMCcuts;
   
   Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
   Bool_t IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
   Bool_t IsPairSelected(AliReducedPairInfo* pair, Float_t* values=0x0);
   void CreateFilteredEvent();
   Bool_t CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack);
 Bool_t CheckReconstructedLegMCTruth(AliReducedBaseTrack* track);
  void    FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label);
  AliReducedTrackInfo* FindMCtruthTrackByLabel(Int_t label);
  void    LoopOverMCTracks(Int_t trackArray =1);
  UInt_t CheckMotherMCTruth(AliReducedTrackInfo* mother);
  UInt_t CheckDaughterMCTruth(AliReducedTrackInfo* daughter); 

  Bool_t TrackIsCandidateLeg(AliReducedBaseTrack* track);
   void WriteFilteredPairs();
   void WriteFilteredTracks(Int_t array=1);
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
