//
// Creation date: 2016/09/06
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#ifndef ALIREDUCEDANALYSISPSI2EE_H
#define ALIREDUCEDANALYSISPSI2EE_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedCaloClusterTrackMatcher.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"


//________________________________________________________________
class AliReducedAnalysisPsi2ee : public AliReducedAnalysisTaskSE, public AliAnalysisTaskSE {
  
public:
  AliReducedAnalysisPsi2ee();
  AliReducedAnalysisPsi2ee(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisPsi2ee();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();
  
  // setters
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddClusterCut(AliReducedInfoCut* cut) {fClusterCuts.Add(cut); fFillCaloClusterHistograms=kTRUE; }
  void AddTrackCut(AliReducedInfoCut* cut);
  void AddPrefilterTrackCut(AliReducedInfoCut* cut) {fPreFilterTrackCuts.Add(cut);}
  void AddPairCut(AliReducedInfoCut* cut);
  void AddPrefilterPairCut(AliReducedInfoCut* cut) {fPreFilterPairCuts.Add(cut);}
  void SetRunEventMixing(Bool_t option) {fOptionRunMixing = option;};
  void SetRunPairing(Bool_t option) {fOptionRunPairing = option;};
  void SetRunOverMC(Bool_t option) {fOptionRunOverMC = option;};
  void SetRunLikeSignPairing(Bool_t option) {fOptionRunLikeSignPairing = option;}
  void SetLoopOverTracks(Bool_t option) {
     fOptionLoopOverTracks = option; 
     if(!fOptionLoopOverTracks) {fOptionRunPairing = kFALSE; fOptionRunMixing = kFALSE; fOptionRunLikeSignPairing = kFALSE;}     
  }
  void SetRunPrefilter(Bool_t option) {fOptionRunPrefilter = option;}
  void SetStoreccbarCandidates(Bool_t option) {fOptionStoreccbarCandidates = option;}
  void SetMCJpsiPtWeights(TH1F* weights) {fMCJpsiPtWeights = weights;}
  void SetMCPsiPtWeights(TH1F* weights) {fMCPsiPtWeights = weights;}
  void SetFillCaloClusterHistograms(Bool_t option) {fFillCaloClusterHistograms = option;}
  void SetClusterTrackMatcher(AliReducedCaloClusterTrackMatcher* matcher) {fClusterTrackMatcher = matcher;}
  void SetStoreTree(Bool_t option) {fOptionStoreTree = option;}

  void AddLegCandidateMCcut(AliReducedInfoCut* cut, Bool_t sameMother=kTRUE) {
     if(fLegCandidatesMCcuts.GetEntries()>=32) return;
     fLegCandidatesMCcuts.Add(cut);
     fLegCandidatesMCcuts_RequestSameMother[fLegCandidatesMCcuts.GetEntries()-1] = sameMother;
  }
  void AddJpsiMotherMCCut(AliReducedInfoCut* cutMother, AliReducedInfoCut* cutElectron) {
     if(fJpsiMotherMCcuts.GetEntries()>=32) return;
     fJpsiMotherMCcuts.Add(cutMother);
     fJpsiElectronMCcuts.Add(cutElectron);
  }
  void AddPsiMotherMCCut(AliReducedInfoCut* cutMother, AliReducedInfoCut* cutElectron) {
     if(fPsiMotherMCcuts.GetEntries()>=32) return;
     fPsiMotherMCcuts.Add(cutMother);
     fPsiElectronMCcuts.Add(cutElectron);
  }
  void DefineTree();
  
  // getters
  virtual AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  virtual AliMixingHandler* GetMixingHandler() const {return fMixingHandler;}
  virtual AliReducedCaloClusterTrackMatcher* GetClusterTrackMatcher() const {return fClusterTrackMatcher;}
  Int_t GetNClusterCuts() const {return fClusterCuts.GetEntries();}
  const Char_t* GetClusterCutName(Int_t i) const {return (i<fClusterCuts.GetEntries() ? fClusterCuts.At(i)->GetName() : "");}
  AliReducedInfoCut* GetCaloClusterCut(Int_t i) const {return (i<fClusterCuts.GetEntries() ? (AliReducedInfoCut*)fClusterCuts.At(i) : NULL);}
  Int_t GetNTrackCuts() const {return fTrackCuts.GetEntries();}
  const Char_t* GetTrackCutName(Int_t i) const {return (i<fTrackCuts.GetEntries() ? fTrackCuts.At(i)->GetName() : "");}
  AliReducedInfoCut* GetTrackCut(Int_t i) const {return (i<fTrackCuts.GetEntries() ? (AliReducedInfoCut*)fTrackCuts.At(i) : NULL);}
  Int_t GetNPairCuts() const {return fPairCuts.GetEntries();}
  const Char_t* GetPairCutName(Int_t i) const {return (i<fPairCuts.GetEntries() ? fPairCuts.At(i)->GetName() : "");}
  AliReducedInfoCut* GetPairCut(Int_t i) const {return (i<fPairCuts.GetEntries() ? (AliReducedInfoCut*)fPairCuts.At(i) : NULL);}
  Bool_t GetRunOverMC() const {return fOptionRunOverMC;};
  Bool_t GetRunLikeSignPairing() const {return fOptionRunLikeSignPairing;}
  Bool_t GetRunEventMixing() const {return fOptionRunMixing;}
  Bool_t GetRunPairing() const {return fOptionRunPairing;}
  Bool_t GetLoopOverTracks() const {return fOptionLoopOverTracks;}
  Bool_t GetRunPrefilter() const {return fOptionRunPrefilter;}
  Bool_t GetStoreccbarCandidates() const {return fOptionStoreccbarCandidates;}
  Bool_t GetFillCaloClusterHistograms() const {return fFillCaloClusterHistograms;}
  Int_t GetNLegCandidateMCcuts() const {return fLegCandidatesMCcuts.GetEntries();}
  const Char_t* GetLegCandidateMCcutName(Int_t i) const {return (i<fLegCandidatesMCcuts.GetEntries() ? fLegCandidatesMCcuts.At(i)->GetName() : "");}
  Int_t GetNJpsiMotherMCCuts() const {return fJpsiMotherMCcuts.GetEntries();}
  Int_t GetNPsiMotherMCCuts() const {return fPsiMotherMCcuts.GetEntries();}
  const Char_t* GetJpsiMotherMCcutName(Int_t i) const {return (i<fJpsiMotherMCcuts.GetEntries() ? fJpsiMotherMCcuts.At(i)->GetName() : "");}
  const Char_t* GetPsiMotherMCcutName(Int_t i) const {return (i<fPsiMotherMCcuts.GetEntries() ? fPsiMotherMCcuts.At(i)->GetName() : "");}

  
protected:
   AliHistogramManager*               fHistosManager;       // Histogram manager
   AliMixingHandler*                  fMixingHandler;       // mixing handler
   AliReducedCaloClusterTrackMatcher* fClusterTrackMatcher; // cluster-track matcher
   
   Bool_t fOptionRunMixing;    // true: run event mixing, false: no event mixing
   Bool_t fOptionRunPairing;    // true: run pairing, false: only apply the track cuts
   Bool_t fOptionRunOverMC;  // true: trees contain MC info -> fill histos to compute efficiencies, false: run normally as on data
   Bool_t fOptionRunLikeSignPairing;   // true (default): performs the like sign pairing in addition to the opposite pairing
   Bool_t fOptionLoopOverTracks;       // true (default); if false do not loop over tracks and consequently no pairing
   Bool_t fOptionRunPrefilter;        // true (default); if false do not run the prefilter
   Bool_t fOptionStoreccbarCandidates;   // false (default); if true, store the same event jpsi candidates in a TList 
   Bool_t fFillCaloClusterHistograms;   // false (default); if true, fill calorimeter cluster histograms
   Bool_t fOptionStoreTree; 

   TList fEventCuts;               // array of event cuts
   TList fClusterCuts;             // array of cluster cuts
   TList fTrackCuts;               // array of track cuts
   TList fPreFilterTrackCuts;  // track cuts to be used at the prefilter stage
   TList fPairCuts;                  // array of pair cuts
   TList fPreFilterPairCuts;     // pair cuts to be used at the prefilter stage

   TList fClusters;               // list of selected clusters
   TList fPosTracks;               // list of selected positive tracks in the current event
   TList fNegTracks;              // list of selected negative tracks in the current event
   TList fPrefilterPosTracks;  // list of prefilter selected positive tracks in the current event
   TList fPrefilterNegTracks; // list of prefilter selected negative tracks in the current event
   TList fccbarCandidates;       // list of Jpsi candidates --> to be used in analyses inheriting from this 
   
   // selection based on the MC truth information of the reconstructed leg candidates
   // NOTE:    The list is a list of AliReducedInfoCut objects which can be used to 
   //              apply cuts on the MC flags of the tracks.
   // NOTE: The names of the cuts are used in the naming of the histogram classes
   TList fLegCandidatesMCcuts;
   Bool_t fLegCandidatesMCcuts_RequestSameMother[32];
   
   // selection cuts for the pure MC truth (select the J/psi from stack)
   // the list should contains cuts which can be applied to a pure MC truth particle (no reconstructed information)
   //  e.g. cuts on the MC flags and on kinematics
   //  For each selection, a separate histogram directory will be created
   TList fJpsiMotherMCcuts;
   TList fPsiMotherMCcuts;
   
   // Selection on the MC truth of the electrons from the jpsi decay
   //  Tipically, here one can specify the kinematic selection on the electrons from jpsi decay
   //       so dividing the jpsi yield at this step by the yield of jpsi selected by the fJpsiMotherMCcuts, one can obtain the
   //       acceptance efficiency.
   //  NOTE: The number of selections on the jpsi electron needs to be the same and in sync with the number of fJpsiMotherMCcuts cuts
   TList fJpsiElectronMCcuts;
   TList fPsiElectronMCcuts;

   TTree* fPairTreeOutput;
   TTree* fLeg1TreeOutput;
   TTree* fLeg2TreeOutput;
   TTree* fMCTruthTreeOutput;

  Float_t* fTreeVariable;
  Float_t* fPairTreeVariable;
  Float_t* fLeg1TreeVariable;
  Float_t* fLeg2TreeVariable;
  Float_t* fMCTruthTreeVariable;
  std::vector <TString> motherTags = {"MCRecJpsiElectron", "MCRecJpsiElectronNonPrompt", "MCRecJpsiElectronPrompt", "MCRecPsiElectron", "MCRecPsiElectronNonPrompt", "MCRecPsiElectronPrompt"};
   
  Float_t fValuesTrk1[AliReducedVarManager::kNVars] = {-9999.,};
  Float_t fValuesTrk2[AliReducedVarManager::kNVars] = {-9999.,};
  
  Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t IsClusterSelected(AliReducedCaloClusterInfo* cluster, Float_t* values=0x0);
  Bool_t IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsTrackPrefilterSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  ULong_t IsPairSelected(Float_t* values);
  Bool_t IsPairPreFilterSelected(Float_t* values);
  UInt_t CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack);
  UInt_t CheckReconstructedLegMCTruth(AliReducedBaseTrack* track);
  void    FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label, Int_t& leg2Label);
  AliReducedTrackInfo* FindMCtruthTrackByLabel(Int_t label);
  void    LoopOverMCTracks(Int_t trackArray =1);
  UInt_t CheckMotherMCTruth(AliReducedTrackInfo* mother, Bool_t checkReweight=kFALSE, Int_t type = AliReducedPairInfo::kJpsiToEE);
  UInt_t CheckDaughterMCTruth(AliReducedTrackInfo* daughter, Int_t type = AliReducedPairInfo::kJpsiToEE);
  
  void RunPrefilter();
  void RunSameEventPairing(TString pairClass = "PairSE");
  void RunTrackSelection();
  void RunClusterSelection();
  void LoopOverTracks(Int_t arrayOption=1);
  void FillTrackHistograms(TString trackClass = "Track");
  void FillTrackHistograms(AliReducedBaseTrack* track, TString trackClass = "Track");
  void FillPairHistograms(ULong_t trackMask, ULong_t pairMask, Int_t pairType, TString pairClass = "PairSE", UInt_t mcDecisions = 0);
  void FillClusterHistograms(TString clusterClass="CaloCluster");
  void FillClusterHistograms(AliReducedCaloClusterInfo* cluster, TString clusterClass="CaloCluster");
  void FillMCTruthHistograms();
  void FillMCTruthTree(Float_t* fValues, Int_t type);
  void FillLegTree(Float_t* fValues, Int_t trkTag = 1, UInt_t mcDecisions = 0);
  void FillPairTree(Float_t* fValues, UInt_t mcDecisions = 0);

  TList*  fClusterTrackMatcherHistograms;             // list of cluster-track matcher histograms
  TH1I*   fClusterTrackMatcherMultipleMatchesBefore;  // multiple matches of tracks to same cluster before matching
  TH1I*   fClusterTrackMatcherMultipleMatchesAfter;   // multiple matches of tracks to same cluster after matching

  Bool_t fSkipMCEvent;          // decision to skip MC event
  TH1F*  fMCJpsiPtWeights;            // weights vs pt to reject events depending on the jpsi true pt (needed to re-weights jpsi Pt distribution)
  TH1F*  fMCPsiPtWeights;            // weights vs pt to reject events depending on the jpsi true pt (needed to re-weights psi Pt distribution)
  
  ClassDef(AliReducedAnalysisPsi2ee,13);
};

#endif
