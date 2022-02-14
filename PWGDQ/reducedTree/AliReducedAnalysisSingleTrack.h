///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Analysis task for single tracks                                       //
//                                                                       //
// creation date: 10/10/2018                                             //
// author: Lucas Altenkamper, lucas.altenkamper@cern.ch                  //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIREDUCEDANALYSISSINGLETRACK_H
#define ALIREDUCEDANALYSISSINGLETRACK_H

#include <TList.h>

#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedCaloClusterTrackMatcher.h"
#include "AliHistogramManager.h"
#include "AliMixingHandler.h"

//________________________________________________________________
class AliReducedAnalysisSingleTrack : public AliReducedAnalysisTaskSE {
  
public:
  AliReducedAnalysisSingleTrack();
  AliReducedAnalysisSingleTrack(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisSingleTrack();
  
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();

  // setters
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrackCut(AliReducedInfoCut* cut) {fTrackCuts.Add(cut);}
  void AddClusterCut(AliReducedInfoCut* cut) {fClusterCuts.Add(cut); fOptionRunOverCaloCluster = kTRUE;}
  void AddMCSignalCut(AliReducedInfoCut* cut) {if (fMCSignalCuts.GetEntries()>=32) return; fMCSignalCuts.Add(cut);}
  void SetRunOverMC(Bool_t option) {fOptionRunOverMC = option;}
  void SetRunOverCaloCluster(Bool_t option) {fOptionRunOverCaloCluster = option;}
  void SetClusterTrackMatcher(AliReducedCaloClusterTrackMatcher* matcher) {fClusterTrackMatcher = matcher;}

  //getters
  virtual AliHistogramManager*  GetHistogramManager() const {return fHistosManager;}
  virtual AliReducedCaloClusterTrackMatcher* GetClusterTrackMatcher() const {return fClusterTrackMatcher;}
  Int_t                         GetNTrackCuts() const {return fTrackCuts.GetEntries();}
  const Char_t*                 GetTrackCutName(Int_t i) const {return (i<fTrackCuts.GetEntries() ? fTrackCuts.At(i)->GetName() : "");}
  Bool_t                        GetRunOverCaloCluster() const {return fOptionRunOverCaloCluster;};
  Int_t                         GetNClusterCuts() const {return fClusterCuts.GetEntries();}
  const Char_t*                 GetClusterCutName(Int_t i) const {return (i<fClusterCuts.GetEntries() ? fClusterCuts.At(i)->GetName() : "");}
  Bool_t                        GetRunOverMC() const {return fOptionRunOverMC;};
  Int_t                         GetNMCSignalCuts() const {return fMCSignalCuts.GetEntries();}
  const Char_t*                 GetMCSignalCutName(Int_t i) const {return (i<fMCSignalCuts.GetEntries() ? fMCSignalCuts.At(i)->GetName() : "");}
  
protected:
  AliHistogramManager*                fHistosManager;         // Histogram manager
  AliReducedCaloClusterTrackMatcher*  fClusterTrackMatcher;   // cluster-track matcher
  
  Bool_t fOptionRunOverMC;          // true: trees contain MC info -> fill histos to compute efficiencies, false: run normally as on data
  Bool_t fOptionRunOverCaloCluster; // true: run over calorimeter clusters and fill histograms

  TList fEventCuts;     // array of event cuts
  TList fTrackCuts;     // array of track cuts
  TList fClusterCuts;   // array of cluster cuts
  TList fMCSignalCuts;  // array of MC signal cuts

  TList fTracks;        // list of selected tracks
  TList fClusters;      // list of selected clusters

  TList*  fClusterTrackMatcherHistograms;             // list of cluster-track matcher histograms
  TH1I*   fClusterTrackMatcherMultipleMatchesBefore;  // multiple matches of tracks to same cluster before matching
  TH1I*   fClusterTrackMatcherMultipleMatchesAfter;   // multiple matches of tracks to same cluster after matching

  Bool_t  IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);
  Bool_t  IsTrackSelected(AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t  IsClusterSelected(AliReducedCaloClusterInfo* cluster, Float_t* values=0x0);
  UInt_t  CheckTrackMCTruth(AliReducedBaseTrack* track);
  void    FillMCTruthHistograms();
  void    RunTrackSelection();
  void    RunClusterSelection();
  void    FillTrackHistograms(TString trackClass="Track");
  void    FillTrackHistograms(AliReducedBaseTrack* track, TString trackClass="Track");
  void    FillClusterHistograms(TString clusterClass="CaloCluster");
  void    FillClusterHistograms(AliReducedCaloClusterInfo* cluster, TString clusterClass="CaloCluster");

  ClassDef(AliReducedAnalysisSingleTrack,3);
};

#endif

