/*
 * AliFemtoDreamAnalysis.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIFEMTODREAMANALYSIS_H_
#define ALIFEMTODREAMANALYSIS_H_

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamControlSample.h"
#include "Rtypes.h"
class AliFemtoDreamAnalysis {
 public:
  AliFemtoDreamAnalysis();
  AliFemtoDreamAnalysis(const AliFemtoDreamAnalysis& analysis);
  AliFemtoDreamAnalysis& operator=(const AliFemtoDreamAnalysis& analysis);
  virtual ~AliFemtoDreamAnalysis();
  void SetMVPileUp(bool mvPileUp) {
    fMVPileUp = mvPileUp;
  }
  ;
  void SetEvtCutQA(bool setQA) {
    fEvtCutQA = setQA;
  }
  ;
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) {
    fEvtCuts = cuts;
  }
  ;
  TList *GetEventCutHists() {
    return fEvtCuts->GetHistList();
  }
  ;
  void SetTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCuts = cuts;
  }
  ;
  TList *GetTrackCutHists() {
    return fTrackCuts->GetQAHists();
  }
  ;
  TList *GetTrackCutHistsMC() {
    return fTrackCuts->GetMCQAHists();
  }
  ;
  void SetAntiTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiTrackCuts = cuts;
  }
  ;
  TList *GetAntitrackCutHists() {
    return fAntiTrackCuts->GetQAHists();
  }
  ;
  TList *GetAntitrackCutHistsMC() {
    return fAntiTrackCuts->GetMCQAHists();
  }
  ;
  void Setv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fv0Cuts = cuts;
  }
  ;
  TList *Getv0CutHist() {
    return fv0Cuts->GetQAHists();
  }
  ;
  TList *Getv0MCHist() {
    return fv0Cuts->GetMCQAHists();
  }
  ;
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fAntiv0Cuts = cuts;
  }
  ;
  TList *GetAntiv0CutHist() {
    return fAntiv0Cuts->GetQAHists();
  }
  ;
  TList *GetAntiv0MCHist() {
    return fAntiv0Cuts->GetMCQAHists();
  }
  ;
  void SetCascadeCuts(AliFemtoDreamCascadeCuts *cuts) {
    fCascCuts = cuts;
  }
  ;
  TList *GetCascadeCutHist() {
    return fCascCuts->GetQAHists();
  }
  ;
  TList *GetCascadeMCHist() {
    return fCascCuts->GetMCQAHists();
  }
  ;
  void SetAntiCascadeCuts(AliFemtoDreamCascadeCuts *cuts) {
    fAntiCascCuts = cuts;
  }
  ;
  TList *GetAntiCascadeMCHist() {
    return fAntiCascCuts->GetMCQAHists();
  }
  ;
  TList *GetAntiCascadeCutHist() {
    return fAntiCascCuts->GetQAHists();
  }
  ;
  TList *GetResultList() {
    return fPartColl ? fPartColl->GetHistList() : nullptr;
  }
  ;
  TList *GetResultQAList() {
    return fPartColl ? fPartColl->GetQAList() : nullptr;
  }
  ;
  TList *GetResultSampleList() {
    return fControlSample ? fControlSample->GetHistList() : nullptr;
  }
  ;
  TList *GetResultSampleQAList() {
    return fControlSample ? fControlSample->GetQAList() : nullptr;
  }
  ;
  TList *GetQAList() {
    return fQA;
  }
  ;
  void SetTrackBufferSize(int size) {
    fTrackBufferSize = size;
  }
  ;
  void SetCollectionConfig(AliFemtoDreamCollConfig *conf) {
    fConfig = conf;
  }
  ;
  void Init(bool isMonteCarlo, UInt_t trigger);
  TString ClassName() {
    return "AliFemtoDreamAnalysis";
  }
  ;
  void Make(AliAODEvent *evt);
  void Make(AliESDEvent *evt, AliMCEvent *mcEvent);
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fMVPileUp;                           //!
  bool fEvtCutQA;                           //!
  TList *fQA;                               //!
  AliFemtoDreamTrack *fFemtoTrack;          //!
  AliFemtoDreamv0 *fFemtov0;                //!
  AliFemtoDreamCascade *fFemtoCasc;         //!
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamEventCuts *fEvtCuts;         //!
  AliFemtoDreamTrackCuts *fTrackCuts;       //!
  AliFemtoDreamTrackCuts *fAntiTrackCuts;   //!
  AliFemtoDreamv0Cuts *fv0Cuts;             //!
  AliFemtoDreamv0Cuts *fAntiv0Cuts;         //!
  AliFemtoDreamCascadeCuts *fCascCuts;      //!
  AliFemtoDreamCascadeCuts *fAntiCascCuts;  //!
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamControlSample *fControlSample;   //!
  int fTrackBufferSize;
  AliAODTrack **fGTI;             //!
  AliFemtoDreamCollConfig *fConfig;         //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  bool fIsMC;                               //!
ClassDef(AliFemtoDreamAnalysis,3)
};

#endif /* ALIFEMTODREAMANALYSIS_H_ */
