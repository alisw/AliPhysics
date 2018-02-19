/*
 * AliFemtoDreamAnalysis.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIFEMTODREAMANALYSIS_H_
#define ALIFEMTODREAMANALYSIS_H_

#include "AliAODEvent.h"
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
#include "Rtypes.h"
class AliFemtoDreamAnalysis {
 public:
  AliFemtoDreamAnalysis();
  void SetMVPileUp(bool mvPileUp){fMVPileUp=mvPileUp;};
  void SetEvtCutQA(bool setQA){fEvtCutQA=setQA;};
  void SetEventCuts(AliFemtoDreamEventCuts *cuts){fEvtCuts=cuts;};
  TList *GetEventCutHists(){return fEvtCuts->GetHistList();};
  void SetTrackCuts(AliFemtoDreamTrackCuts *cuts){fTrackCuts=cuts;};
  TList *GetTrackCutHists(){return fTrackCuts->GetQAHists();};
  TList *GetTrackCutHistsMC(){return fTrackCuts->GetMCQAHists();};
  void SetAntiTrackCuts(AliFemtoDreamTrackCuts *cuts){fAntiTrackCuts=cuts;};
  TList *GetAntitrackCutHists(){return fAntiTrackCuts->GetQAHists();};
  TList *GetAntitrackCutHistsMC(){return fAntiTrackCuts->GetMCQAHists();};
  void Setv0Cuts(AliFemtoDreamv0Cuts *cuts){fv0Cuts=cuts;};
  TList *Getv0CutHist(){return fv0Cuts->GetQAHists();};
  TList *Getv0MCHist(){return fv0Cuts->GetMCQAHists();};
//  TList *Getv0PosDaugHist(){return fv0Cuts->GetQAHistsPosDaug();};
//  TList *Getv0PosDaugMCHist(){return fv0Cuts->GetMCQAHistsPosDaug();};
//  TList *Getv0NegDaugHist(){return fv0Cuts->GetQAHistsNegDaug();};
//  TList *Getv0NegDaugMCHist(){return fv0Cuts->GetMCQAHistsNegDaug();};
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts *cuts){fAntiv0Cuts=cuts;};
  TList *GetAntiv0CutHist(){return fAntiv0Cuts->GetQAHists();};
  TList *GetAntiv0MCHist(){return fAntiv0Cuts->GetMCQAHists();};
//  TList *GetAntiv0PosDaugHist(){return fAntiv0Cuts->GetQAHistsPosDaug();};
//  TList *GetAntiv0PosDaugMCHist(){return fAntiv0Cuts->GetMCQAHistsPosDaug();};
//  TList *GetAntiv0NegDaugHist(){return fAntiv0Cuts->GetQAHistsNegDaug();};
//  TList *GetAntiv0NegDaugMCHist(){return fAntiv0Cuts->GetMCQAHistsNegDaug();};
  void SetCascadeCuts(AliFemtoDreamCascadeCuts *cuts){fCascCuts=cuts;};
  TList *GetCascadeCutHist(){return fCascCuts->GetQAHists();};
  void SetAntiCascadeCuts(AliFemtoDreamCascadeCuts *cuts){fAntiCascCuts=cuts;};
  TList *GetAntiCascadeCutHist(){return fAntiCascCuts->GetQAHists();};
  TList *GetResultList() {return fPartColl->GetHistList();};
  TList *GetResultQAList() {return fPartColl->GetQAList();};
  TList *GetQAList() {return fQA;};
  void SetTrackBufferSize(int size){fTrackBufferSize=size;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *conf) {fConfig=conf;};
  void Init(bool isMonteCarlo);
  TString ClassName() {return "AliFemtoDreamAnalysis";};
  void Make(AliAODEvent *evt);
  virtual ~AliFemtoDreamAnalysis();
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
  int fTrackBufferSize;
  AliAODTrack           **fGTI;             //!
  AliFemtoDreamCollConfig *fConfig;         //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  ClassDef(AliFemtoDreamAnalysis,1)
};

#endif /* ALIFEMTODREAMANALYSIS_H_ */
