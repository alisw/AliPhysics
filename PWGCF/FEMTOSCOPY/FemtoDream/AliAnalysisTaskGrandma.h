/*
 * AliAnalysisTaskGrandma.h
 *
 *  Created on: Sep 11, 2018
 *      Author: hohlweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGRANDMA_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGRANDMA_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamAnalysis.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamCollConfig.h"
#include "TList.h"
class AliAnalysisTaskGrandma : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGrandma();
  AliAnalysisTaskGrandma(const char *name, bool isMC);
  virtual ~AliAnalysisTaskGrandma();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {
  }  //
  ;
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
  void SetTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCuts = cuts;
  }
  ;
  void SetAntiTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiTrackCuts = cuts;
  }
  ;
  void Setv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fv0Cuts = cuts;
  }
  ;
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts *cuts){
    fAntiv0Cuts = cuts;
  }
  ;
  void SetXiCuts(AliFemtoDreamCascadeCuts* cuts) {
    fXi = cuts;
  }
  ;
  void SetAntiXiCuts(AliFemtoDreamCascadeCuts* cuts) {
    fAntiXi = cuts;
  }
  ;
  void SetTrackBufferSize(int trackBuffer) {
    fTrackBufferSize = trackBuffer;
  }
  ;
  void StoreGlobalTrackReference(AliAODTrack *track);
  void ResetGlobalTrackReference();
  void SetCollectionConfig(AliFemtoDreamCollConfig* conf) {
    fConfig = conf;
  }
  ;
 private:
  AliAnalysisTaskGrandma(const AliAnalysisTaskGrandma &task);
  AliAnalysisTaskGrandma &operator=(const AliAnalysisTaskGrandma &task);
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  TList *fQA;                               //!
  bool fMinBookingME;                       //
  bool fMinBookingSample;                   //
  bool fMVPileUp;                           //
  bool fEvtCutQA;                           //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamEventCuts *fEvtCuts;         // stream this one!
  TList *fEvtHistList;                      //!
  AliFemtoDreamTrack *fFemtoTrack;          //!
  AliFemtoDreamTrackCuts *fTrackCuts;       //
  TList *fTrackCutHistList;                 //!
  TList *fTrackCutHistMCList;               //!
  AliFemtoDreamTrackCuts *fAntiTrackCuts;   //
  TList *fAntiTrackCutHistList;             //!
  TList *fAntiTrackCutHistMCList;           //!
  AliFemtoDreamv0 *fFemtov0;                //!
  AliFemtoDreamv0Cuts *fv0Cuts;             //
  TList *fv0CutHistList;                 //!
  TList *fv0CutHistMCList;               //!
  AliFemtoDreamv0Cuts *fAntiv0Cuts;         //
  TList *fAntiv0CutHistList;             //!
  TList *fAntiv0CutHistMCList;           //!
  AliFemtoDreamCascade* fCascade;//!
  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;//!
  TList* fXiMCList;//!
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;//!
  TList* fAntiXiMCList;//!
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliFemtoDreamCollConfig *fConfig;         //
  TList *fResultList;                       //!
  TList *fResultQAList;                     //!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  AliAODTrack **fGTI;                       //!
ClassDef(AliAnalysisTaskGrandma,4)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGRANDMA_H_ */
