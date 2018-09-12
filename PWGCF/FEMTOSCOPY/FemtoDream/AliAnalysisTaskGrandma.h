/*
 * AliAnalysisTaskGrandma.h
 *
 *  Created on: Sep 11, 2018
 *      Author: hohlweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGRANDMA_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGRANDMA_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "TList.h"
class AliAnalysisTaskGrandma : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGrandma();
  AliAnalysisTaskGrandma(const char *name, bool isMC);
  virtual ~AliAnalysisTaskGrandma();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {
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
  //sadfsdaf
 private:
  bool  fIsMC;                              //
  TList *fQA;                               //!
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
ClassDef(AliAnalysisTaskGrandma,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGRANDMA_H_ */
