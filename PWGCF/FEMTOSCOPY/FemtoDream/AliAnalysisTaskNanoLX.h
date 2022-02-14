/*
 * AliAnalysisTaskNanoLX.h
 *
 *  Created on: June 26, 2019
 *      Author: me
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLX_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLX_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskNanoLX : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoLX();
  AliAnalysisTaskNanoLX(const char* name);
  //AliAnalysisTaskNanoXioton(const AliAnalysisTaskNanoXioton& analysis) = default;
  //AliAnalysisTaskNanoXioton& operator=(const AliAnalysisTaskNanoXioton& analysis) = default;
  virtual ~AliAnalysisTaskNanoLX();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();                  //concerns protons only (?)
  void StoreGlobalTrackReference(AliVTrack *track);
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fLambda = v0Cuts;
  }
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fAntiLambda = v0Cuts;
  }
  void SetXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fXi = cascCuts;
  }
  void SetAntiXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiXi = cascCuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }
 private:
  AliAnalysisTaskNanoLX(const AliAnalysisTaskNanoLX &task);
  AliAnalysisTaskNanoLX &operator=(const AliAnalysisTaskNanoLX &task);
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambda;//
  TList* fLambdaList;
  AliFemtoDreamv0Cuts* fAntiLambda;//
  TList* fAntiLambdaList;
  AliFemtoDreamCascade* fCascade;//!
  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!

  ClassDef(AliAnalysisTaskNanoLX,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLX_H_ */


