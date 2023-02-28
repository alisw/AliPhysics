/*
 * AliAnalysisTaskNanoBBar.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIPI_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIPI_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"
#include "AliFemtoDreamBaseDump.h"


class AliAnalysisTaskNanoXiPi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoXiPi();
  AliAnalysisTaskNanoXiPi(const char* name, bool isMC);
  virtual ~AliAnalysisTaskNanoXiPi();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void SetPionCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fPion = trkCuts;
  }
  void SetAntiPionCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiPion = trkCuts;
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
  AliAnalysisTaskNanoXiPi(const AliAnalysisTaskNanoXiPi &task);
  AliAnalysisTaskNanoXiPi &operator=(const AliAnalysisTaskNanoXiPi &task);
  bool fisLightWeight;//
  bool fIsMC;        //
  TList *fQA;        //!
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fPion;//
  TList* fPionList;//!
  TList* fPionMCList;//!
  AliFemtoDreamTrackCuts* fAntiPion;//
  TList* fAntiPionList;//!
  TList* fAntiPionMCList;//!
  AliFemtoDreamCascade* fCascade;//!
  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;//!
  TList* fXiMCList;//!
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;//!
  TList* fAntiXiMCList;//!
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskNanoXiPi,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIPI_H_ */


