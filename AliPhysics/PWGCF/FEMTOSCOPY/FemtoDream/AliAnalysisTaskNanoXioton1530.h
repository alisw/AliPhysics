/*
 * AliAnalysisTaskNanoXioton.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON1530_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON1530_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskNanoXioton1530 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoXioton1530();
  AliAnalysisTaskNanoXioton1530(const char* name);
  //AliAnalysisTaskNanoXioton(const AliAnalysisTaskNanoXioton& analysis) = default;
  //AliAnalysisTaskNanoXioton& operator=(const AliAnalysisTaskNanoXioton& analysis) = default;
  virtual ~AliAnalysisTaskNanoXioton1530();
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
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
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
  void SetXi1530Cuts(AliFemtoDreamv0Cuts* cascCuts) {
    fXi1530Cuts = cascCuts;
  }
  void SetAntiXi1530Cuts(AliFemtoDreamv0Cuts* cascCuts) {
    fAntiXi1530Cuts = cascCuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }
 private:
  AliAnalysisTaskNanoXioton1530(const AliAnalysisTaskNanoXioton1530 &task);
  AliAnalysisTaskNanoXioton1530 &operator=(const AliAnalysisTaskNanoXioton1530 &task);
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  AliFemtoDreamTrackCuts* fPion;//
  TList* fPionList;//!
  AliFemtoDreamTrackCuts* fAntiPion;//
  TList* fAntiPionList;//!
  AliFemtoDreamCascade* fCascade;//!
  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliFemtoDreamv0 *fXi1530;         //!
  AliFemtoDreamv0Cuts *fXi1530Cuts;//
  TList* fXi1530List;
  AliFemtoDreamv0Cuts *fAntiXi1530Cuts;//
  TList* fAntiXi1530List;
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!

  ClassDef(AliAnalysisTaskNanoXioton1530,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON1530_H_ */


