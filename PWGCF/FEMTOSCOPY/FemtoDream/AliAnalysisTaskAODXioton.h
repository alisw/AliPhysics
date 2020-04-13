/*
 * AliAnalysisTaskNanoXioton.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKAODXIOTON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKAODXIOTON_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskAODXioton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAODXioton();
  AliAnalysisTaskAODXioton(const char* name, bool isMC);
  //AliAnalysisTaskNanoXioton(const AliAnalysisTaskNanoXioton& analysis) = default;
  //AliAnalysisTaskNanoXioton& operator=(const AliAnalysisTaskNanoXioton& analysis) = default;
  virtual ~AliAnalysisTaskAODXioton();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
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
  AliAnalysisTaskAODXioton(const AliAnalysisTaskAODXioton &task);
  AliAnalysisTaskAODXioton &operator=(const AliAnalysisTaskAODXioton &task);
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  TList* fProtonMCList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  AliFemtoDreamCascade* fCascade;//!
  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;
  TList* fXiMCList;
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  TList* fAntiXiMCList;
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  AliAODTrack **fGTI;  //!

  ClassDef(AliAnalysisTaskAODXioton,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKAODXIOTON_H_ */


