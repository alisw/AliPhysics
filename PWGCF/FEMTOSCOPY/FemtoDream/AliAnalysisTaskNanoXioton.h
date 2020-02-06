/*
 * AliAnalysisTaskNanoXioton.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON_H_
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
#include "AliFemtoDreamBaseDump.h"

class AliAnalysisTaskNanoXioton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoXioton();
  AliAnalysisTaskNanoXioton(const char* name, bool isMC);
  //AliAnalysisTaskNanoXioton(const AliAnalysisTaskNanoXioton& analysis) = default;
  //AliAnalysisTaskNanoXioton& operator=(const AliAnalysisTaskNanoXioton& analysis) = default;
  virtual ~AliAnalysisTaskNanoXioton();
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
  AliAnalysisTaskNanoXioton(const AliAnalysisTaskNanoXioton &task);
  AliAnalysisTaskNanoXioton &operator=(const AliAnalysisTaskNanoXioton &task);
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
  AliFemtoDreamDump *fProtonProtonDump; //!
  AliFemtoDreamDump *fAntiProtonAntiProtonDump; //!
  TList* fDumpster; //!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!

  ClassDef(AliAnalysisTaskNanoXioton,2)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON_H_ */


