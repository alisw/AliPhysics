/*
 * AliAnalysisTaskNanoBBar.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOBBAR_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOBBAR_H_
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


class AliAnalysisTaskNanoBBar : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoBBar();
  AliAnalysisTaskNanoBBar(const char* name, bool isMC);
  virtual ~AliAnalysisTaskNanoBBar();
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
  void SetUseDumpster(bool use) {
    fUseDumpster = use;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
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
  AliAnalysisTaskNanoBBar(const AliAnalysisTaskNanoBBar &task);
  AliAnalysisTaskNanoBBar &operator=(const AliAnalysisTaskNanoBBar &task);
  bool fisLightWeight;//
  bool fIsMC;        //
  bool fUseDumpster;  //
  TList *fQA;        //!
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
  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambda;//
  TList* fLambdaList;//!
  TList* fLambdaMCList;//!
  AliFemtoDreamv0Cuts* fAntiLambda;//
  TList* fAntiLambdaList;//!
  TList* fAntiLambdaMCList;//!
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
  AliFemtoDreamDump *fProtonAntiProtonDump; //!
  AliFemtoDreamDump *fProtonAntiLambdaDump; //!
  AliFemtoDreamDump *fAntiProtonLambdaDump; //!
  AliFemtoDreamDump *fLambdaAntiLambdaDump; //!
  AliFemtoDreamDump *fProtonAntiXiDump; //!
  AliFemtoDreamDump *fAntiProtonXiDump; //!
  AliFemtoDreamDump *fLambdaAntiXiDump; //!
  AliFemtoDreamDump *fAntiLambdaXiDump; //!
  AliFemtoDreamDump *fXiAntiXiDump; //!
  TList* fDumpster; //!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskNanoBBar,4)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOBBAR_H_ */


