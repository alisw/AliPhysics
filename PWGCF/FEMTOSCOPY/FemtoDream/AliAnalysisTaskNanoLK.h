/*
 * AliAnalysisTaskNanoLK.h
 *
 *  Created on: Nov 12, 2021
 *      Author: Valentina Mantovani Sarti
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLK_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLK_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"

class AliAnalysisTaskNanoLK : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoLK();
  AliAnalysisTaskNanoLK(const char* name, bool isMC);
  virtual ~AliAnalysisTaskNanoLK();
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

  void SetPosKaonCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fPosKaonCuts = trkCuts;
  }

  void SetNegKaonCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fNegKaonCuts = trkCuts;
  }

  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fLambdaCuts = v0Cuts;
  }
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fAntiLambdaCuts = v0Cuts;
  }

  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }

 private:
  AliAnalysisTaskNanoLK(const AliAnalysisTaskNanoLK &task);
  AliAnalysisTaskNanoLK &operator=(const AliAnalysisTaskNanoLK &task);
  bool fisLightWeight;//
  bool fIsMC;        //
  TList *fQA;        //!
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fPosKaonCuts;//
  TList* fPosKaonList;//!
  TList* fPosKaonMCList;//!
  AliFemtoDreamTrackCuts* fNegKaonCuts;//
  TList* fNegKaonList;//!
  TList* fNegKaonMCList;//!
  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambdaCuts;//
  TList* fLambdaList;//!
  TList* fLambdaMCList;//!
  AliFemtoDreamv0Cuts* fAntiLambdaCuts;//
  TList* fAntiLambdaList;//!
  TList* fAntiLambdaMCList;//!
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
  ClassDef(AliAnalysisTaskNanoLK,4)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLK_H_ */
