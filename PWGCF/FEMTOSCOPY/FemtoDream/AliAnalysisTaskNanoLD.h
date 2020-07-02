/*
 * AliAnalysisTaskNanoLD.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 *  Modifications for LambdaDeuteron
 *      Author: stheckel
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLD_H_
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

class AliAnalysisTaskNanoLD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoLD();
  AliAnalysisTaskNanoLD(const char* name);
  virtual ~AliAnalysisTaskNanoLD();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  Float_t CalculateMassSqTOF(AliFemtoDreamTrack *track);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void SetDeuteronCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fDeuteron = trkCuts;
  }
  void SetAntiDeuteronCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiDeuteron = trkCuts;
  }
  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fLambda = v0Cuts;
  }
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fAntiLambda = v0Cuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }
 private:
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  TH1F* fSimpleEventCounter;//!
  TH1F* fSimpleParticleCounter;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fDeuteron;//
  TList* fDeuteronList;//!
  TH2F  *fDeuteronMassSqTOF; //!
  AliFemtoDreamTrackCuts* fAntiDeuteron;//
  TList* fAntiDeuteronList;//!
  TH2F  *fAntiDeuteronMassSqTOF; //!
  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambda;//
  TList* fLambdaList;//!
  AliFemtoDreamv0Cuts* fAntiLambda;//
  TList* fAntiLambdaList;//!
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskNanoLD,2)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLD_H_ */

