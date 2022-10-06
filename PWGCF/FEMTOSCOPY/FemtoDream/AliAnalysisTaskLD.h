/*
 * AliAnalysisTaskLD.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 *  Modifications for LambdaDeuteron
 *      Author: stheckel
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLD_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"

class AliAnalysisTaskLD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLD();
  AliAnalysisTaskLD(const char* name);
  virtual ~AliAnalysisTaskLD();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  Float_t CalculateMassSqTOF(AliFemtoDreamTrack *track);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetPairCleanerUsage(int cleanPairs) {
    fPairCleanerSettings = cleanPairs;
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
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig = config;
  }
 private:
  bool fisLightWeight;//
  int fPairCleanerSettings;// 0=no, 1=ld,ll (default), 2=ld,ll,pl
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
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  AliAODTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskLD,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLD_H_ */

