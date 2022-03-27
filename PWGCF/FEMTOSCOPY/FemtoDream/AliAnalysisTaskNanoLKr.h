/*
 * AliAnalysisTaskNanoLKr.h
 *
 *  Created on: November 21, 2021
 *      Author: rossana
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLKR_H_ ///these 2 rows are include guards, to avoid double inclusion
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLKR_H_
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

///AliAnalysisTaskNanoLKr is our class, 
class AliAnalysisTaskNanoLKr : public AliAnalysisTaskSE { //it inherits AliAnalysisTaskSE members
 public: 
  AliAnalysisTaskNanoLKr();  ///constructor 1. Not parameters 
  AliAnalysisTaskNanoLKr(const char* name, bool isMC); ///constructor 2
  virtual ~AliAnalysisTaskNanoLKr(); ///destructor

  ///prototypes of the functions:
  virtual void UserCreateOutputObjects(); ///define whataver output we want to our root files
  virtual void UserExec(Option_t *option); ///function called for each single event over which the analysis runs: is the "event loop". It checks if an event has properties you want to have in the analysis (if not, return); it access the physics objects specific to your analysis; it fills the histograms etc

  ///define a series of void functions
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);


  ///function for the output histograms. Useful. For fast runs
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }


  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }

  void SetPosKaonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fPosKaon = trkCuts;
  }
  void SetNegKaonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fNegKaon = trkCuts;
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

  void SetTrigger(UInt_t trigger) { fTrigger = trigger; } ///ask

 private:
 ///two standard rows, don t modify them
  AliAnalysisTaskNanoLKr(const AliAnalysisTaskNanoLKr &task);
  AliAnalysisTaskNanoLKr &operator=(const AliAnalysisTaskNanoLKr &task);

  ///define every kind of variables. some were already in the void functions
  bool fisLightWeight;//
  bool fIsMC;        //  ///switch Montecarlo on or off
  UInt_t fTrigger; //
  TList *fQA;        //!
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!

  AliFemtoDreamTrack* fTrack;//!

  AliFemtoDreamTrackCuts* fPosKaon;//
  TList* fPosKaonList;//!
  TList* fPosKaonMCList;//!
  AliFemtoDreamTrackCuts* fNegKaon;//
  TList* fNegKaonList;//!
  TList* fNegKaonMCList;//!

  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambda;//
  TList* fLambdaList;//!
  TList* fLambdaMCList;//!
  AliFemtoDreamv0Cuts* fAntiLambda;//
  TList* fAntiLambdaList;//!
  TList* fAntiLambdaMCList;//!

  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!

  TList *fResults;//!
  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  ///TList *fResultsSample;//!
  ///TList *fResultsSampleQA;//!
  int fTrackBufferSize;//   ///vector to fill with particles, to create pairs
  AliVTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskNanoLKr,4)
};

#endif 
///third include guard, to avoid double inclusion
