/*
 * AliAnalysisTaskFemtoDreamDeuteron.h
 *
 *  Created on: 21 Mar 2018
 *      Author: bernhardhohlweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMDEUTERON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMDEUTERON_H_
#include "Rtypes.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskFemtoDreamDeuteron : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamDeuteron();
  AliAnalysisTaskFemtoDreamDeuteron(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoDreamDeuteron();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
  void SetTrackCutsPart1(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPart1=trkCuts;};
  void SetTrackCutsPart2(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPart2=trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;};
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                               //
  TList *fOutput;                           //!  No need to stream this, are compiled
  AliFemtoDreamEvent *fEvent;               //!  on Runtime
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //   Stream these bad boys
  AliFemtoDreamTrackCuts *fTrackCutsPart1;  //
  AliFemtoDreamTrackCuts *fTrackCutsPart2;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoDreamDeuteron,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOTUTORIAL_H_ */
