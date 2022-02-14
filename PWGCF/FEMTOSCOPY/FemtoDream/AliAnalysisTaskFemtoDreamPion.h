/*
 * AliAnalysisTaskFemtoDreamPion.cxx
 *
 *  Created on: 11 Oct 2019
 *      Author: M. Korwieser
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMPION_H_
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

class AliAnalysisTaskFemtoDreamPion : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamPion();
  AliAnalysisTaskFemtoDreamPion(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoDreamPion();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
  void SetTrackCutsPosPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPosPion=trkCuts;};
  void SetTrackCutsNegPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsNegPion=trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;};
  void SetTrigger(UInt_t trigger) { fTrigger = trigger;};
  void SetIsMC(bool isMC) { fIsMC = isMC;};
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                               //
  TList *fOutput;                           //!
  UInt_t fTrigger;                          //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsPosPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsNegPion;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoDreamPion,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMPION_H_ */
