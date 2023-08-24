/*
 * AliAnalysisTaskFemtoDreamSigPi.cxx
 *
 *  Created on: 30 Jan 2020
 *      Author: M. Korwieser
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMSIGPI_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMSIGPI_H_
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

class AliAnalysisTaskFemtoDreamSigPi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamSigPi();
  AliAnalysisTaskFemtoDreamSigPi(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoDreamSigPi();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
  void SetTrackCutsSigmaPlus(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsSigmaPlus=trkCuts;};
  void SetTrackCutsNegPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsNegPion=trkCuts;};
  void SetTrackCutsSigmaMinus(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsSigmaMinus=trkCuts;};
  void SetTrackCutsPosPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPosPion=trkCuts;};
  void SetTrackCutsSigma0(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsSigma0=trkCuts;};
  void SetTrackCutsPion0(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPion0=trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;};
  void SetIsMC(bool isMC) { fIsMC = isMC;};
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                               //
  TList *fOutput;                           //!
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsPosPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsNegPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsPion0;  //
  AliFemtoDreamTrackCuts *fTrackCutsSigmaPlus;  //
  AliFemtoDreamTrackCuts *fTrackCutsSigmaMinus;  //
  AliFemtoDreamTrackCuts *fTrackCutsSigma0;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoDreamSigPi,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMSIGPI_H_ */
