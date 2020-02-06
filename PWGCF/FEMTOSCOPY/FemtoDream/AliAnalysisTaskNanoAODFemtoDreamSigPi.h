/*
 * AliAnalysisTaskNanoAODFemtoDreamSigPi.cxx
 *
 *  Created on: 30 Jan 2020
 *      Author: M. Korwieser
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOAODFEMTODREAMSIGPI_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOAODFEMTODREAMSIGPI_H_
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

class AliAnalysisTaskNanoAODFemtoDreamSigPi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoAODFemtoDreamSigPi();
  AliAnalysisTaskNanoAODFemtoDreamSigPi(const char *name, bool isMC);
  virtual ~AliAnalysisTaskNanoAODFemtoDreamSigPi();
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
  AliAnalysisTaskNanoAODFemtoDreamSigPi(const AliAnalysisTaskNanoAODFemtoDreamSigPi &);
  AliAnalysisTaskNanoAODFemtoDreamSigPi &operator=(const AliAnalysisTaskNanoAODFemtoDreamSigPi &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
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
  AliVTrack** fGTI;           //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskNanoAODFemtoDreamSigPi,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOAODFEMTODREAMSIGPI_H_ */
