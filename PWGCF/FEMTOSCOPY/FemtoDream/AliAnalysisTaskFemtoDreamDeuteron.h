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
  Float_t GetMass2sq(AliFemtoDreamTrack *track);
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
  void SetTrackCutsDeuteronDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsDeuteronDCA=trkCuts;};
  void SetTrackCutsDeuteronMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsDeuteronMass=trkCuts;};
  void SetTrackCutsAntiDeuteronDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiDeuteronDCA=trkCuts;};
  void SetTrackCutsAntiDeuteronMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiDeuteronMass=trkCuts;};
  void SetTrackCutsProtonDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProtonDCA=trkCuts;};
  void SetTrackCutsProtonMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProtonMass=trkCuts;};
  void SetTrackCutsAntiProtonDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProtonDCA=trkCuts;};
  void SetTrackCutsAntiProtonMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProtonMass=trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;};
 private:
  AliAnalysisTaskFemtoDreamDeuteron(
      const AliAnalysisTaskFemtoDreamDeuteron &task);
  AliAnalysisTaskFemtoDreamDeuteron &operator=(
      const AliAnalysisTaskFemtoDreamDeuteron &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                               //
  TList *fOutput;                           //!  No need to stream this, are compiled
  AliFemtoDreamEvent *fEvent;               //!  on Runtime
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //   Stream these bad boys
  AliFemtoDreamTrackCuts *fTrackCutsDeuteronDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsProtonDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsProtonMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProtonDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProtonMass;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  TH2F* fDRestMass; //!
  TH2F* fAntiDRestMass; //!
  TH2F* fPRestMass; //!
  TH2F* fAntiPRestMass; //! 
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoDreamDeuteron,3)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOTUTORIAL_H_ */
