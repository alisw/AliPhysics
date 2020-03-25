/*
 * AliAnalysisTaskNanoAODFemtoDreamPion.cxx
 *
 *  Created on: 11 Oct 2019
 *      Author: M. Korwieser
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOAODFEMTODREAMPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOAODFEMTODREAMPION_H_
#include "Rtypes.h"

#include "AliAnalysisTaskSE.h"
#include "AliVTrack.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskNanoAODFemtoDreamPion : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoAODFemtoDreamPion();
  AliAnalysisTaskNanoAODFemtoDreamPion(const char *name, bool isMC);
  virtual ~AliAnalysisTaskNanoAODFemtoDreamPion();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
  void SetTrackCutsPosPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPosPion=trkCuts;};
  void SetTrackCutsNegPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsNegPion=trkCuts;};
  //void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton=trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;};
  void SetTrigger(UInt_t trigger) { fTrigger = trigger;};
 private:
  AliAnalysisTaskNanoAODFemtoDreamPion(const AliAnalysisTaskNanoAODFemtoDreamPion &);
  AliAnalysisTaskNanoAODFemtoDreamPion &operator=(const AliAnalysisTaskNanoAODFemtoDreamPion &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;                               //
  TList *fOutput;                           //!
  UInt_t fTrigger;                          //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsPosPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsNegPion;  //
  //AliFemtoDreamTrackCuts *fTrackCutsProton;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliVTrack** fGTI;           //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskNanoAODFemtoDreamPion,2)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOAODFEMTODREAMPION_H_ */
