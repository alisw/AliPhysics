/*
 * AliAnalysisTaskFemtoSMI.h
 * 
 *   Created on: 2020-08-17
 *   Author: S. Schneider
 *
 * This file is based on a copy of:
 *   AliAnalysisTaskFemtoTutorial.h
 * with some parts taken from:
 *   AliAnalysisTaskAODLoton.h
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOSMI_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOSMI_H_
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

// new additions following AliAnalysisTaskAODLoton.h:
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamPartCollection.h" // necessary here?
#include "AliFemtoDreamControlSample.h" // necessary here?

class AliAnalysisTaskFemtoSMI : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoSMI();
  AliAnalysisTaskFemtoSMI(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoSMI();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){}; //NIC
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts=evtCuts;};
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {fProton = trkCuts;};
  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {fLambda = v0Cuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig=config;}; 
  
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }

 private:
  void ResetGlobalTrackReference(); //NIC
  void StoreGlobalTrackReference(AliAODTrack *track); //NIC
  bool fIsMC;                               // NIC, but have fisLightWeight
  TList *fOutput;                           //!  No need to stream this, are compiled
  AliFemtoDreamEvent *fEvent;               //!  on Runtime
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //   Stream these bad boys
  AliFemtoDreamTrackCuts *fProton;			//
  AliFemtoDreamv0 *fv0;						//! 		
  AliFemtoDreamv0Cuts *fLambda;				//

  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoSMI,1)

  UInt_t fTrigger;							//

};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOSMI_H_ */
