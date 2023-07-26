/*
 * AliAnalysisTaskFemtoDreamRho.cxx
 *
 *  Created on: 19 Jul 2023
 *      Author: M. Korwieser
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamRho_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamRho_H_
#include "Rtypes.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"

class AliAnalysisTaskFemtoDreamRho : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamRho();
  AliAnalysisTaskFemtoDreamRho(const char *name, bool isMC, bool doCleaning);
  virtual ~AliAnalysisTaskFemtoDreamRho();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };
  void SetPosPionCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fPosPionCuts = trkCuts;
  }
  void SetNegPionCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fNegPionCuts = trkCuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fPosProtonCuts = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fNegProtonCuts = trkCuts;
  }
  void SetRhoCuts(AliFemtoDreamv0Cuts *rhoCuts) { 
    fRhoCuts = rhoCuts; 
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger;};
  void SetIsMC(bool isMC) { fIsMC = isMC;};
  void SetDoCleaning(bool doCleaning) { fDoCleaning = doCleaning;};
 private:
  AliAnalysisTaskFemtoDreamRho(const AliAnalysisTaskFemtoDreamRho &);
  AliAnalysisTaskFemtoDreamRho &operator=(const AliAnalysisTaskFemtoDreamRho &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                               //
  TList *fOutput;                           //!
  UInt_t fTrigger;                          //
  bool fDoCleaning;                         //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamv0 *fRhoParticle;            //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fPosPionCuts;     //
  AliFemtoDreamTrackCuts *fNegPionCuts;     //
  AliFemtoDreamTrackCuts *fPosProtonCuts;   //
  AliFemtoDreamTrackCuts *fNegProtonCuts;   //
  AliFemtoDreamv0Cuts *fRhoCuts;            //

  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack **fGTI;                       //!
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoDreamRho, 1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamRho_H_ */
