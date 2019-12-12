/*
 * AliAnalysisTaskFemtoTutorial.h
 *
 *  Created on: 11 Dec 2019
 *      Author: Boris Bajtl
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliFemtoDreamTrackCuts.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

#include "AliPIDResponse.h"
#include "AliAODInputHandler.h"


#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"

class AliAnalysisTaskPOmegaPenne : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPOmegaPenne();
  AliAnalysisTaskPOmegaPenne(const char *name, bool isMC);
  virtual ~AliAnalysisTaskPOmegaPenne();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts              *evtCuts){fEventCuts            =   evtCuts;};
  void SetTrackCutsProton(AliFemtoDreamTrackCuts        *trkCuts){fTrackCutsProton      =   trkCuts;};
  void SetTrackCutsAntiProton(AliFemtoDreamTrackCuts    *trkCuts){fTrackCutsAntiProton  =   trkCuts;};
  void SetTrackCutsKaon(AliFemtoDreamTrackCuts          *trkCuts){fTrackCutsKaon        =   trkCuts;};
  void SetTrackCutsAntiKaon(AliFemtoDreamTrackCuts      *trkCuts){fTrackCutsAntiKaon    =   trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig      *config) {fConfig               =   config;};
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                                               //
  TList                             *fOutput;               //!
  AliFemtoDreamEvent                *fEvent;                //!
  AliFemtoDreamTrack                *fTrack;                //!
  AliFemtoDreamEventCuts            *fEventCuts;            //
  AliFemtoDreamTrackCuts            *fTrackCutsProton;      //
  AliFemtoDreamTrackCuts            *fTrackCutsAntiProton;  //
  AliFemtoDreamTrackCuts            *fTrackCutsKaon;        //
  AliFemtoDreamTrackCuts            *fTrackCutsAntiKaon;    //
  AliFemtoDreamCollConfig           *fConfig;               //
  AliFemtoDreamPairCleaner          *fPairCleaner;          //!
  AliFemtoDreamPartCollection       *fPartColl;             //!
  AliAODTrack                       **fGTI;                 //!
  int                               fTrackBufferSize;       //
  ClassDef(AliAnalysisTaskPOmegaPenne,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_ */
