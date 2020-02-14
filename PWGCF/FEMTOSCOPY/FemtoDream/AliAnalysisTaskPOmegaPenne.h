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
#include "AliLog.h"
#include "AliAODTrack.h"
#include "TROOT.h"
#include "TSystem.h"
#include "AliPIDResponse.h"
#include "AliAODInputHandler.h"

#include "AliFemtoDreamBasePart.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamCascadeCuts.h"

class AliAnalysisTaskPOmegaPenne : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPOmegaPenne();
  AliAnalysisTaskPOmegaPenne(const char *name, bool isMC);
  virtual ~AliAnalysisTaskPOmegaPenne();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(            AliFemtoDreamEventCuts   *evtCuts )  { fEventCuts            =   evtCuts;  };
  void SetTrackCutsProton(      AliFemtoDreamTrackCuts   *trkCuts )  { fTrackCutsProton      =   trkCuts;  };
  void SetTrackCutsAntiProton(  AliFemtoDreamTrackCuts   *trkCuts )  { fTrackCutsAntiProton  =   trkCuts;  };
  void SetTrackCutsXion(        AliFemtoDreamCascadeCuts *cascCuts)  { fCascadeCutsXion      =   cascCuts; };
  void SetTrackCutsAntiXion(    AliFemtoDreamCascadeCuts *cascCuts)  { fCascadeCutsAntiXion  =   cascCuts; };
  void SetCollectionConfig(     AliFemtoDreamCollConfig  *config  )  { fConfig               =   config;   };
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool                                fIsMC;                 //
  AliAODEvent                        *Event;                 //      UserExec:Current Event
  AliAODTrack                        *track;                 //      UserExec:Current Track
  TList                              *fOutput;               //!
  AliFemtoDreamEvent                 *fEvent;                //!
  AliFemtoDreamTrack                 *fTrack;                //!
  AliFemtoDreamCascade               *fCascade;              //!
  AliFemtoDreamEventCuts             *fEventCuts;            //
  AliFemtoDreamTrackCuts             *fTrackCutsProton;      //
  AliFemtoDreamTrackCuts             *fTrackCutsAntiProton;  //
  AliFemtoDreamCascadeCuts           *fCascadeCutsXion;      //
  AliFemtoDreamCascadeCuts           *fCascadeCutsAntiXion;  //
  AliFemtoDreamCollConfig            *fConfig;               //
  AliFemtoDreamPairCleaner           *fPairCleaner;          //!
  AliFemtoDreamPartCollection        *fPartColl;             //!
  AliAODTrack                       **fGTI;                  //!
  int                                 fTrackBufferSize;      //
  ClassDef(AliAnalysisTaskPOmegaPenne,5)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_POMEGA_PENNE_H_ */
