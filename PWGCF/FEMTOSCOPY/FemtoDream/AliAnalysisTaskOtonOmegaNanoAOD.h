/*
 * AliAnalysisTaskOtonOmegaNanoAOD.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONOMEGANANOAOD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONOMEGANANOAOD_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliOtonOmegaCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"

class AliAnalysisTaskOtonOmegaNanoAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskOtonOmegaNanoAOD();
  AliAnalysisTaskOtonOmegaNanoAOD(const char* name);
  //AliAnalysisTaskOtonOmegaNanoAOD(const AliAnalysisTaskOtonOmegaNanoAOD& analysis) = default;
  //AliAnalysisTaskOtonOmegaNanoAOD& operator=(const AliAnalysisTaskOtonOmegaNanoAOD& analysis) = default;
  virtual ~AliAnalysisTaskOtonOmegaNanoAOD();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
  }

  void SetXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fXi = cascCuts;
  }
  void SetAntiXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiXi = cascCuts;
  }
  void SetOmegaCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fOmega = cascCuts;
  }
  void SetAntiOmegaCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiOmega = cascCuts;
  }
/*
  void SetXiCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fXi = cascCuts;
  }
  void SetAntiXiCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fAntiXi = cascCuts;
  }
  void SetOmegaCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fOmega = cascCuts;
  }
  void SetAntiOmegaCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fAntiOmega = cascCuts;
  }
*/
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }
 private:
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  AliFemtoDreamCascade* fCascade;//!

  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliFemtoDreamCascadeCuts* fOmega;//
  TList* fOmegaList;
  AliFemtoDreamCascadeCuts* fAntiOmega;//
  TList* fAntiOmegaList;

/*
  AliOtonOmegaCascadeCuts* fXi;//
  TList* fXiList;
  AliOtonOmegaCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliOtonOmegaCascadeCuts* fOmega;//
  TList* fOmegaList;
  AliOtonOmegaCascadeCuts* fAntiOmega;//
  TList* fAntiOmegaList;
*/

  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!

  ClassDef(AliAnalysisTaskOtonOmegaNanoAOD,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskOtonOmegaNanoAOD_H_ */


