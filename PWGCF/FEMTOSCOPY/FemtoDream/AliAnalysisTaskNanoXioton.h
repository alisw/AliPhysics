/*
 * AliAnalysisTaskNanoXioton.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"

class AliAnalysisTaskNanoXioton : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoXioton();
  AliAnalysisTaskNanoXioton(const char* name);
  virtual ~AliAnalysisTaskNanoXioton();
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

  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!

  ClassDef(AliAnalysisTaskNanoXioton,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOXIOTON_H_ */


