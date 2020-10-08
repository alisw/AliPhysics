#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALINANALYSISTASKNANOPPCOALESCENCE_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALINANALYSISTASKNANOPPCOALESCENCE_H_

#include "AliAnalysisTaskSE.h"
#include "AliVTrack.h"
#include "AliMCEvent.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMultSelection.h"
#include "AliNanoAODTrack.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include "TChain.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoPPCoalescence : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoPPCoalescence();
  AliAnalysisTaskNanoPPCoalescence(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskNanoPPCoalescence();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetEventCuts(AliFemtoDreamEventCuts *cuts) {
    fEvtCuts = cuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fProtonTrack = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiProtonTrack = cuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }
  void SetMCTruth(bool mct) {
    fIsMCTruth = mct;
  }

 private:
  AliAnalysisTaskNanoPPCoalescence(
    const AliAnalysisTaskNanoPPCoalescence &task);
  AliAnalysisTaskNanoPPCoalescence &operator=(
    const AliAnalysisTaskNanoPPCoalescence &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;//
  bool fIsMCTruth;//

  AliVEvent *fInputEvent;//!
  AliFemtoDreamEvent *fEvent;//!
  AliFemtoDreamEventCuts *fEvtCuts;//
  AliFemtoDreamTrack *fTrack;//
  AliFemtoDreamTrackCuts *fProtonTrack;//
  AliFemtoDreamTrackCuts *fAntiProtonTrack;//


  int fTrackBufferSize; //
  AliVTrack **fGTI;  //!

  TList *fEvtList;//!
  TList *fProtonList;//!
  TList *fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList *fAntiProtonMCList;//!

  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!

  TList *fResults;                          //!
  TList *fResultsQA;                        //!

  ClassDef(AliAnalysisTaskNanoPPCoalescence, 1)
};
#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALINANALYSISTASKNANOPPCOALESCENCE_H_ */
