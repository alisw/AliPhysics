#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamPhi_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamPhi_H_
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

class AliAnalysisTaskFemtoDreamPhi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamPhi();
  AliAnalysisTaskFemtoDreamPhi(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoDreamPhi();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };
  void SetPosKaonCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fPosKaonCuts = trkCuts;
  }
  void SetNegKaonCuts(AliFemtoDreamTrackCuts *trkCuts) {
    fNegKaonCuts = trkCuts;
  }
  void SetPhiCuts(AliFemtoDreamv0Cuts *phiCuts) { fPhiCuts = phiCuts; }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }

 private:
  AliAnalysisTaskFemtoDreamPhi(const AliAnalysisTaskFemtoDreamPhi &);
  AliAnalysisTaskFemtoDreamPhi &operator=(const AliAnalysisTaskFemtoDreamPhi &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                            //
  TList *fOutput;                        //!
  AliFemtoDreamEvent *fEvent;            //!
  AliFemtoDreamTrack *fTrack;            //!
  AliFemtoDreamv0 *fPhiParticle;         //!
  AliFemtoDreamEventCuts *fEventCuts;    //
  AliFemtoDreamTrackCuts *fPosKaonCuts;  //
  AliFemtoDreamTrackCuts *fNegKaonCuts;  //
  AliFemtoDreamv0Cuts *fPhiCuts;         //

  AliFemtoDreamCollConfig *fConfig;        //
  AliFemtoDreamPairCleaner *fPairCleaner;  //!
  AliFemtoDreamPartCollection *fPartColl;  //!
  AliAODTrack **fGTI;                      //!
  int fTrackBufferSize;                    //
  ClassDef(AliAnalysisTaskFemtoDreamPhi, 1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamPhi_H_ */
