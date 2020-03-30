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

  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartProton = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartAntiProton = cuts;
  }

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
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetOEventMixing(bool mix) {  fUseOMixing = mix; }
  void SetMCTruth(bool mct) {  fIsMCTruth = mct; }

 private:
  AliAnalysisTaskFemtoDreamPhi(const AliAnalysisTaskFemtoDreamPhi &);
  AliAnalysisTaskFemtoDreamPhi &operator=(const AliAnalysisTaskFemtoDreamPhi &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                            //
  bool fIsMCTruth;                      //
  bool fUseOMixing;                     //
  UInt_t fTrigger;                       //
  TList *fOutput;                        //!
  AliFemtoDreamEvent *fEvent;            //!
  AliFemtoDreamTrack *fTrack;            //!
  AliFemtoDreamv0 *fPhiParticle;         //!
  AliFemtoDreamEventCuts *fEventCuts;    //
  AliFemtoDreamTrackCuts *fPosKaonCuts;  //
  AliFemtoDreamTrackCuts *fNegKaonCuts;  //
  AliFemtoDreamv0Cuts *fPhiCuts;         //
  AliFemtoDreamTrackCuts *fTrackCutsPartProton;      //
  AliFemtoDreamTrackCuts *fTrackCutsPartAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;        //
  AliFemtoDreamPairCleaner *fPairCleaner;  //!
  AliFemtoDreamPartCollection *fPartColl;  //!
  AliAODTrack **fGTI;                      //!
  int fTrackBufferSize;                    //
  ClassDef(AliAnalysisTaskFemtoDreamPhi, 54)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamPhi_H_ */
