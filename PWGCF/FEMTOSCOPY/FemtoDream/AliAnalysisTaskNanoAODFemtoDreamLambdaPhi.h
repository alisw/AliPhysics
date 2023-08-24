#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskNanoAODFemtoDreamLambdaPhi_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskNanoAODFemtoDreamLambdaPhi_H_
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
#include "AliFemtoDreamControlSample.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoAODFemtoDreamLambdaPhi : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskNanoAODFemtoDreamLambdaPhi();
  AliAnalysisTaskNanoAODFemtoDreamLambdaPhi(const char *name, bool isMC);
  virtual ~AliAnalysisTaskNanoAODFemtoDreamLambdaPhi();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };

  void SetLambdaCuts(AliFemtoDreamv0Cuts *cuts)
  {
    fLambdaCuts = cuts;
  }
  void SetAntiLambdaCuts(AliFemtoDreamv0Cuts *cuts)
  {
    fAntiLambdaCuts = cuts;
  }
  void SetPosKaonCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fPosKaonCuts = trkCuts;
  }
  void SetNegKaonCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fNegKaonCuts = trkCuts;
  }
  void SetPhiCuts(AliFemtoDreamv0Cuts *phiCuts) { fPhiCuts = phiCuts; }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config)
  {
    fConfig = config;
  }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }

private:
  AliAnalysisTaskNanoAODFemtoDreamLambdaPhi(const AliAnalysisTaskNanoAODFemtoDreamLambdaPhi &);
  AliAnalysisTaskNanoAODFemtoDreamLambdaPhi &operator=(const AliAnalysisTaskNanoAODFemtoDreamLambdaPhi &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;               //
  bool fUseOMixing;         //
  float fInvMassCutSBdown;  //
  float fInvMassCutSBup;    //
  UInt_t fTrigger;          //
  TList *fResults;          //!
  TList *fResultsQA;        //!
  TList *fQA;               //!
  TList *fEvtList;          //!
  TList *fLambdaList;       //!
  TList *fLambdaMCList;     //!
  TList *fAntiLambdaList;   //!
  TList *fAntiLambdaMCList; //!
  TList *fKaonPlusList;     //!
  TList *fKaonPlusMCList;   //!
  TList *fKaonMinusList;    //!
  TList *fKaonMinusMCList;  //!
  TList *fPhiList;          //!
  TList *fPhiMCList;        //!
  AliVEvent *fInputEvent;                 //! current event
  AliFemtoDreamEvent *fEvent;             //!
  AliFemtoDreamv0 *fLambda;               //!
  AliFemtoDreamv0 *fPhiParticle;          //!
  AliFemtoDreamTrack *fTrack;             //!
  AliFemtoDreamEventCuts *fEventCuts;     //
  AliFemtoDreamTrackCuts *fPosKaonCuts;   //
  AliFemtoDreamTrackCuts *fNegKaonCuts;   //
  AliFemtoDreamv0Cuts *fPhiCuts;          //
  AliFemtoDreamv0Cuts *fLambdaCuts;       //
  AliFemtoDreamv0Cuts *fAntiLambdaCuts;   //
  AliFemtoDreamCollConfig *fConfig;       //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamControlSample *fSample;    //!
  AliVTrack **fGTI;                       //!
  int fTrackBufferSize;                   //
  ClassDef(AliAnalysisTaskNanoAODFemtoDreamLambdaPhi, 1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskNanoAODFemtoDreamPhi_H_ */
