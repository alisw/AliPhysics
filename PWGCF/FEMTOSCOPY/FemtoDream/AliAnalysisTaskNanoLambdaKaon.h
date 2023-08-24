#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLAMBDAKAON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLAMBDAKAON_H_
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

class AliAnalysisTaskNanoLambdaKaon : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskNanoLambdaKaon();
  AliAnalysisTaskNanoLambdaKaon(const char *name, bool isMC, bool isNewPC);
  virtual ~AliAnalysisTaskNanoLambdaKaon();
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
  void SetCollectionConfig(AliFemtoDreamCollConfig *config)
  {
    fConfig = config;
  }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }

private:
  AliAnalysisTaskNanoLambdaKaon(const AliAnalysisTaskNanoLambdaKaon &);
  AliAnalysisTaskNanoLambdaKaon &operator=(const AliAnalysisTaskNanoLambdaKaon &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;                             //
  bool fUseOMixing;                       //
  bool fIsNewPC;                         //
  UInt_t fTrigger;                        //
  TList *fResults;                        //!
  TList *fResultsQA;                      //!
  TList *fQA;                             //!
  TList *fEvtList;                        //!
  TList *fLambdaList;                     //!
  TList *fLambdaMCList;                   //!
  TList *fAntiLambdaList;                 //!
  TList *fAntiLambdaMCList;               //!
  TList *fKaonPlusList;                   //!
  TList *fKaonPlusMCList;                 //!
  TList *fKaonMinusList;                  //!
  TList *fKaonMinusMCList;                //!
  AliVEvent *fInputEvent;                 //! current event
  AliFemtoDreamEvent *fEvent;             //!
  AliFemtoDreamv0 *fLambda;               //!
  AliFemtoDreamTrack *fTrack;             //!
  AliFemtoDreamEventCuts *fEventCuts;     //
  AliFemtoDreamTrackCuts *fPosKaonCuts;   //
  AliFemtoDreamTrackCuts *fNegKaonCuts;   //
  AliFemtoDreamv0Cuts *fLambdaCuts;       //
  AliFemtoDreamv0Cuts *fAntiLambdaCuts;   //
  AliFemtoDreamCollConfig *fConfig;       //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamControlSample *fSample;    //!
  AliVTrack **fGTI;                       //!
  int fTrackBufferSize;                   //
  ClassDef(AliAnalysisTaskNanoLambdaKaon, 1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLAMBDAKAON_H_ */
