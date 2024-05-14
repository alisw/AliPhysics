#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLAMBDAPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLAMBDAPION_H_
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
class AliAODTrack;

class AliAnalysisTaskLambdaPion : public AliAnalysisTaskSE
{
public:
  enum PCSettings {NoPC, OldPC, NewPC};
  AliAnalysisTaskLambdaPion();
  AliAnalysisTaskLambdaPion(const char *name, bool isMC, PCSettings pcsettings, bool usenolambdaevt);
  virtual ~AliAnalysisTaskLambdaPion();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };

  void SetSphericity(double min, double max)
  {
    fEventCuts->SetSphericityCuts(min, max);
  }

  void SetPairCleaner(AliAnalysisTaskLambdaPion::PCSettings pairCleaner)
  {
    fPCSettings = pairCleaner;
  }
  void SetLambdaCuts(AliFemtoDreamv0Cuts *cuts)
  {
    fLambdaCuts = cuts;
  }
  void SetAntiLambdaCuts(AliFemtoDreamv0Cuts *cuts)
  {
    fAntiLambdaCuts = cuts;
  }
  void SetPosPionCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fPosPionCuts = trkCuts;
  }
  void SetNegPionCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fNegPionCuts = trkCuts;
  }
  void SetMixingEvt(bool usenolambdaevt)
  {
    fUseEvtNoLambda = usenolambdaevt;
  }
  void SetExcludeDausOf(std::vector<UInt_t> motherList) {
    fExcludedMothers = motherList;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config)
  {
    fConfig = config;
  }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }

private:
  AliAnalysisTaskLambdaPion(const AliAnalysisTaskLambdaPion &);
  AliAnalysisTaskLambdaPion &operator=(const AliAnalysisTaskLambdaPion &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                             //
  bool fUseOMixing;                       //
  bool fUseEvtNoLambda;                   //
  std::vector<UInt_t> fExcludedMothers;   // Only valid if run on MC
  PCSettings fPCSettings;                 //
  UInt_t fTrigger;                        //
  TList *fResults;                        //!
  TList *fResultsQA;                      //!
  TList *fQA;                             //!
  TList *fEvtList;                        //!
  TList *fLambdaList;                     //!
  TList *fLambdaMCList;                   //!
  TList *fAntiLambdaList;                 //!
  TList *fAntiLambdaMCList;               //!
  TList *fPionPlusList;                   //!
  TList *fPionPlusMCList;                 //!
  TList *fPionMinusList;                  //!
  TList *fPionMinusMCList;                //!
  AliVEvent *fInputEvent;                 //! current event
  AliFemtoDreamEvent *fEvent;             //!
  AliFemtoDreamv0 *fLambda;               //!
  AliFemtoDreamTrack *fTrack;             //!
  AliFemtoDreamEventCuts *fEventCuts;     //
  AliFemtoDreamTrackCuts *fPosPionCuts;   //
  AliFemtoDreamTrackCuts *fNegPionCuts;   //
  AliFemtoDreamv0Cuts *fLambdaCuts;       //
  AliFemtoDreamv0Cuts *fAntiLambdaCuts;   //
  AliFemtoDreamCollConfig *fConfig;       //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamControlSample *fSample;    //!
  AliAODTrack **fGTI;                     //!
  int fTrackBufferSize;                   //
  ClassDef(AliAnalysisTaskLambdaPion, 2)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKLAMBDAPION_H_ */
