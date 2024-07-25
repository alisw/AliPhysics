#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFROMAODLAMBDAPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFROMAODLAMBDAPION_H_
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

class AliAnalysisTaskNanoFromAODLambdaPion : public AliAnalysisTaskSE
{
public:
  enum PCSettings {NoPC, OldPC, NewPC};
  AliAnalysisTaskNanoFromAODLambdaPion();
  AliAnalysisTaskNanoFromAODLambdaPion(const char *name, bool isMC, PCSettings pcsettings, bool usenolambdaevt);
  virtual ~AliAnalysisTaskNanoFromAODLambdaPion();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };

  void SetSphericity(double min, double max)
  {
    fEventCuts->SetSphericityCuts(min, max);
  }
  void SetPairCleaner(AliAnalysisTaskNanoFromAODLambdaPion::PCSettings pairCleaner)
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
  void SetRequiredPDG(int pdg, int identifiedAsPDG)
  {
    fRequiredPDG[pdg] = identifiedAsPDG;
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
  AliAnalysisTaskNanoFromAODLambdaPion(const AliAnalysisTaskNanoFromAODLambdaPion &);
  AliAnalysisTaskNanoFromAODLambdaPion &operator=(const AliAnalysisTaskNanoFromAODLambdaPion &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;                             //
  bool fUseOMixing;                       //
  bool fUseEvtNoLambda;                   //
  std::map<int, int> fRequiredPDG;        // first pdg is the expected one, second is the target. Target=0 (default) accepts any pdg 
  std::vector<UInt_t> fExcludedMothers;   // Only valid if run on MC
  PCSettings fPCSettings;                 //
  UInt_t fTrigger;                        //
  TList *fResults;                        //!
  TList *fResultsQA;                      //!
  TList *fResultsMCGen;                   //!
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
  AliFemtoDreamPairCleaner *fPairCleanerMCGen; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamPartCollection *fPartCollMCGen; //!
  AliFemtoDreamControlSample *fSample;    //!
  AliVTrack **fGTI;                     //!
  int fTrackBufferSize;                   //
  ClassDef(AliAnalysisTaskNanoFromAODLambdaPion, 3)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFROMAODLAMBDAPION_H_ */
