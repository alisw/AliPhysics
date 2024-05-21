#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLAMBDAPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLAMBDAPION_H_
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

class AliAnalysisTaskNanoLambdaPion : public AliAnalysisTaskSE
{
public:
  enum PCSettings {NoPC, OldPC, NewPC};
  AliAnalysisTaskNanoLambdaPion();
  AliAnalysisTaskNanoLambdaPion(const char *name, bool isMC, PCSettings pcsettings, bool usenolambdaevt);
  virtual ~AliAnalysisTaskNanoLambdaPion();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };

  void SetSphericity(double min, double max)
  {
    fEventCuts->SetSphericityCuts(min, max);
  }
  void SetPairCleaner(AliAnalysisTaskNanoLambdaPion::PCSettings pairCleaner)
  {
    fPCSettings = pairCleaner;
  }
  void SetMinvLambdaKStar(bool doIt)
  {
    cout << "Invariant mass as a function of kstar will be computed" << endl;
    fMinvLambdaKStar = doIt;
  }
  void SetPionsForMinvLambdaKStar(std::vector<AliFemtoDreamBasePart> pionsPlus,
                                  std::vector<AliFemtoDreamBasePart> pionsMinus)
  {
    cout << "Setting pions for kstar of Lambda candidates" << endl; 
    fPionsPlusForMinvV0KStar = pionsPlus;
    fPionsMinusForMinvV0KStar = pionsMinus;
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
  AliAnalysisTaskNanoLambdaPion(const AliAnalysisTaskNanoLambdaPion &);
  AliAnalysisTaskNanoLambdaPion &operator=(const AliAnalysisTaskNanoLambdaPion &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  std::vector<AliFemtoDreamBasePart> fPionsPlusForMinvV0KStar;
  std::vector<AliFemtoDreamBasePart> fPionsMinusForMinvV0KStar;
  bool fIsMC;                             //
  bool fUseOMixing;                       //
  bool fUseEvtNoLambda;                   //
  bool fMinvLambdaKStar;                  //
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
  AliVTrack **fGTI;                       //!
  int fTrackBufferSize;                   //
  ClassDef(AliAnalysisTaskNanoLambdaPion, 2)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOLAMBDAPION_H_ */
