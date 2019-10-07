#ifndef AliAnalysisTaskNanoPt_H
#define AliAnalysisTaskNanoPt_H

#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamControlSample.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TChain.h"
class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoPt : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoPt();
  AliAnalysisTaskNanoPt(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskNanoPt();
  Float_t GetMass2sq(AliFemtoDreamTrack *track);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) { fEvtCuts = cuts; }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fProtonTrack = cuts;
  }
 
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiProtonTrack = cuts;
  }
  void SetDeuteronCuts(AliFemtoDreamTrackCuts *cuts) {
    fDeuteronTrack = cuts;
  }
 
  void SetAntiDeuteronCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiDeuteronTrack = cuts;
  }
 
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }

 private:
  AliAnalysisTaskNanoPt(
      const AliAnalysisTaskNanoPt &task);
  AliAnalysisTaskNanoPt &operator=(
      const AliAnalysisTaskNanoPt &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);

  AliVEvent *fInputEvent;                        //! current event
  AliMCEvent *fMCEvent;                          //! corresponding MC event
  

  TRandom3 *fRandom;  //!

  AliFemtoDreamEvent *fEvent;                        //!
  AliFemtoDreamEventCuts *fEvtCuts;                  //
  AliFemtoDreamTrack *fTrack;                  //
  AliFemtoDreamTrackCuts *fProtonTrack;              //
  AliFemtoDreamTrackCuts* fAntiProtonTrack;          //
  AliFemtoDreamTrackCuts *fDeuteronTrack;              //
  AliFemtoDreamTrackCuts* fAntiDeuteronTrack;          //
  AliFemtoDreamCollConfig *fConfig;                  //
  AliFemtoDreamControlSample *fSample;               //!

  bool fIsMC;              //
  UInt_t fTrigger;         //
//  bool fIsLightweight;     //
 
  int fTrackBufferSize; //
  AliVTrack **fGTI;  //!
  TList* fEvtList;//!
  TList* fProtonList;//! 
  TList* fAntiProtonList;//!   
  TList* fDeuteronList;//! 
  TList* fAntiDeuteronList;//! 

  TH2F* fProtonRestMass; //!
  TH2F* fAntiProtonRestMass; //! 
  TH2F* fDeuteronRestMass; //!
  TH2F* fAntiDeuteronRestMass; //!
  ClassDef(AliAnalysisTaskNanoPt,5)
};
#endif
