#ifndef AliAnalysisTaskHMTFMCMultEst_cxx
#define AliAnalysisTaskHMTFMCMultEst_cxx

#include "AliAnalysisTaskSE.h"

#include "AliEventClassifierBase.h"
#include "AliObservableBase.h"

class AliAnalysisTaskHMTFMCMultEst : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHMTFMCMultEst();
  AliAnalysisTaskHMTFMCMultEst(const char *name );
  virtual ~AliAnalysisTaskHMTFMCMultEst() {};

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // set the name of the global trigger in the AddTask marcro through this function
  void SetGlobalTrigger(Int_t triggerEnum) {fGlobalTrigger = triggerEnum;};
  void SetGlobalSystem(Int_t systemEnum) {fGlobalSystem = systemEnum;}
 private:
  TList *fMyOut;                          // Output list
  std::vector<AliEventClassifierBase*> fClassifiers;
  std::vector<AliObservableBase*> fObservables;

  Int_t fGlobalTrigger;
  enum {kINEL, kINELGT0, kV0AND};
  Int_t fGlobalSystem;
  enum {kPP, kPPB, kPBPB};

  void SetupInelAsGlobalTrigger();
  void SetupInelGt0AsGlobalTrigger(AliEventClassifierBase* etaLt1);
  void SetupV0ANDAsGlobalTrigger(AliEventClassifierBase* V0A, AliEventClassifierBase* V0C);

  Bool_t IsInel(AliMCEvent *event, AliStack *stack);
  Bool_t IsInelGt0(AliMCEvent *event, AliStack *stack);
  Bool_t IsV0AND(AliMCEvent *event, AliStack *stack);
  // vector to save the classifiers used in the global trigger
  std::vector<AliEventClassifierBase*> fGlobalTriggerClassifiers;

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskHMTFMCMultEst(const AliAnalysisTaskHMTFMCMultEst&); // not implemented
  AliAnalysisTaskHMTFMCMultEst& operator=(const AliAnalysisTaskHMTFMCMultEst&); // not implemented

  ClassDef(AliAnalysisTaskHMTFMCMultEst, 2); // example of analysis
};

#endif
