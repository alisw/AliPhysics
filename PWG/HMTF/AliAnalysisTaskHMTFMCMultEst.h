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
  void SetGlobalTrigger(const char* globalTrigger) {fGlobalTriggerName = globalTrigger;};

  std::vector<AliEventClassifierBase*> fClassifiers;
  std::vector<AliObservableBase*> fObservables;

 private:
  TList *fMyOut;                          // Output list
  TString fGlobalTriggerName;             // eg. INEL>0
  Float_t fGlobalTriggerMinValue;         // Smallest value of the trigger classifier to accept an event
  AliEventClassifierBase *fGlobalTrigger; // A classifier used to extract the value triggered on

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskHMTFMCMultEst(const AliAnalysisTaskHMTFMCMultEst&); // not implemented
  AliAnalysisTaskHMTFMCMultEst& operator=(const AliAnalysisTaskHMTFMCMultEst&); // not implemented

  ClassDef(AliAnalysisTaskHMTFMCMultEst, 1); // example of analysis
};

#endif
