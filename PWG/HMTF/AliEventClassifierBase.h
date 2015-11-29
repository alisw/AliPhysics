#ifndef AliEventClassifierBase_cxx
#define AliEventClassifierBase_cxx

#include "TList.h"
#include "TNamed.h"

#include "AliMCEvent.h"
#include "AliStack.h"

class AliEventClassifierBase : public TNamed {
 public:
  AliEventClassifierBase();
  AliEventClassifierBase(const char* name, const char* title, TList *fTaskOutputList);
  virtual ~AliEventClassifierBase() {}

  Float_t GetClassifierValue(AliMCEvent *event, AliStack *stack);
  void ResetClassifier() {fClassifierValueIsCached = false;}
  TList* GetClassifierOutputList() {return fClassifierOutputList;}
  Int_t GetExpectedMinValue() {return fExpectedMinValue;}
  Int_t GetExpectedMaxValue() {return fExpectedMaxValue;}

 protected:
  virtual void CalculateClassifierValue(AliMCEvent *event, AliStack *stack) = 0;
  Bool_t fClassifierValueIsCached;    // Is the classifier value already computed?
  Float_t fClassifierValue;           // The value for this classifier for the current event
  Int_t fExpectedMinValue;            // The expected min value produced by this estimator, used for hists
  Int_t fExpectedMaxValue;            // The expected max value produced by this estimator, used for hists
  
  TList *fClassifierOutputList;  // The "folder" in which the hists binned in this classifier a saved
  TList *fTaskOutputList;        // The list for the entire task

  // Some generator have funky particles in the stack which break the execution down the line.
  // Check if the particles are of a safe type (should cover like 99%) of the particles present
  std::vector< Int_t > fSafePdgCodes;

  ClassDef(AliEventClassifierBase, 1);
};

#endif
  
