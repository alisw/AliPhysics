#ifndef AliEventClassifierBase_cxx
#define AliEventClassifierBase_cxx

#include "TList.h"
#include "TNamed.h"

#include "AliMCEvent.h"
#include "AliStack.h"

class AliEventClassifierBase : public TNamed {
 public:
  AliEventClassifierBase();
  AliEventClassifierBase(const char* name, const char* title, TList *fTaskOutputList, Int_t fCollisionSystem = 0);
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
  Int_t fCollisionSystem;             // 0 = pp, 1 = ppb, 2 = pbpb. Only used for some derived classes (e.g. multiplicity classifier) , some others (e.g. spherocity) are system 'agnostic' by definition
  
  TList *fClassifierOutputList;  // The "folder" in which the hists binned in this classifier a saved
  TList *fTaskOutputList;        // The list for the entire task

  ClassDef(AliEventClassifierBase, 2);
};

#endif
  
