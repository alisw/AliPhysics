#ifndef AliObservableCorrelationsOfClassifiers_cxx
#define AliObservableCorrelationsOfClassifiers_cxx

#include "AliObservableBase.h"

class AliObservableCorrelationsOfClassifiers : public AliObservableBase {
 public:
  AliObservableCorrelationsOfClassifiers();
  AliObservableCorrelationsOfClassifiers(AliEventClassifierBase *classifier0, AliEventClassifierBase *classifier1);
  ~AliObservableCorrelationsOfClassifiers() {};

  void Fill(AliMCEvent *event, AliStack *stack);

 private:
  TH2F *fhistogram;
  AliEventClassifierBase *fclassifier0;
  AliEventClassifierBase *fclassifier1;

  ClassDef(AliObservableCorrelationsOfClassifiers, 1);
};

#endif
