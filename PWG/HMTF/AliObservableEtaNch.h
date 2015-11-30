#ifndef AliObservableEtaNch_cxx
#define AliObservableEtaNch_cxx

#include "AliObservableBase.h"

class AliObservableEtaNch : public AliObservableBase {
 public:
  AliObservableEtaNch();
  AliObservableEtaNch(AliEventClassifierBase *classifier);
  ~AliObservableEtaNch() {};

  void Fill(AliMCEvent *event, AliStack *stack);
 private:
  TH2F *fhistogram;
  AliEventClassifierBase *fclassifier;

  ClassDef(AliObservableEtaNch, 1);
};

#endif
