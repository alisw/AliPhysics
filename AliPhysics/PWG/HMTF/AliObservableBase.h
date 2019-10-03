#ifndef AliObservableBase_cxx
#define AliObservableBase_cxx

#include "AliEventClassifierBase.h"

class AliObservableBase : public TNamed {
 public:
  AliObservableBase()
    : TNamed() {}
  AliObservableBase(const char* name, const char* title)
    : TNamed(name, title) {}
  ~AliObservableBase() {}

  virtual void Fill(AliMCEvent *event, AliStack *stack) = 0;

  ClassDef(AliObservableBase, 1);
};

#endif

