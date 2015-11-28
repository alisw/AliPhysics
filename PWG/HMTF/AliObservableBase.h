#ifndef AliObservableBase_cxx
#define AliObservableBase_cxx

#include "AliEventClassifierBase.h"

class AliObservableBase : public TNamed {
 public:
  AliObservableBase();
  AliObservableBase(const char* name, const char* title);
  ~AliObservableBase() {};
  virtual void Fill(AliMCEvent *event, AliStack *stack) = 0;

 protected:
  // Some generator have funky particles in the stack which break the execution down the line.
  // Check if the particles are of a safe type (should cover like 99%) of the particles present
  std::vector< Int_t > fSafePdgCodes;

  ClassDef(AliObservableBase, 1);
};

#endif

