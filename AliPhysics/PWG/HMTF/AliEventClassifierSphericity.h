#ifndef AliEventClassifierSphericity_cxx
#define AliEventClassifierSphericity_cxx

#include "AliEventClassifierBase.h"

class AliEventClassifierSphericity : public AliEventClassifierBase {
 public:
  AliEventClassifierSphericity()
    : AliEventClassifierBase() {}
  AliEventClassifierSphericity(const char* name, const char* title,
			TList *taskOutputList);
  virtual ~AliEventClassifierSphericity() {}

 private:
  void CalculateClassifierValue(AliMCEvent *event, AliStack *stack);
  
  ClassDef(AliEventClassifierSphericity, 1);
};

#endif
  
