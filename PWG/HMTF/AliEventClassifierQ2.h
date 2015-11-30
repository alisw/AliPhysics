#ifndef AliEventClassifierQ2_cxx
#define AliEventClassifierQ2_cxx

#include "AliEventClassifierBase.h"

class AliEventClassifierQ2 : public AliEventClassifierBase {
 public:
  AliEventClassifierQ2()
    : AliEventClassifierBase() {}
  AliEventClassifierQ2(const char* name, const char* title,
			TList *taskOutputList);
  virtual ~AliEventClassifierQ2() {}

 private:
  void CalculateClassifierValue(AliMCEvent *event, AliStack *stack);
  
  ClassDef(AliEventClassifierQ2, 1);
};

#endif
  
