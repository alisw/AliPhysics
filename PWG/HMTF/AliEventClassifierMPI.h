#ifndef AliEventClassifierMPI_cxx
#define AliEventClassifierMPI_cxx

#include "AliEventClassifierBase.h"

class AliEventClassifierMPI : public AliEventClassifierBase {
 public:
  AliEventClassifierMPI()
    : AliEventClassifierBase() {}
  AliEventClassifierMPI(const char* name, const char* title,
			TList *taskOutputList);
  virtual ~AliEventClassifierMPI() {}

 private:
  void CalculateClassifierValue(AliMCEvent *event, AliStack *stack);
  
  ClassDef(AliEventClassifierMPI, 1);
};

#endif
  
