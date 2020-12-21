#ifndef AliPHOSEventCounter_h
#define AliPHOSEventCounter_h

#include "AliAnalysisTaskSE.h"

class AliPHOSEventCounter : public AliAnalysisTaskSE {
public:
  AliPHOSEventCounter(const char *name = "AliPHOSEventCounter");
  virtual ~AliPHOSEventCounter() {}
  
  //Standard methods
  void  UserCreateOutputObjects();
  void  UserExec(Option_t *option);
  void  Terminate(Option_t *);

  Int_t GetAnalyzedEvent() {return fEventCounter;}


private:

	Int_t fEventCounter;

  AliPHOSEventCounter(const AliPHOSEventCounter&); // not implemented
  AliPHOSEventCounter& operator=(const AliPHOSEventCounter&); // not implemented

  ClassDef(AliPHOSEventCounter, 1); // PHOS analysis task
};

#endif
