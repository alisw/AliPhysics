#ifndef ALITPC_EXB
#define ALITPC_EXB

#include "AliCorrector.h"

class AliTPCExB:public AliCorrector {
public:
  virtual ~AliTPCExB() {};
  void SetDriftVelocity(Double_t driftVelocity) {
    fDriftVelocity=driftVelocity;
  };
protected:
  Double_t fDriftVelocity; // The electron drift velocity.
  
  ClassDef(AliTPCExB,1)
};

#endif
