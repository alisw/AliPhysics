#ifndef ALITPC_EXB
#define ALITPC_EXB

#include "TObject.h"

class AliTPCExB:public TObject {
public:
  virtual ~AliTPCExB() {};
  virtual void Correct(const Double_t *position,Double_t *corrected)=0;
  virtual void CorrectInverse(const Double_t *position,Double_t *corrected) {
    Correct(position,corrected);
    for (Int_t i=0;i<3;++i)
      corrected[i]=position[i]-(corrected[i]-position[i]);
  }

  ClassDef(AliTPCExB,0)
};

#endif
