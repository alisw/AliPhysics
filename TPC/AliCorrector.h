#ifndef ALI_CORRECTOR_H
#define ALI_CORRECTOR_H

#include "TObject.h"

class AliCorrector:public TObject {
public:
  virtual ~AliCorrector() {};
  virtual void Correct(const Double_t *position,Double_t *corrected)=0;
  virtual void CorrectInverse(const Double_t *position,Double_t *corrected) {
    Correct(position,corrected);
    for (Int_t i=0;i<3;++i)
      corrected[i]=position[i]-(corrected[i]-position[i]);
  }
  ClassDef(AliCorrector,1)
};

#endif
