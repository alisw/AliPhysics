#ifndef ALIRHOPARAMETER_H
#define ALIRHOPARAMETER_H

// $Id: $

class TString;
class TF1;

#include <TParameter.h>

class AliRhoParameter : public TParameter<Double_t> {
 public: 
  AliRhoParameter();
  AliRhoParameter(const char *name, Double_t val=0);
  void        Clear(Option_t *option="");

private:
  AliRhoParameter(const AliRhoParameter&);             // not implemented
  AliRhoParameter& operator=(const AliRhoParameter&);  // not implemented
  
  ClassDef(AliRhoParameter, 1); // Rho parameter
};
#endif
