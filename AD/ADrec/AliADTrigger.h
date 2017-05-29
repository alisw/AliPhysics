// -*- C++ -*-
#ifndef ALIADTrigger_H
#define ALIADTrigger_H
// ---------------------
// Class AliADTrigger
// ---------------------
// Top class to simulate the AD trigger response
// This class is only used for interface with AliTriggerDetector
// Its create and Set  Inputs of the CTP
// The Calculation of the trigger response is done into AliADTriggerSimulator
//


#include "AliTriggerDetector.h"

class AliADTrigger : public AliTriggerDetector
{
public:
  AliADTrigger();
  virtual        ~AliADTrigger() {}
  virtual void    CreateInputs();
  virtual void    Trigger();

  ClassDef( AliADTrigger, 1);  // AD Trigger Detector class
};

#endif // AliADTrigger_H


