#ifndef ALIVZEROTrigger_H
#define ALIVZEROTrigger_H
// ---------------------
// Class AliVZEROTrigger
// ---------------------
// Top class to simulate the VZERO trigger response
// This class is only used for interface with AliTriggerDetector
// Its create and Set  Inputs of the CTP
// The Calculation of the trigger response is done into AliVZEROTriggerSimulator
//


#include "AliTriggerDetector.h"

class AliVZEROTrigger : public AliTriggerDetector
{
 public:
                   AliVZEROTrigger();   // constructor
   virtual        ~AliVZEROTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

   ClassDef( AliVZEROTrigger, 2 )  // VZERO Trigger Detector class
};

#endif // AliVZEROTrigger_H


