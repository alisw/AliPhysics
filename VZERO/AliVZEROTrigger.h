#ifndef ALIVZEROTrigger_H
#define ALIVZEROTrigger_H

///_________________________________________________________________________
///
///  Class for making  VZERO Trigger
///_________________________________________________________________________   


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

