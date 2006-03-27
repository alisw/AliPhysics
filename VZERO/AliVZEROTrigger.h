#ifndef ALIVZEROTrigger_H
#define ALIVZEROTrigger_H

///_________________________________________________________________________
///
///  Class for making  VZERO Trigger
///_________________________________________________________________________   


#include "AliTriggerDetector.h"
#include "AliTriggerInput.h"

#include "AliVZEROLoader.h"
#include "AliVZEROdigit.h"

#include "AliLog.h"


class AliVZEROTrigger : public AliTriggerDetector
{
 public:
                   AliVZEROTrigger();   // constructor
   virtual        ~AliVZEROTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

   void            SetAdcThreshold(Float_t t=62.5) 
     {fAdcThresHold=t; AliDebug(1,Form("ADC threshold set to %0.2f", fAdcThresHold));}
   
   void            SetTimeWindowWidth(Float_t w=2) {fTimeWindowWidth=w;}


private:

   Float_t fAdcThresHold;
   Float_t fTimeWindowWidth; // 

   ClassDef( AliVZEROTrigger, 1 )  // VZERO Trigger Detector class
};

#endif // AliVZEROTrigger_H
