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

   void            SetAdcThreshold(Float_t t=55.0) 
     {fAdcThresHold=t; 
     AliDebug(1,Form("ADC threshold set to %0.2f", fAdcThresHold));}
   
   void            SetTimeWindowWidth(Float_t w=50.0) 
     {fTimeWindowWidthBBA=fTimeWindowWidthBGA
	=fTimeWindowWidthBBC=fTimeWindowWidthBGC=w;}
   void            SetTimeWindowWidthBBA(Float_t w=50.0)
     {fTimeWindowWidthBBA=w;}
   void            SetTimeWindowWidthBBC(Float_t w=50.0)
     {fTimeWindowWidthBBC=w;}
   void            SetTimeWindowWidthBGA(Float_t w=20.0) 
     {fTimeWindowWidthBGA=w;}
   void            SetTimeWindowWidthBGC(Float_t w=20.0) 
     {fTimeWindowWidthBGC=w;}

private:

   Float_t fAdcThresHold;
   Float_t fTimeWindowWidthBBA; // 
   Float_t fTimeWindowWidthBGA; // 
   Float_t fTimeWindowWidthBBC; // 
   Float_t fTimeWindowWidthBGC; // 

   ClassDef( AliVZEROTrigger, 1 )  // VZERO Trigger Detector class
};

#endif // AliVZEROTrigger_H
