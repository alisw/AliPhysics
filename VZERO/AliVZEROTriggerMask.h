#ifndef ALIVZEROTriggerMask_H
#define ALIVZEROTriggerMask_H

///_________________________________________________________________________
///
///  Auxiliary classs to compute the  VZERO Trigger
///_________________________________________________________________________   

#include <TObject.h>
#include <TTree.h>
#include <TClonesArray.h>

class AliVZEROTriggerMask : public TObject
{
 public:
                   AliVZEROTriggerMask();   // constructor
   virtual        ~AliVZEROTriggerMask(){}  // destructor

   void FillMasks(TTree* vzeroDigitsTree,
		  TClonesArray* vzeroDigits);
   Double_t GetZPosition(const char* symname);

   void            SetAdcThreshold(Float_t t=55.0) 
     {fAdcThresHold=t;}
   void            SetTimeWindowWidthBBA(Float_t w=50.0)
     {fTimeWindowWidthBBA=w;}
   void            SetTimeWindowWidthBBC(Float_t w=50.0)
     {fTimeWindowWidthBBC=w;}
   void            SetTimeWindowWidthBGA(Float_t w=20.0) 
     {fTimeWindowWidthBGA=w;}
   void            SetTimeWindowWidthBGC(Float_t w=20.0) 
     {fTimeWindowWidthBGC=w;}

   UInt_t GetBBtriggerV0A() { return fBBtriggerV0A;}
   UInt_t GetBGtriggerV0A() { return fBGtriggerV0A;}
   UInt_t GetBBtriggerV0C() { return fBBtriggerV0C;}
   UInt_t GetBGtriggerV0C() { return fBGtriggerV0C;}

private:

   Float_t fAdcThresHold;
   Float_t fTimeWindowWidthBBA; // 
   Float_t fTimeWindowWidthBGA; // 
   Float_t fTimeWindowWidthBBC; // 
   Float_t fTimeWindowWidthBGC; // 
   UInt_t fBBtriggerV0A; // bit mask for Beam-Beam trigger in V0A
   UInt_t fBGtriggerV0A; // bit mask for Beam-Gas trigger in V0A
   UInt_t fBBtriggerV0C; // bit mask for Beam-Beam trigger in V0C
   UInt_t fBGtriggerV0C; // bit mask for Beam-Gas trigger in V0C


   ClassDef( AliVZEROTriggerMask, 1 )  // VZERO Trigger Detector class
};

#endif // AliVZEROTriggerMask_H
