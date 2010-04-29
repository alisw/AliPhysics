#ifndef ALIVZEROTRIGGERMASK_H
#define ALIVZEROTRIGGERMASK_H

///_________________________________________________________________________
///
///  Auxiliary classs to compute the  VZERO Trigger
///_________________________________________________________________________   

#include <TObject.h>
class TTree;
class TClonesArray;
class AliESDVZERO;
class AliVZEROCalibData;

class AliVZEROTriggerMask : public TObject
{
 public:
                   AliVZEROTriggerMask();   // constructor
   virtual        ~AliVZEROTriggerMask(){}  // destructor

   void FillMasks(AliESDVZERO *esdV0,
		  AliVZEROCalibData *cal,
		  TF1 *slewing);
   Double_t GetZPosition(const char* symname);

   void            SetAdcThreshold(Float_t thr) 
     {fAdcThresHold=thr;}
   void            SetTimeWindowBBA(Float_t wlow,Float_t wup)
     {fTimeWindowBBALow=wlow; fTimeWindowBBAUp=wup;}
   void            SetTimeWindowBBC(Float_t wlow,Float_t wup)
     {fTimeWindowBBCLow=wlow; fTimeWindowBBCUp=wup;}
   void            SetTimeWindowBGA(Float_t wlow,Float_t wup) 
     {fTimeWindowBGALow=wlow; fTimeWindowBGAUp=wup;}
   void            SetTimeWindowBGC(Float_t wlow,Float_t wup) 
     {fTimeWindowBGCLow=wlow; fTimeWindowBGCUp=wup;}

private:

   Float_t fAdcThresHold; // Threshold on the ADC
   Float_t fTimeWindowBBALow;  // BBA window (lower cut)
   Float_t fTimeWindowBBAUp;   // BBA window (upper cut)
   Float_t fTimeWindowBGALow;  // BGA window (lower cut)
   Float_t fTimeWindowBGAUp;   // BGA window (upper cut)
   Float_t fTimeWindowFakeALow;// Fake V0A window (lower cut)
   Float_t fTimeWindowFakeAUp; // Fake V0A window (upper cut)
   Float_t fTimeWindowBBCLow;  // BBC window (lower cut)
   Float_t fTimeWindowBBCUp;   // BBC window (upper cut)
   Float_t fTimeWindowBGCLow;  // BGC window (lower cut)
   Float_t fTimeWindowBGCUp;   // BGC window (upper cut)
   Float_t fTimeWindowFakeCLow;// Fake V0C window (lower cut)
   Float_t fTimeWindowFakeCUp; // Fake V0C window (upper cut)

   Float_t fV0ADist;     // Z position of V0A
   Float_t fV0CDist;     // Z position of V0C


   ClassDef( AliVZEROTriggerMask, 2 )  // VZERO Trigger Detector class
};

#endif // ALIVZEROTRIGGERMASK_H
