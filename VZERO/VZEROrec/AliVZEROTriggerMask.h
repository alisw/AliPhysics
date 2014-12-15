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
class AliVZERORecoParam;

class AliVZEROTriggerMask : public TObject
{
 public:
                   AliVZEROTriggerMask();   // constructor
   virtual        ~AliVZEROTriggerMask(){}  // destructor

   void FillMasks(AliESDVZERO *esdV0,
		  AliVZEROCalibData *cal,
		  TF1 *slewing);
   Double_t GetZPosition(const char* symname);

   void SetRecoParam(const AliVZERORecoParam *param) { fRecoParam = param; }
   const AliVZERORecoParam* GetRecoParam() const
   {
     if (!fRecoParam) {
       AliFatal("Reco-param object is not set!");
       return NULL;
     }
     return fRecoParam;
   }

private:
   AliVZEROTriggerMask(const AliVZEROTriggerMask& mask);
   AliVZEROTriggerMask& operator = (const AliVZEROTriggerMask& mask);

   Float_t fV0ADist;     // Z position of V0A
   Float_t fV0CDist;     // Z position of V0C
   const AliVZERORecoParam* fRecoParam; //! Pointer to VZERO reco-param object

   ClassDef( AliVZEROTriggerMask, 3 )  // VZERO Trigger Detector class
};

#endif // ALIVZEROTRIGGERMASK_H
