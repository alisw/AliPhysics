#ifndef ALIADDECISION_H
#define ALIADDECISION_H

///_________________________________________________________________________
///
///  Auxiliary classs to compute the AD Trigger
///_________________________________________________________________________   

#include <TObject.h>
#include "AliLog.h"

class TTree;
class TClonesArray;
class TF1;
class AliESDAD;
class AliADCalibData;
class AliADRecoParam;

class AliADDecision : public TObject
{
 public:
  AliADDecision();   // constructor
  virtual        ~AliADDecision();

   void FillDecisions(AliESDAD *esdAD);
   Double_t GetZPosition(const char* symname);

   void SetRecoParam(const AliADRecoParam *param) { fRecoParam = param; }
   const AliADRecoParam* GetRecoParam() const
   {
     if (!fRecoParam) {
       AliError("Reco-param object is not set!");
       return NULL;
     }
     return fRecoParam;
   }

private:
   AliADDecision(const AliADDecision& mask);
   AliADDecision& operator = (const AliADDecision& mask);

   Float_t fADADist;     // Z position of ADA
   Float_t fADCDist;     // Z position of ADC
   const AliADRecoParam* fRecoParam; //! Pointer to AD reco-param object
   TF1 *fEarlyHitCutShape; //! Shape of cut on early hits

   ClassDef( AliADDecision, 2)  // AD Offline trigger class
};

#endif // ALIADDECISION_H
