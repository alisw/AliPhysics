//
// *** Class AliRsnCutValue ***
//
// This cut implementation can be used to cut generically on
// any value which can be computed from AliRsnValue class.
// Since that value is implemented always as a Double_t one,
// then this cut operates only with the Double_t data members
// of the AliRsnCut base class.
// It allows to customize the reference AliRsnValue object by
// means of a getter that returns a pointer to it.
// This cut can apply to any kind of object, but the type of
// target must be one of those for which the chosen value type
// makes sense to be computed
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTVALUE_H
#define ALIRSNCUTVALUE_H

#include "AliRsnCut.h"
#include "AliRsnValue.h"

class AliRsnCutValue : public AliRsnCut {
public:

   AliRsnCutValue();
   AliRsnCutValue(const char *name, Double_t min, Double_t max, Bool_t isMC);
   AliRsnCutValue(const AliRsnCutValue& copy);
   AliRsnCutValue& operator=(const AliRsnCutValue& copy);
   virtual ~AliRsnCutValue() { }

   Double_t       GetComputedValue()              {if (fValue) return fValue->GetComputedValue(); return -1E20;}
   AliRsnValue*   GetValueObj()                   {return fValue;}
   void           SetValueObj(AliRsnValue *value) {fValue = value; SetTargetType(value->GetTargetType());}
   Bool_t         IsUsingMC()                     {return fUseMC;}
   void           UseMC(Bool_t yn = kTRUE)        {fUseMC = yn;}

   virtual Bool_t IsSelected(TObject *object);
   virtual void   Print(const Option_t *option = "") const;

protected:

   Bool_t       fUseMC;
   AliRsnValue *fValue;

   ClassDef(AliRsnCutValue, 1)
};

#endif
