//
// *** Class AliRsnCutValue ***
//
// This cut implementation can be used to cut generically on 
// any value which can be computed from AliRsnValue class.
// Since that value is implemented always as a Double_t one,
// then this cut operates only with the Double_t data members
// of the AliRsnCut base class.
// It allows to cusomize the reference AliRsnValue object by
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

class AliRsnPairDef;

class AliRsnCutValue : public AliRsnCut
{
  public:

    AliRsnCutValue();
    AliRsnCutValue(const char *name, AliRsnValue::EValueType type, Double_t min, Double_t max, AliRsnPairDef *pd = 0x0);
    AliRsnCutValue(const AliRsnCutValue& copy);
    AliRsnCutValue& operator=(const AliRsnCutValue& copy);
    virtual ~AliRsnCutValue() { }

    void           SetPairDef(AliRsnPairDef *pd) {fPairDef = pd;}
    AliRsnPairDef* GetPairDef()  {return fPairDef;}
    Double_t       GetValue()    {return fValue.GetComputedValue();}
    AliRsnValue*   GetValueObj() {return &fValue;}
    
    virtual Bool_t IsSelected(TObject *object);
    virtual void   Print(const Option_t *option = "") const;

  protected:
  
    AliRsnValue    fValue;
    AliRsnPairDef *fPairDef;

    ClassDef(AliRsnCutValue, 1)
};

#endif
