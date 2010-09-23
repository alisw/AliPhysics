//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTVALUE_H
#define ALIRSNCUTVALUE_H

#include "AliRsnCut.h"

class AliRsnDaughter;
class AliRsnMother;
class AliRsnEvent;

class AliRsnCutValue : public AliRsnCut
{
  public:

    AliRsnCutValue();
    AliRsnCutValue(const char *name, ETarget target, Double_t min, Double_t max, AliRsnPairDef *pd = 0x0);
    AliRsnCutValue(const AliRsnCutValue& copy);
    AliRsnCutValue& operator=(const AliRsnCutValue& copy);
    virtual ~AliRsnCutValue() { }

    void           SetPairDef(AliRsnPairDef *pd) {fPairDef = pd;}
    AliRsnValue*   GetRsnValue() {return &fValue;}
    Double_t       GetCutValue() {return fValue.GetValue();}
    virtual Bool_t IsSelected(TObject *obj1, TObject *obj2 = 0x0);

  protected:
  
    AliRsnValue    fValue;
    AliRsnPairDef *fPairDef;

    ClassDef(AliRsnCutValue, 1)
};

#endif
