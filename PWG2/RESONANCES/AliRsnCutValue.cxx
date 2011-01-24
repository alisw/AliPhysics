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

#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnEvent.h"
#include "AliRsnPairDef.h"

#include "AliRsnCutValue.h"

ClassImp(AliRsnCutValue)

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue() :
  AliRsnCut(),
  fValue(),
  fPairDef(0x0)
{
//
// Default constructor.
//

  SetTargetType(fValue.GetTargetType());
}

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue
(const char *name, AliRsnValue::EValueType type, Double_t min, Double_t max, AliRsnPairDef *pd) :
  AliRsnCut(name, AliRsnTarget::kTargetTypes, min, max),
  fValue(Form("val_%s", name), type),
  fPairDef(pd)
{
//
// Main constructor.
// Recalls the setter for the value type of the AliRsnValue data member,
// which determines also the type of target to be expected
//

  SetTargetType(fValue.GetTargetType());
}

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue(const AliRsnCutValue& copy) :
  AliRsnCut(copy),
  fValue(copy.fValue),
  fPairDef(copy.fPairDef)
{
//
// Copy constructor.
// Does not duplicate memory allocation.
//

  SetTargetType(fValue.GetTargetType());
}

//_________________________________________________________________________________________________
AliRsnCutValue& AliRsnCutValue::operator=(const AliRsnCutValue& copy)
{
//
// Assignment operator.
// Does not duplicate memory allocation.
//

  AliRsnCut::operator=(copy);
  
  fValue   = copy.fValue;
  fPairDef = copy.fPairDef;
  SetTargetType(fValue.GetTargetType());
  
  return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutValue::IsSelected(TObject *object)
{
//
// Checks the cut.
// Calls the AliRsnValue::Eval() method and then checks its output.
//

  // make sure that target of this object matches that
  // of the inserted value object
  SetTargetType(fValue.GetTargetType());
  
  // try to compute values
  Bool_t success = fValue.Eval(object);
  
  // check success
  if (!success)
  {
    AliWarning(Form("[%s] Failed to compute value", GetName()));
    return kFALSE;
  }
  
  // check in range
  fCutValueD = fValue.GetComputedValue();
  return OkRangeD();
}

//_________________________________________________________________________________________________
void AliRsnCutValue::Print(const Option_t *) const
{
//
// Print information on this cut
//

  AliInfo(Form("Cut name   : %s", GetName()));
  AliInfo(Form("Cut value  : %s", fValue.GetValueTypeName()));
  AliInfo(Form("Cut range  : %f - %f", fMinD, fMaxD));
}
