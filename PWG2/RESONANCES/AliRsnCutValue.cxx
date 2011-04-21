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
   fValue(0x0)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue
(const char *name, Double_t min, Double_t max) :
   AliRsnCut(name, AliRsnTarget::kTargetTypes, min, max),
   fValue(0x0)
{
//
// Main constructor.
// Sets the AliRsnValue data member accordingly to arguments passed here.
// NOTE: if the value needs a support object, it must be passed separately
//       using the GetValueObje() of this class
//
}

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue(const AliRsnCutValue& copy) :
   AliRsnCut(copy),
   fValue(copy.fValue)
{
//
// Copy constructor.
// Does not duplicate memory allocation.
//
}

//_________________________________________________________________________________________________
AliRsnCutValue& AliRsnCutValue::operator=(const AliRsnCutValue& copy)
{
//
// Assignment operator.
// Does not duplicate memory allocation.
//

   AliRsnCut::operator=(copy);
   fValue = copy.fValue;

   return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutValue::IsSelected(TObject *object)
{
//
// Checks the cut.
// Calls the AliRsnValue::Eval() method and then checks its output.
//

   // skip cut if value is not initialized
   if (!fValue) return kTRUE;
   
   // match target types
   SetTargetType(fValue->GetTargetType());

   // try to compute values
   Bool_t success = fValue->Eval(object);

   // check success
   if (!success) {
      AliWarning(Form("[%s] Failed to compute value", GetName()));
      return kFALSE;
   }

   // check in range
   fCutValueD = fValue->GetComputedValue();
   return OkRangeD();
}

//_________________________________________________________________________________________________
void AliRsnCutValue::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name   : %s", GetName()));
   AliInfo(Form("Cut value  : %s", fValue->GetName()));
   AliInfo(Form("Cut range  : %f - %f", fMinD, fMaxD));
}
