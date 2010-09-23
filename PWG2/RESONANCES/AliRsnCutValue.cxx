//
// Class AliRsnCutValue
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <TMath.h>
#include <TLorentzVector.h>

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnMother.h"
#include "AliRsnEvent.h"
#include "AliRsnValue.h"

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
}

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue
(const char *name, ETarget target, Double_t min, Double_t max, AliRsnPairDef *pd) :
  AliRsnCut(name, target, min, max),
  fValue(Form("val_%s", name), AliRsnValue::kValueTypes),
  fPairDef(pd)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutValue::AliRsnCutValue(const AliRsnCutValue& copy) :
  AliRsnCut(copy),
  fValue(copy.fValue),
  fPairDef(copy.fPairDef)
{
//
// Copy constructor
//
}

//_________________________________________________________________________________________________
AliRsnCutValue& AliRsnCutValue::operator=(const AliRsnCutValue& copy)
{
//
// Assignment operator
//

  (*this)  = copy;
  fValue   = copy.fValue;
  fPairDef = copy.fPairDef;
  
  return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutValue::IsSelected(TObject *obj1, TObject * /*obj2*/)
{
//
// Checks the cut.
// Calls the appropriate AliRsnValue::Eval() method
// depending on the type of passed object.
// It is up to the user to be sure that the association is meaningful
//

  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(obj1);
  AliRsnMother   *mother   = dynamic_cast<AliRsnMother*>(obj1);
  
  if (daughter)
  {
    if (!fValue.Eval(daughter, fEvent)) return kFALSE;
    fCutValueD = fValue.GetValue();
  }
  else if (mother)
  {
    if (!fValue.Eval(mother, fPairDef, fEvent)) return kFALSE;
    fCutValueD = fValue.GetValue();
  }
  
  return OkRangeD();
}
