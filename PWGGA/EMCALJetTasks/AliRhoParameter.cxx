// $Id: $
//
// Parameter class with a Clear.
//
// Author: C.Loizides

#include "AliRhoParameter.h"

ClassImp(AliRhoParameter)

//________________________________________________________________________
AliRhoParameter::AliRhoParameter() : 
  TParameter<Double_t>()
{
  // Dummy constructor.
}

//________________________________________________________________________
AliRhoParameter::AliRhoParameter(const char *name, Double_t val) :
  TParameter<Double_t>(name,val)
{
  // Constructor.
}

//________________________________________________________________________
void AliRhoParameter::Clear(Option_t * /*option*/) 
{ 
  // Clear.

  SetVal(0);
}
