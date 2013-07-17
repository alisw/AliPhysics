// $Id$
//
// Named string object.
//
// Author: S.Aiola

#include "AliNamedString.h"

ClassImp(AliNamedString)

//________________________________________________________________________
AliNamedString::AliNamedString() : 
  TObjString(),
  fName()
{
  // Dummy constructor.

}

//________________________________________________________________________
AliNamedString::AliNamedString(const char *name, const char *string) :
  TObjString(string),
  fName(name)
{
  // Standard constructor.

}
