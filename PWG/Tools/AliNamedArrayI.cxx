// $Id: AliNamedArrayI.cxx  $
//
// Named integer array.
//
// Author: S.Aiola

#include "AliNamedArrayI.h"

ClassImp(AliNamedArrayI)

//________________________________________________________________________
AliNamedArrayI::AliNamedArrayI() : 
  TNamed("AliNamedArrayI","AliNamedArrayI"),
  TArrayI()
{
  // Dummy constructor.

}

//________________________________________________________________________
AliNamedArrayI::AliNamedArrayI(const char *name, Int_t n) :
  TNamed(name,name),
  TArrayI(n)
{
  // Standard constructor.
  Clear();
}

//________________________________________________________________________
AliNamedArrayI::AliNamedArrayI(const char *name, Int_t n, const Int_t* array) :
  TNamed(name,name),
  TArrayI(n, array)
{
  // TArrayI copy c-style array constructor.

}

//________________________________________________________________________
AliNamedArrayI::AliNamedArrayI(const char *name, const TArrayI& array) :
  TNamed(name,name),
  TArrayI(array)
{
  // TArrayI copy constructor.
  
}

//________________________________________________________________________
void AliNamedArrayI::Clear(Option_t * /*option*/) 
{ 
  // Clear.

  Reset(-1);
}
