// $Id$
// Category: basic
//
// Class AliMpTransformPadIterator
// -------------------------------
// Composite of iterator and transformer.
// Transforms returned pad. 
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpTransformPadIterator.h"
#include "AliMpTransformer.h"

ClassImp(AliMpTransformPadIterator)

//______________________________________________________________________________
AliMpTransformPadIterator::AliMpTransformPadIterator(AliMpVPadIterator* it, 
                                         const AliMpTransformer* transformer)
  : AliMpVPadIterator(),
    fIterator(it),
    fkTransformer(transformer)
{
// Standard constructor
}

//______________________________________________________________________________
AliMpTransformPadIterator::AliMpTransformPadIterator(
                                         const AliMpTransformPadIterator& right)
  : AliMpVPadIterator(right)
{
// Copy constructor
  
  *this = right; 
}

//______________________________________________________________________________
AliMpTransformPadIterator::AliMpTransformPadIterator()
  : AliMpVPadIterator(),
    fIterator(0),
    fkTransformer(0)
{
// Default constructor
}

//______________________________________________________________________________
AliMpTransformPadIterator::~AliMpTransformPadIterator()
{
// Destructor
// Not owner of its components, does not delete them.
}

//
// operators
//

//______________________________________________________________________________
AliMpTransformPadIterator& 
AliMpTransformPadIterator::operator = (const AliMpTransformPadIterator& right)
{
// assignement operator

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  TObject::operator=(right);
  
  fIterator = right.fIterator;
  fkTransformer = right.fkTransformer;

  return *this;
} 

// 
// public methods 
//

//______________________________________________________________________________
void AliMpTransformPadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the area

  fIterator->First();
}  

//______________________________________________________________________________
void AliMpTransformPadIterator::Next()
{
// Move the iterator to the next valid pad.
// ---
 
  fIterator->Next();
}  

//______________________________________________________________________________
Bool_t AliMpTransformPadIterator::IsDone() const
{
// 

  return fIterator->IsDone();
}  


//______________________________________________________________________________
AliMpPad AliMpTransformPadIterator::CurrentItem() const
{
// Returns current pad with applied transformation.
// ---

  return fkTransformer->Transform(fIterator->CurrentItem());
}  

//______________________________________________________________________________
void AliMpTransformPadIterator::Invalidate()
{
// Set iterator to invalid state.
// ---

  fIterator->Invalidate();
}  

