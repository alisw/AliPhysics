// $Id$
// Category: plane
//
// Class AliMpPlaneAreaPadIterator
// -------------------------------
// Class, which defines an iterator over the pads 
// inside a given area in a plane in horizontal direction.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TVector2.h>

#include "AliMpPlaneAreaPadIterator.h"
#include "AliMpTransformPadIterator.h"
#include "AliMpPlaneSegmentation.h"
#include "AliMpSectorSegmentation.h"

ClassImp(AliMpPlaneAreaPadIterator)

//______________________________________________________________________________
AliMpPlaneAreaPadIterator::AliMpPlaneAreaPadIterator(
                                const AliMpPlaneSegmentation* segmentation,
                                const AliMpArea& area) 
 : AliMpVPadIterator(),
   fkPlaneSegmentation(segmentation),
   fkArea(area),
   fPadIterators()
{
// Normal constructor, start in invalid position
 
  DecomposeArea();

  fCurrentIterator = fPadIterators.end();
}

//______________________________________________________________________________
AliMpPlaneAreaPadIterator::AliMpPlaneAreaPadIterator(
                                const AliMpPlaneAreaPadIterator& right)
  : AliMpVPadIterator(right)
{
// Copy constructor
 
  Fatal("Copy constructor", "Not implemented");
}

//______________________________________________________________________________
AliMpPlaneAreaPadIterator::AliMpPlaneAreaPadIterator()
 : AliMpVPadIterator(),
   fkPlaneSegmentation(0),
   fkArea(AliMpArea()),
   fPadIterators()
{
// Dummy default constructor.
}

//______________________________________________________________________________
AliMpPlaneAreaPadIterator::~AliMpPlaneAreaPadIterator()
{
// Destructor

  // delete created iterators here
}

//
// operators
//

//______________________________________________________________________________
AliMpPlaneAreaPadIterator& 
AliMpPlaneAreaPadIterator::operator = (const AliMpPlaneAreaPadIterator& right)
{
// Assignement operator

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");

  return *this;
} 

//
// private methods
//

//______________________________________________________________________________
void AliMpPlaneAreaPadIterator::DecomposeArea()
{
// Decompose the area into areas belonging to the quadrants.
// --

  for (Int_t i=0; i<fkPlaneSegmentation->GetNofTransformers(); i++) {
  
    AliMpTransformer* transformer = fkPlaneSegmentation->GetTransformer(i);
    AliMpArea area = transformer->CutArea(fkArea);

    if (area.IsValid()) {
    
      AliMpSectorSegmentation* segmentation 
	= fkPlaneSegmentation->GetSectorSegmentation(transformer->GetScale());
	  
      AliMpVPadIterator* sectorIt 
	= segmentation->CreateIterator(area);
	    
      fPadIterators.push_back(
        new AliMpTransformPadIterator(sectorIt, transformer));
    }	
  }
}

//
// public methods
//

//______________________________________________________________________________
void AliMpPlaneAreaPadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the area
// ---
  if (fPadIterators.size()==0) return;

  fCurrentIterator = fPadIterators.begin();
  (*fCurrentIterator)->First();

  while ( fCurrentIterator != fPadIterators.end() &&
          (*fCurrentIterator)->IsDone()) {
	 
    fCurrentIterator++;
    if (fCurrentIterator != fPadIterators.end()) {
      (*fCurrentIterator)->First();
    }  	 
  }
}

//______________________________________________________________________________
void AliMpPlaneAreaPadIterator::Next()
{
// Move the iterator to the next valid pad.
// ---

  (*fCurrentIterator)->Next();
  
  while ( fCurrentIterator != fPadIterators.end() &&
          (*fCurrentIterator)->IsDone()) {
	 
    fCurrentIterator++;
    if (fCurrentIterator != fPadIterators.end()) {
      (*fCurrentIterator)->First();
    }  	 
  }
}

//______________________________________________________________________________
Bool_t AliMpPlaneAreaPadIterator::IsDone() const
{
// 
  return  fCurrentIterator == fPadIterators.end();
}

//______________________________________________________________________________
AliMpPad AliMpPlaneAreaPadIterator::CurrentItem() const 
{
// Returns the current pad.
// ---

  if (fCurrentIterator != fPadIterators.end())
    return (*fCurrentIterator)->CurrentItem();
  else
    return AliMpPad::Invalid();  
}

//______________________________________________________________________________
void AliMpPlaneAreaPadIterator::Invalidate()
{
// Invalidates all sector iterators and sets the current
// iterator to invalid position.
// ---
 
  PadIteratorVectorIterator it;
  for (it=fPadIterators.begin(); it !=fPadIterators.end(); it++) {
    (*it)->Invalidate(); 
  }
  
  fCurrentIterator = fPadIterators.end();
}

