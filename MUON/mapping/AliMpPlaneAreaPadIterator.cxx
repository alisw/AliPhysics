// $Id$
// Category: plane
//
// Class AliMpPlaneAreaPadIterator
// -------------------------------
// Class, which defines an iterator over the pads 
// inside a given area in a plane in horizontal direction.
// Included in AliRoot: 2003/05/02
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

#ifdef WITH_STL
  fCurrentIterator = fPadIterators.end();
#endif

#ifdef WITH_ROOT
  fCurrentIterator = fPadIterators.GetEntriesFast();
#endif
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
	    
#ifdef WITH_STL
      fPadIterators.push_back(
        new AliMpTransformPadIterator(sectorIt, transformer));
#endif

#ifdef WITH_ROOT
      fPadIterators.Add(
        new AliMpTransformPadIterator(sectorIt, transformer));
#endif
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

#ifdef WITH_STL
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
#endif

#ifdef WITH_ROOT
  if (fPadIterators.GetEntriesFast()==0) return;

  fCurrentIterator = 0;
  ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))->First();

  while ( fCurrentIterator != fPadIterators.GetEntriesFast() &&
          ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))
	    ->IsDone()) {

    fCurrentIterator++;
    if (fCurrentIterator != fPadIterators.GetEntriesFast()) {
      ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))->First();
    }  	 
  }
#endif
}

//______________________________________________________________________________
void AliMpPlaneAreaPadIterator::Next()
{
// Move the iterator to the next valid pad.
// ---

#ifdef WITH_STL
  (*fCurrentIterator)->Next();
  
  while ( fCurrentIterator != fPadIterators.end() &&
          (*fCurrentIterator)->IsDone()) {
	 
    fCurrentIterator++;
    if (fCurrentIterator != fPadIterators.end()) {
      (*fCurrentIterator)->First();
    }  	 
  }
#endif

#ifdef WITH_ROOT
  ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))->Next();

  while ( fCurrentIterator != fPadIterators.GetEntriesFast() &&
          ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))
	    ->IsDone()) {
	 
    fCurrentIterator++;
    if (fCurrentIterator != fPadIterators.GetEntriesFast()) {
      ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))->First();
    }  	 
  }
#endif
}

//______________________________________________________________________________
Bool_t AliMpPlaneAreaPadIterator::IsDone() const
{
// 
#ifdef WITH_STL
  return  fCurrentIterator == fPadIterators.end();
#endif

#ifdef WITH_ROOT
  return  fCurrentIterator == fPadIterators.GetEntriesFast();
#endif
}

//______________________________________________________________________________
AliMpPad AliMpPlaneAreaPadIterator::CurrentItem() const 
{
// Returns the current pad.
// ---

#ifdef WITH_STL
  if (fCurrentIterator != fPadIterators.end())
    return (*fCurrentIterator)->CurrentItem();
  else
    return AliMpPad::Invalid();  
#endif

#ifdef WITH_ROOT
  if (fCurrentIterator != fPadIterators.GetEntriesFast())
    return ((AliMpTransformPadIterator*)fPadIterators.At(fCurrentIterator))
             ->CurrentItem();
  else
    return AliMpPad::Invalid();  
#endif
}

//______________________________________________________________________________
void AliMpPlaneAreaPadIterator::Invalidate()
{
// Invalidates all sector iterators and sets the current
// iterator to invalid position.
// ---
 
#ifdef WITH_STL
  PadIteratorVectorIterator it;
  for (it=fPadIterators.begin(); it !=fPadIterators.end(); it++) {
    (*it)->Invalidate(); 
  }
  
  fCurrentIterator = fPadIterators.end();
#endif

#ifdef WITH_ROOT
  for (Int_t i=0; i<fPadIterators.GetEntriesFast(); i++) {
    ((AliMpTransformPadIterator*)fPadIterators.At(i))->Invalidate();
  }  

  fCurrentIterator = fPadIterators.GetEntriesFast();
#endif
}

