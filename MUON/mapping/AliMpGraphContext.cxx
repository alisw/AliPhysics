// $Id$
// Category: graphics
//
// Class AliMpGraphContext
// -----------------------
// Class describing a the correspondance between a given area
// in pad, and a zone of real (cm) position
// Included in AliRoot: 2003/05/02
// Author: David GUEZ, IPN Orsay

#include <TError.h>

#include "AliMpGraphContext.h"

ClassImp(AliMpGraphContext)

AliMpGraphContext *AliMpGraphContext::fgInstance = 0;
GraphContextVector AliMpGraphContext::fgStack;
#ifdef WITH_ROOT
Int_t              AliMpGraphContext::fgStackSize = 0;
#endif

//_____________________________________________________________________________
AliMpGraphContext::AliMpGraphContext():
  TObject(),
  fPadPosition(TVector2(0.5,0.5)),
  fPadDimensions(TVector2(0.49,0.49)),
  fRealPosition(TVector2(0.,0.)),
  fRealDimensions(TVector2(1,1))
{
// private constructor

  fColor = 20;
  // default constructor (private)
}

//_____________________________________________________________________________
AliMpGraphContext::AliMpGraphContext(const AliMpGraphContext& right) 
  : TObject(right) 
{
// protected copy constructor

  Fatal("AliMpGraphContext", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpGraphContext& 
AliMpGraphContext::operator=(const AliMpGraphContext& right)
{
// protected assignement operator

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//_____________________________________________________________________________
AliMpGraphContext *AliMpGraphContext::Instance()
{
  // return or create a unique instance of this class
  if (fgInstance) return fgInstance;
  fgInstance = new AliMpGraphContext;
  return fgInstance;
}

//_____________________________________________________________________________
TVector2 AliMpGraphContext::RealToPad(const TVector2 &position) const
{
  // transform a real position into its equivalent position in the pad
  Double_t x=position.X();
  Double_t y=position.Y();
  x-= (fRealPosition.X()-fRealDimensions.X());
  x/=fRealDimensions.X();
  x*=fPadDimensions.X();
  x+= (fPadPosition.X()-fPadDimensions.X() );

  y-= (fRealPosition.Y()-fRealDimensions.Y());
  y/=fRealDimensions.Y();
  y*=fPadDimensions.Y();
  y+= (fPadPosition.Y()-fPadDimensions.Y() );

  return TVector2(x,y);
}

//_____________________________________________________________________________
void AliMpGraphContext::RealToPad(const TVector2 &position,
			      const TVector2 &dimensions,
			      TVector2 &padPosition,
			      TVector2 &padDimensions) const
{
  // transform the real area (position,dimensions) to
  // its equivalent pad area
  padPosition = RealToPad(position);
  padDimensions = 
    TVector2(dimensions.X()*fPadDimensions.X()/fRealDimensions.X(),
	     dimensions.Y()*fPadDimensions.Y()/fRealDimensions.Y());

}

//_____________________________________________________________________________
void AliMpGraphContext::SetPadPosForReal(const TVector2 &position,
				     const TVector2 &dimensions)
{
  // Set the pad area from the actual one
  // corresponding to the given real area.
  RealToPad(position,dimensions,fPadPosition,fPadDimensions);
}

//_____________________________________________________________________________
void AliMpGraphContext::Push() const
{
  // Store the current configuration
  AliMpGraphContext *save = new AliMpGraphContext(*this);

#ifdef WITH_STL
  fgStack.push_back(save);
#endif

#ifdef WITH_ROOT
  fgStack.AddAt(save, fgStackSize++);
#endif
}

//_____________________________________________________________________________
void AliMpGraphContext::Pop()
{
// Pops object from the stack.
#ifdef WITH_STL
  // restore the last saved configuration
  if (!fgStack.empty()){
    AliMpGraphContext *obj = fgStack.back();
    *this = *obj;
    fgStack.pop_back();
    delete obj;
  }
#endif

#ifdef WITH_ROOT
  // restore the last saved configuration
  if ( fgStackSize ){
    AliMpGraphContext *obj 
      = (AliMpGraphContext*)fgStack.At(--fgStackSize);
    *this = *obj;
    fgStack.RemoveAt(fgStackSize);
    delete obj;
  }
#endif
}
