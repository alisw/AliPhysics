/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
// $MpId: AliMpGraphContext.cxx,v 1.9 2006/03/02 16:37:24 ivana Exp $
// Category: graphics
//
// Class AliMpGraphContext
// -----------------------
// Class describing a the correspondance between a given area
// in pad, and a zone of real (cm) position
// Included in AliRoot: 2003/05/02
// Author: David GUEZ, IPN Orsay

#include "AliMpGraphContext.h"

#include <TError.h>

AliMpGraphContext *AliMpGraphContext::fgInstance = 0;
AliMpGraphContext::GraphContextVector AliMpGraphContext::fgStack;
#ifdef WITH_ROOT
Int_t              AliMpGraphContext::fgStackSize = 0;
#endif

ClassImp(AliMpGraphContext)

//_____________________________________________________________________________
AliMpGraphContext::AliMpGraphContext():
  TObject(),
  fPadPosition(TVector2(0.5,0.5)),
  fPadDimensions(TVector2(0.49,0.49)),
  fRealPosition(TVector2(0.,0.)),
  fRealDimensions(TVector2(1,1))
{
/// Private constructor

  fColor = 20;
  // default constructor (private)
}

//_____________________________________________________________________________
AliMpGraphContext::AliMpGraphContext(const AliMpGraphContext& right) 
  : TObject(right),
    fColor(right.fColor),
    fPadPosition(right.fPadPosition),
    fPadDimensions(right.fPadDimensions),
    fRealPosition(right.fRealPosition),
    fRealDimensions(right.fRealDimensions)     
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMpGraphContext& 
AliMpGraphContext::operator=(const AliMpGraphContext& right)
{
/// Protected assignment operator

  // check assignment to self
  if (this == &right) return *this;

  fColor = right.fColor;
  fPadPosition = right.fPadPosition;
  fPadDimensions = right.fPadDimensions;
  fRealPosition = right.fRealPosition;
  fRealDimensions = right.fRealDimensions;
    
  return *this;  
}    

//_____________________________________________________________________________
AliMpGraphContext *AliMpGraphContext::Instance()
{
  /// Return or create a unique instance of this class
  
  if (fgInstance) return fgInstance;
  fgInstance = new AliMpGraphContext;
  return fgInstance;
}

//_____________________________________________________________________________
TVector2 AliMpGraphContext::RealToPad(const TVector2 &position) const
{
  /// Transform a real position into its equivalent position in the pad
  
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
  // Transform the real area (position,dimensions) to
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
  /// Set the pad area from the actual one
  /// corresponding to the given real area.

  RealToPad(position,dimensions,fPadPosition,fPadDimensions);
}

//_____________________________________________________________________________
void AliMpGraphContext::Push() const
{
  /// Store the current configuration

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
/// Pop an object from the stack.

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
