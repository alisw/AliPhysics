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
// $MpId: AliMpMotifTypePadIterator.cxx,v 1.6 2006/05/24 13:58:41 ivana Exp $
// Category: motif

//-----------------------------------------------------------------------------
// Class AliMpMotifTypePadIterator
// -------------------------------
// Class, which defines an iterator over the pads of a given motif type
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpMotifTypePadIterator.h"
#include "AliMpMotifType.h"
#include "AliMpEncodePair.h"

/// \cond CLASSIMP
ClassImp(AliMpMotifTypePadIterator)
/// \endcond

//______________________________________________________________________________
AliMpMotifTypePadIterator::AliMpMotifTypePadIterator():
    AliMpVPadIterator(),
    fkMotifType(0),
    fCurrentIx(-1),
    fCurrentIy(-1)
{
/// Default constructor, set the current position to "invalid"
}

//______________________________________________________________________________
AliMpMotifTypePadIterator::AliMpMotifTypePadIterator( 
                                const AliMpMotifType* motifType)
  : AliMpVPadIterator(),
    fkMotifType(motifType),
    fCurrentIx(-1),
    fCurrentIy(-1)
{
/// Standard constructor, let *this to invalid position
}

//______________________________________________________________________________
AliMpMotifTypePadIterator::AliMpMotifTypePadIterator(
                                const AliMpMotifTypePadIterator& right)
  : AliMpVPadIterator(right),
    fkMotifType(right.fkMotifType),
    fCurrentIx(right.fCurrentIx),
    fCurrentIy(right.fCurrentIy)
    
{
/// Copy constructor
}

//______________________________________________________________________________
AliMpMotifTypePadIterator::~AliMpMotifTypePadIterator()
{
/// Destructor
}

// operators

//______________________________________________________________________________
AliMpMotifTypePadIterator& 
AliMpMotifTypePadIterator::operator = (const AliMpMotifTypePadIterator& right)
{
/// Assignment operator.                                                      \n
/// If the right hand iterator isn't of good type
/// the current operator is invalidated

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
  AliMpVPadIterator::operator=(right);

  fkMotifType = right.fkMotifType;
  fCurrentIx = right.fCurrentIx;
  fCurrentIy = right.fCurrentIy;

  return *this;
}  

//
//private methods
//

//______________________________________________________________________________
Bool_t
AliMpMotifTypePadIterator::FindFirstPadInLine(Int_t ix, Int_t iy, 
                                              Int_t& newIx, Int_t& newIy) const
{
/// Find the indices of the first pad in the same line
/// as the \a indices, and in column, at least equal, to the
/// one of \a indices

    if ( ! fkMotifType ) {
      newIx = -1;
      newIy = -1;
      return false;
    }  

    while ( ix < fkMotifType->GetNofPadsX() ) {
      if ( fkMotifType->HasPadByLocalIndices(ix, iy) ) {
        newIx = ix;
        newIy = iy;
        return true;
      }  
      ix++;
    }
    
    newIx = -1;
    newIy = -1;
    return false;
} 

//______________________________________________________________________________
Bool_t AliMpMotifTypePadIterator::IsValid() const
{
/// Is the iterator in a valid position?

    return fkMotifType!=0 && fCurrentIx >=0 &&  fCurrentIy >=0;
} 

//
//public methods
//

//______________________________________________________________________________
void AliMpMotifTypePadIterator::First()
{
/// Reset the iterator, so that it points to the first available
/// pad in the motif type

    if ( ! fkMotifType ) {
      Invalidate();
      return ;
    }
    
    fCurrentIx = 0;
    fCurrentIy = 0;
    if ( fkMotifType->HasPadByLocalIndices(fCurrentIx, fCurrentIy) ) return;
    
    
    // if (0,0) is not available
    // place itself to the first avalable motif after (0,0) (if exists)
    // ++(*this);
    Next();
    
    return;
}

//______________________________________________________________________________
void AliMpMotifTypePadIterator::Next()
{
/// Move the iterator to the next valid pad.

    if ( ! IsValid() ) return;

    while ( fCurrentIy < fkMotifType->GetNofPadsY() ) {
      Int_t nextTryIx, nextTryIy;
      Bool_t result 
        = FindFirstPadInLine(fCurrentIx+1, fCurrentIy, nextTryIx, nextTryIy);

      if ( result ){
        fCurrentIx = nextTryIx;
        fCurrentIy = nextTryIy;
        return;
      }
      fCurrentIx = -1;
      fCurrentIy++;
    }
    
    // if the loop is finished, there's not available pads at all...
    Invalidate();
    return;
}

//______________________________________________________________________________
Bool_t AliMpMotifTypePadIterator::IsDone() const
{
/// Is the iterator in the end ?
 
  return ! IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpMotifTypePadIterator::CurrentItem() const 
{
/// Return current pad.

    if ( ! fkMotifType )
      return AliMpPad::Invalid();
    else
      return AliMpPad(0, 0,
	              fCurrentIx, fCurrentIy,
                      TVector2(),TVector2());
}

//______________________________________________________________________________
void AliMpMotifTypePadIterator::Invalidate()
{
/// Let the iterator point to the invalid position

    fCurrentIx = -1;
    fCurrentIy = -1;
} 

