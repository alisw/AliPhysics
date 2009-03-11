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
// $MpId: AliMpNeighboursPadIterator.cxx,v 1.12 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpNeighboursPadIterator
// --------------------------------
// Class, which defines an iterator over the pads surrounding a given pad
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpNeighboursPadIterator.h"
#include "AliMpIntPair.h"
#include "AliMpVSegmentation.h"
#include "AliMpRow.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <TVector2.h>

/// \cond CLASSIMP
ClassImp(AliMpNeighboursPadIterator)
/// \endcond

const UInt_t AliMpNeighboursPadIterator::fgkInvalidIndex = 9999; 
                                                   //never so much neighbours...

//______________________________________________________________________________
AliMpNeighboursPadIterator::AliMpNeighboursPadIterator()
  : AliMpVPadIterator(),
    fkSegmentation(0),
    fCenterPad(AliMpPad::Invalid()),
    fPads(),
    fIndex(fgkInvalidIndex)
{
/// Default constructor, set the current position to "invalid"
}

//______________________________________________________________________________
AliMpNeighboursPadIterator::AliMpNeighboursPadIterator(
                                 const AliMpVSegmentation* segmentation,
                                 const AliMpPad& centerPad,
                                 Bool_t includeCenter)
  : AliMpVPadIterator(),
    fkSegmentation(segmentation),
    fCenterPad(centerPad),
    fPads(),
    fIndex(fgkInvalidIndex)
{
/// Standard constructor, set *this to invalid position

    FillPadsVector(includeCenter);
}

//______________________________________________________________________________
AliMpNeighboursPadIterator::AliMpNeighboursPadIterator(
                                 const AliMpNeighboursPadIterator& right)
  : AliMpVPadIterator(right),
    fkSegmentation(0),
    fCenterPad(AliMpPad::Invalid()),
    fPads(),
    fIndex(fgkInvalidIndex)
{
/// Copy constructor

  *this = right;
}

//______________________________________________________________________________
AliMpNeighboursPadIterator::~AliMpNeighboursPadIterator()
{
/// Destructor

  fPads.Delete();
}

// operators

//______________________________________________________________________________
AliMpNeighboursPadIterator& 
AliMpNeighboursPadIterator::operator = (const AliMpNeighboursPadIterator& right)
{
/// Assignment operator.                                                     \n
/// If the right hand iterator isn't of a good type
/// the current operator is invalidated                                      \n

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
  AliMpVPadIterator::operator=(right);

  AliFatal("Not allowed assignment for TObjArray");

  return *this;
} 

//
// private methods
//

//______________________________________________________________________________
Bool_t AliMpNeighboursPadIterator::IsNeighbour(const AliMpPad& pad) const
{
/// Return true if the pad located by \a pad is a neighbour of those
/// located at a\ fCenterPad

    if ( !pad.IsValid() ) return kFALSE;    

    TVector2 relPos  = pad.Position()   - fCenterPad.Position();
    TVector2 bounds  = pad.Dimensions() + fCenterPad.Dimensions();
    return (TMath::Abs(relPos.X())- bounds.X()<AliMpConstants::LengthTolerance()) && 
           (TMath::Abs(relPos.Y())- bounds.Y()<AliMpConstants::LengthTolerance());

}

//______________________________________________________________________________
TObjArray* 
AliMpNeighboursPadIterator::PadVectorLine(const AliMpPad& from,
                                          const AliMpIntPair& direction) const
{
/// Fill  a new vector with all pads which have common
/// parts with the pad located at \a fCenterPad, in a given line
/// starting from \a from and moving by \a direction

    AliMpPad current = from;
    TObjArray* ans = new TObjArray();
    Bool_t cont=kTRUE;
    do {
        if (IsNeighbour(current))
            ans->Add(new AliMpPad(current));
        else
            cont=kFALSE;
        TVector2 nextPos = current.Position() + TVector2(
          current.Dimensions().X()*(AliMpConstants::LengthStep()+1.)*direction.GetFirst(),
          current.Dimensions().Y()*(AliMpConstants::LengthStep()+1.)*direction.GetSecond());
        current = fkSegmentation->PadByPosition(nextPos, false);
    } while (cont);
    return ans;
}

//______________________________________________________________________________
void  AliMpNeighboursPadIterator::UpdateTotalSet(TObjArray& setTotal, 
                                                 TObjArray* from) const
{
/// Add pads from pad vector to the total set 
/// only if they are not yet included and deletes the pad vector

    for (Int_t i=0; i<from->GetEntriesFast(); i++) {
      AliMpPad* candidate = (AliMpPad*)from->At(i);
      
      Bool_t isInSetTotal = false;
      for (Int_t j=0; j<setTotal.GetEntriesFast(); j++) {
         AliMpPad* pad = (AliMpPad*)setTotal.At(j);
	 
	 if (pad->GetIndices() == candidate->GetIndices()) {
	   isInSetTotal = true;
	   break;
	 }       
      }
      if (!isInSetTotal) 
        setTotal.Add(candidate);
      else
        delete candidate;	
    }
    delete from;
} 

//______________________________________________________________________________
void AliMpNeighboursPadIterator::FillPadsVector(Bool_t includeCenter)
{
/// Fill the indices vector with all indices of pads which have common
/// parts with the pad located at \a fCenterPad

    if (!fkSegmentation || !fCenterPad.IsValid()) return;
    
    
    AliMpPad from;
    AliMpIntPair direction;
    TObjArray* found;
    
    // repare a unique simple associative container
    // --> no doublons, rapid insersion
    TObjArray setTotal;

  /////////////  Left side
  
  ////////////////// up direction
    
    from = fkSegmentation->PadsLeft(fCenterPad).GetFirst();
    direction = AliMpIntPair(0,1);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);

  ////////////////// down direction

    from = fkSegmentation->PadsDown(from).GetFirst(); // the Pad down is already added
    direction = AliMpIntPair(0,-1);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);
    
  /////////////  Up side
  
  ////////////////// right direction

    from = fkSegmentation->PadsUp(fCenterPad).GetFirst();
    direction = AliMpIntPair(1,0);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);
    
  ////////////////// left direction

    from = fkSegmentation->PadsLeft(from).GetFirst(); // the pad up is already added
    direction = AliMpIntPair(-1,0);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);
    
  /////////////  Right side
  
  ////////////////// Up direction
    
    from = fkSegmentation->PadsRight(fCenterPad).GetFirst();
    direction = AliMpIntPair(0,1);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);
    
  ////////////////// down direction

    from = fkSegmentation->PadsDown(from).GetFirst(); // the pad right is already added
    direction = AliMpIntPair(0,-1);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);
    
  /////////////  Down side
  
  ////////////////// Right direction

    from = fkSegmentation->PadsDown(fCenterPad).GetFirst();
    direction = AliMpIntPair(1,0);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);
    
  ////////////////// left direction
    
    from = fkSegmentation->PadsLeft(from).GetFirst(); // the pad down is already added
    direction = AliMpIntPair(-1,0);
    found = PadVectorLine(from,direction);
    UpdateTotalSet(setTotal, found);

    // fill the fIndices vector with the set (-->pass from a rapid insertion,
    // to rapid and indexed access, for the rest of the job)

    fPads.Delete();
    // include the center pad if requiered
    if (includeCenter) fPads.Add(new AliMpPad(fCenterPad));

    for (Int_t i = 0; i<setTotal.GetEntriesFast(); i++)
      fPads.Add(setTotal.At(i));
}

//______________________________________________________________________________
Bool_t AliMpNeighboursPadIterator::IsValid() const
{
/// Is the iterator in a valid position?

    return (fkSegmentation!=0 && fIndex!=fgkInvalidIndex);
} 

//public methods

//______________________________________________________________________________
void AliMpNeighboursPadIterator::First()
{
/// Reset the iterator, so that it points to the first available
/// pad in the sector

    if ((fkSegmentation != 0) && (fPads.GetEntriesFast() != 0)) 
      fIndex=0; 
    else 
      fIndex=fgkInvalidIndex;

}

//______________________________________________________________________________
void AliMpNeighboursPadIterator::Next()
{
/// Pre-increment operator. Should be used by default for iterating over
/// pads


  if (!IsValid()) return;
  
  if (Int_t(fIndex) < fPads.GetEntriesFast()-1) 
    fIndex++; 
  else 
    Invalidate();
}

//______________________________________________________________________________
Bool_t AliMpNeighboursPadIterator::IsDone() const
{
/// Is the iterator in the end?
 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpNeighboursPadIterator::CurrentItem() const 
{
/// Dereferencement function

  if (!IsValid())
    return AliMpPad::Invalid();
  else
    return *((AliMpPad*)fPads[fIndex]);
}

//______________________________________________________________________________
void AliMpNeighboursPadIterator::Invalidate()
{
/// Let the iterator point to the invalid position

    fIndex=fgkInvalidIndex;
}

