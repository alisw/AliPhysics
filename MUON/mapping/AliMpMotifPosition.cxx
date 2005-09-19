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
// $MpId: AliMpMotifPosition.cxx,v 1.7 2005/08/26 15:43:36 ivana Exp $
//
// Class AliMpMotifPosition
// ------------------------
// Class that represents a placed motif.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpMotifPosition.h"
#include "AliMpMotifPositionPadIterator.h"
#include "AliMpMotifType.h"
#include <iostream>

ClassImp(AliMpMotifPosition)

//______________________________________________________________________________
AliMpMotifPosition::AliMpMotifPosition(Int_t id, AliMpVMotif* motif, 
                                       TVector2 position)
  : AliMpVIndexed(),
    fID(id),
    fMotif(motif),
    fPosition(position) 
{
/// Standard constructor
}

//______________________________________________________________________________
AliMpMotifPosition::AliMpMotifPosition()
  : AliMpVIndexed(), 
    fID(0),
    fMotif(0),
    fPosition(TVector2(0.,0.)) 
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpMotifPosition::AliMpMotifPosition(const AliMpMotifPosition& right) 
  : AliMpVIndexed(right) 
{
/// Protected copy constructor (not provided)

  Fatal("AliMpMotifPosition", "Copy constructor not provided.");
}

//______________________________________________________________________________
AliMpMotifPosition::~AliMpMotifPosition()\
{
/// Destructor 
}

// operators

//_____________________________________________________________________________
AliMpMotifPosition& 
AliMpMotifPosition::operator=(const AliMpMotifPosition& right)
{
/// Protected assignment operator (not provided)

  // check assignment to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignment operator not provided.");
    
  return *this;  
}    

//______________________________________________________________________________
AliMpVPadIterator* AliMpMotifPosition::CreateIterator() const
{
/// Return motif position iterator

  return new AliMpMotifPositionPadIterator(this);
}  

//______________________________________________________________________________
Bool_t AliMpMotifPosition::HasPad(const AliMpIntPair& indices) const
{
/// Return true if pad with the specified indices exists in 
/// this motif position.

  if (!HasIndices(indices)) return kFALSE;
  
  if (fMotif->GetMotifType()->IsFull()) return kTRUE;
  
  return fMotif->GetMotifType()->HasPad(indices-GetLowIndicesLimit());
}

//_____________________________________________________________________________
void
AliMpMotifPosition::SetID(Int_t id)
{
/// Set ID

  fID = id;
}

//_____________________________________________________________________________
void
AliMpMotifPosition::SetPosition(const TVector2& pos)
{
/// Set position

  fPosition = pos;
}

//_____________________________________________________________________________
void
AliMpMotifPosition::Print(Option_t* option) const
{
/// Printing

  std::cout << "MOTIFPOSITION " << GetID() << " MOTIF " 
	    << GetMotif()->GetID()
	    << " at (" << Position().X() << "," 
	    << Position().Y() << ") "
	    << " iMin=(" << GetLowIndicesLimit().GetFirst()
	    << "," << GetLowIndicesLimit().GetSecond()
	    << ") iMax=(" << GetHighIndicesLimit().GetFirst()
	    << "," << GetHighIndicesLimit().GetSecond()
	    << ")" << std::endl;

  if ( option && option[0] == 'M' )
    {
      GetMotif()->Print(option+1);
    }
}
