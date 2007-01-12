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
// $MpId: AliMpMotifPosition.cxx,v 1.9 2006/05/24 13:58:41 ivana Exp $
//
// Class AliMpMotifPosition
// ------------------------
// Class that represents a placed motif.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotifPosition.h"
#include "AliMpMotifPositionPadIterator.h"
#include "AliMpMotifType.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpMotifPosition)
/// \endcond

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

//______________________________________________________________________________
AliMpMotifPosition::~AliMpMotifPosition()
{
/// Destructor 
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

  cout << "MOTIFPOSITION " << GetID() << " MOTIF " 
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
