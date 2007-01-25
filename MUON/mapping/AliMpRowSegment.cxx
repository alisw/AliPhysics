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
// $MpId: AliMpRowSegment.cxx,v 1.10 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpRowSegment
// ---------------------
// Class describing a row segment composed of the 
// the identic motifs.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpRowSegment.h"
#include "AliMpRow.h"
#include "AliMpVMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifTypePadIterator.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <TMath.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpRowSegment)
/// \endcond

//_____________________________________________________________________________
AliMpRowSegment::AliMpRowSegment(AliMpRow* row, AliMpVMotif* motif, 
                                 AliMpIntPair padOffset, 
                                 Int_t nofMotifs,
                                 Int_t motifPositionId, Int_t motifPositionDId)
  : AliMpVRowSegment(),
    fNofMotifs(nofMotifs),
    fPadOffset(padOffset),
    fOffset(TVector2()),
    fRow(row),
    fMotif(motif),
    fMotifPositionId(motifPositionId),
    fMotifPositionDId(motifPositionDId)
{
/// Standard constructor
 
  // Keep pad offset in the low indices limits
  SetLowIndicesLimit(padOffset);
}

//_____________________________________________________________________________
AliMpRowSegment::AliMpRowSegment() 
  : AliMpVRowSegment(),
    fNofMotifs(0),
    fPadOffset(AliMpIntPair()),
    fOffset(TVector2()),
    fRow(0),
    fMotif(0),
    fMotifPositionId(0),
    fMotifPositionDId(0)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpRowSegment::~AliMpRowSegment() 
{
/// Destructor  
}

//
// private methods  
//

//_____________________________________________________________________________
Double_t AliMpRowSegment::FirstMotifCenterX() const
{
/// Return the x coordinate of the first motif center
/// in the global coordinate system.

  return fOffset.X();
}  

//_____________________________________________________________________________
Double_t AliMpRowSegment::LastMotifCenterX() const
{
/// Return the x coordinate of the last motif center
/// in the global coordinate system.

  return fOffset.X() + 2.*(fNofMotifs-1)*fMotif->Dimensions().X();
}

//_____________________________________________________________________________
Double_t AliMpRowSegment::MotifCenterX(Int_t motifPositionId) const
{
/// Return the x coordinate of the motif specified with
/// the given position identifier.

  // Check if x is in the row segment range
  if (! HasMotifPosition(motifPositionId)) {
    AliErrorStream() << "Outside row segment region" << endl;
    return 0;
  }
  
  // Find the position number in the segment  
  Int_t num = (motifPositionId - fMotifPositionId) *  fMotifPositionDId;

  return fOffset.X() + num*(fMotif->Dimensions().X() * 2.0);
}

//_____________________________________________________________________________
Double_t AliMpRowSegment::MotifCenterY(Int_t motifPositionId) const
{
/// Return the y coordinate of the motif specified with
/// the given position identifier.

  // Check if x is in the row segment range
  if (! HasMotifPosition(motifPositionId)) {
    AliErrorStream() << "Outside row segment region" << endl;
    return 0;
  }
  
  return GetRow()->Position().Y() + fOffset.Y();
}

//_____________________________________________________________________________
Bool_t AliMpRowSegment::IsInside(const TVector2& position, Bool_t warn) const
{
/// Check if the position is inside some motif of this row segment.

  Double_t minY = GetRow()->Position().Y() + fOffset.Y() - fMotif->Dimensions().Y();
  Double_t maxY = GetRow()->Position().Y() + fOffset.Y() + fMotif->Dimensions().Y();

  if ( position.X() < LeftBorderX() || position.X() > RightBorderX() ||
       position.Y() < minY || position.Y() > maxY ) {

    if (warn)
      AliWarningStream() << "Outside row segment region" << endl;
    return false;
  }
  else
    return true;
}    

//
// public methods  
//

//_____________________________________________________________________________
Double_t  AliMpRowSegment::LeftBorderX() const
{
/// Return the x coordinate of the left row segment border
/// in the global coordinate system.

  return FirstMotifCenterX() - fMotif->Dimensions().X();
}

//_____________________________________________________________________________
Double_t  AliMpRowSegment::RightBorderX() const
{
/// Return the x coordinate of the right row segment border
/// in the global coordinate system.

  return LastMotifCenterX() + fMotif->Dimensions().X();
}

//_____________________________________________________________________________
Double_t  AliMpRowSegment::HalfSizeY() const
{
/// Return the size in y of this row segment.

  return fMotif->Dimensions().Y() + fOffset.Y();
}

//_____________________________________________________________________________
AliMpVMotif*  AliMpRowSegment::FindMotif(const TVector2& position) const
{
/// Return the motif of this row; 

  if (IsInside(position, false))
    return fMotif;
  else  
    return 0;
}  

//_____________________________________________________________________________
Int_t AliMpRowSegment::FindMotifPositionId(const TVector2& position) const
{
/// Return the motif position identified for the given
/// geometric position.

  if (!IsInside(position, false)) return 0;

  // Find the position number in the segment  
  Int_t num 
    = Int_t((position.X() - LeftBorderX()) / (fMotif->Dimensions().X() * 2.0));

  // Calculate the position Id
  return fMotifPositionId + num*fMotifPositionDId;  
}

//_____________________________________________________________________________
Bool_t AliMpRowSegment::HasMotifPosition(Int_t motifPositionId) const
{
/// Return true if the motif specified with the given position identifier
/// is in this segment.

  Int_t minId = TMath::Min(fMotifPositionId, 
                    fMotifPositionId + (fNofMotifs-1)*fMotifPositionDId);
  Int_t maxId = TMath::Max(fMotifPositionId, 
                    fMotifPositionId + (fNofMotifs-1)*fMotifPositionDId);

  if (motifPositionId >= minId && motifPositionId <= maxId) {
    return true;
  }
  else 
    return false;
}

//_____________________________________________________________________________
TVector2 AliMpRowSegment::MotifCenter(Int_t motifPositionId) const
{
/// Return the coordinates of the motif specified with
/// the given position identifier.

  return TVector2(MotifCenterX(motifPositionId), MotifCenterY(motifPositionId));
}

//_____________________________________________________________________________
TVector2 AliMpRowSegment::Position() const
{
/// Return the position of the row segment centre.

  Double_t x = (LeftBorderX() + RightBorderX())/2.;		    
  Double_t y = GetRow()->Position().Y();  
    
  return TVector2(x, y);   
}


//_____________________________________________________________________________
TVector2 AliMpRowSegment::Dimensions() const
{
/// Return the halflengths of the row segment in x, y.
// ---

  Double_t x = (RightBorderX() - LeftBorderX())/2.;		    
  Double_t y = GetRow()->Dimensions().Y();  
    
  return TVector2(x, y);   
}

//_____________________________________________________________________________
void   AliMpRowSegment::SetOffset(const TVector2& offset)
{
/// Calculate offset from given offset and 
/// stored offset in pads.

  AliMpMotifTypePadIterator iter(fMotif->GetMotifType());
  iter.First();
  AliMpIntPair localPos = iter.CurrentItem().GetIndices();
     
  Double_t offsetX 
     = offset.X() 
       + 2.*fPadOffset.GetFirst() * fMotif->GetPadDimensions(localPos).X() 
       + fMotif->Dimensions().X(); 

  Double_t offsetY 
    = offset.Y()
      + fPadOffset.GetSecond() * fMotif->GetPadDimensions(localPos).Y(); 

  fOffset = TVector2(offsetX, offsetY);
}

#include <Riostream.h>
//_____________________________________________________________________________
void AliMpRowSegment::SetGlobalIndices(AliMpRow* /*rowBefore*/)
{
/// Set global indices limits.

  // The low/high indices limits has to be taken as the highest/lowest from all 
  // motif positions
  Int_t ixl = 9999;
  Int_t iyl = 9999;
  Int_t ixh = AliMpConstants::StartPadIndex();
  Int_t iyh = AliMpConstants::StartPadIndex();

  for (Int_t i=0; i<GetNofMotifs(); i++) {
     
     AliMpMotifPosition* mPos 
       = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(i));
       
     // Check if the motif positions has the limits set
     if ( !mPos->HasValidIndices() )
       Fatal("SetGlobalIndices", 
             "Indices of motif positions have to be set first.");
            
     if ( mPos->GetLowIndicesLimit().GetFirst() < ixl ) 
       ixl = mPos->GetLowIndicesLimit().GetFirst();
       
     if ( mPos->GetLowIndicesLimit().GetSecond() < iyl ) 
       iyl = mPos->GetLowIndicesLimit().GetSecond();

     if ( mPos->GetHighIndicesLimit().GetFirst() > ixh ) 
       ixh = mPos->GetHighIndicesLimit().GetFirst();
       
     if ( mPos->GetHighIndicesLimit().GetSecond() > iyh ) 
       iyh = mPos->GetHighIndicesLimit().GetSecond();
  }     

  SetLowIndicesLimit(AliMpIntPair(ixl, iyl));
  SetHighIndicesLimit(AliMpIntPair(ixh, iyh));
}  

//_____________________________________________________________________________
Int_t AliMpRowSegment::SetIndicesToMotifPosition(Int_t i, 
                                 const AliMpIntPair& indices)
{
/// Set global indices to i-th motif position and returns next index
/// in x.

  // Get motif position
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(i));

  // Low limit
  AliMpIntPair low = indices + AliMpIntPair(0, GetLowIndicesLimit().GetSecond());
  motifPosition->SetLowIndicesLimit(low);
	  
  // High limit
  AliMpMotifType* motifType = motifPosition->GetMotif()->GetMotifType();    
  AliMpIntPair high
    = motifPosition->GetLowIndicesLimit() 
      + AliMpIntPair(motifType->GetNofPadsX()-1, motifType->GetNofPadsY()-1);
  motifPosition->SetHighIndicesLimit(high);

  // Return next index in x
  return high.GetFirst()+1;
  // return motifType->GetNofPadsX();
}


//_____________________________________________________________________________
AliMpRow*  AliMpRowSegment::GetRow() const
{
/// Return the row.which this row segment belongs to.

  return fRow;
}  

//_____________________________________________________________________________
Int_t  AliMpRowSegment::GetNofMotifs() const
{
/// Return number of motifs in this this row segment.

  return fNofMotifs;
}  

//_____________________________________________________________________________
Int_t  AliMpRowSegment::GetMotifPositionId(Int_t i) const
{
/// Return number of motifs in this this row segment.

  return fMotifPositionId + i*fMotifPositionDId;
}  

//_____________________________________________________________________________
AliMpVMotif*  AliMpRowSegment::GetMotif(Int_t /*i*/) const
{
/// Return the motif of this row segment.

  return fMotif;
}  
