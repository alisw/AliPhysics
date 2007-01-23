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
// $MpId: AliMpRowSegmentLSpecial.cxx,v 1.7 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpRowSegmentLSpecial
// -----------------------------
// Class describing a special inner row segment composed of the 
// pad rows.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>
#include <TMath.h>

#include "AliMpRowSegmentLSpecial.h"
#include "AliMpRow.h"
#include "AliMpPadRow.h"
#include "AliMpVPadRowSegment.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpConstants.h"

#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMpRowSegmentLSpecial)
/// \endcond

//______________________________________________________________________________
AliMpRowSegmentLSpecial::AliMpRowSegmentLSpecial(AliMpRow* row, Double_t offsetX)
  : AliMpVRowSegmentSpecial(row, offsetX)
{
/// Standard constructor
}

//______________________________________________________________________________
AliMpRowSegmentLSpecial::AliMpRowSegmentLSpecial() 
  : AliMpVRowSegmentSpecial()
{
/// Default constructor
}

//______________________________________________________________________________
AliMpRowSegmentLSpecial::~AliMpRowSegmentLSpecial() 
{
/// Destructor  
}

//
// private methods  
//

//______________________________________________________________________________
AliMpVPadRowSegment*  
AliMpRowSegmentLSpecial::FindMostRightPadRowSegment(Int_t motifPositionId) const
{
/// Find the most right pad row segment with this motifPositionId.

  AliMpVPadRowSegment* found = 0;

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if ( padRowSegment->GetMotifPositionId() == motifPositionId &&
           (!found || padRowSegment->RightBorderX() > found->RightBorderX()))
	   
        found = padRowSegment;  
    }
  }

  return found;   	
}

//
// protected methods  
//

//______________________________________________________________________________
TVector2 AliMpRowSegmentLSpecial::MotifCenterSlow(Int_t motifPositionId) const
{
/// Return the coordinates of the motif specified with
/// the given position identifier.                                           \n
/// !! Applicable only for motifs that have their most down pad in
/// this row segment.

  // Find the first (left, down) pad row segment with this motifPositionId.
  AliMpVPadRowSegment* downPadRowSegment 
    = FindPadRowSegment(motifPositionId);
  AliMpVPadRowSegment* rightPadRowSegment 
    = FindMostRightPadRowSegment(motifPositionId);
  
  // Check if the motifPositionId is present 
  if (!downPadRowSegment || !rightPadRowSegment) {
    AliErrorStream() << "Outside row segment region" << endl;
    return 0;
  }

  // Check if both pad row segments have the same motif 
  if (downPadRowSegment->GetMotif() != rightPadRowSegment->GetMotif()) {
    AliFatal("Outside row segment region");
    return 0;
  }

  // Get position of found row segment
  Double_t x = rightPadRowSegment->RightBorderX();       
  Double_t y = GetRow()->LowBorderY()  ;   
  
  for (Int_t i=0; i<downPadRowSegment->GetPadRow()->GetID(); i++)
    y += GetPadRow(i)->HalfSizeY()*2.;
    
  // Add motifs dimensions
  x -= downPadRowSegment->GetMotif()->Dimensions().X();
  y += downPadRowSegment->GetMotif()->Dimensions().Y();
  
  return TVector2(x, y);
}

//
// public methods  
//

//______________________________________________________________________________
void AliMpRowSegmentLSpecial::UpdatePadsOffset()
{
/// Set low indices limit to the pad offset calculated
/// from the neighbour normal segment.

  // Get the neighbour row segment
  // (the first normal segment)
  AliMpVRowSegment* neighbour = GetRow()->GetRowSegment(1);

  // Get the the pads offset of the neighbour row segment
  // (the first normal segment)
  AliMpIntPair offset = neighbour->GetLowIndicesLimit();
  
  // Find max nof pads in a row
  Int_t maxNofPads = MaxNofPadsInRow();

  // Set limits
  SetLowIndicesLimit(offset - AliMpIntPair(maxNofPads, 0));

  // Reset limits in the neighbour row segment
  // (pad offset is now included in the special segment)  
  neighbour->SetLowIndicesLimit(
    AliMpIntPair(0, neighbour->GetLowIndicesLimit().GetSecond()));
}

//______________________________________________________________________________
Double_t  AliMpRowSegmentLSpecial::LeftBorderX() const
{
/// Return the x coordinate of the left row segment border
/// in the global coordinate system.

  Double_t leftBorder = DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border 
      = padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)->LeftBorderX();
      
    if (border < leftBorder) leftBorder =  border;
  }  
  
  return leftBorder;
}

//______________________________________________________________________________
Double_t  AliMpRowSegmentLSpecial::RightBorderX() const
{
/// Returns the x coordinate of the right row segment border
/// in the global coordinate system.

  Double_t sameBorder = GetOffsetX();

  // Consistence check  
  Double_t rightBorder = -DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border = padRow->GetPadRowSegment(0)->RightBorderX();
    if (border > rightBorder) rightBorder =  border;
  }  

  if (TMath::Abs(GetOffsetX() - rightBorder) > 1.e-04)  {
    AliErrorStream() << "WrongBorder" << endl;
    return sameBorder;
  }  
  
  return rightBorder;
}

//______________________________________________________________________________
TVector2 AliMpRowSegmentLSpecial::Position() const
{
/// Return the position of the row segment centre.
/// The centre is defined as the centre of the rectangular
/// row segment envelope.

  Double_t x = GetOffsetX() - Dimensions().X();		    
  Double_t y = GetRow()->Position().Y();  
    
  return TVector2(x, y);   
}

#include <Riostream.h>
//______________________________________________________________________________
Int_t AliMpRowSegmentLSpecial::SetIndicesToMotifPosition(Int_t i, 
                                         const AliMpIntPair& indices)
{
/// Set global indices to i-th motif position and returns next index in x.

  // Get motif position
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(i));

  // Low limit
  AliMpIntPair low 
    = AliMpIntPair(GetLowIndicesLimit().GetFirst() + AliMpConstants::StartPadIndex(), 
	       indices.GetSecond())
      + FindRelativeLowIndicesOf(GetMotifPositionId(i));
            
  if (! motifPosition->GetHighIndicesLimit().IsValid()) {   
     motifPosition->SetLowIndicesLimit(low);
  } 
  else {
    if (motifPosition->GetLowIndicesLimit().GetFirst() > low.GetFirst())
      motifPosition->SetLowIndicesLimit(
                        AliMpIntPair(low.GetFirst(),
                                 motifPosition->GetLowIndicesLimit().GetSecond()));

    if (motifPosition->GetLowIndicesLimit().GetSecond() > low.GetSecond())
       motifPosition->SetLowIndicesLimit(
                         AliMpIntPair(motifPosition->GetLowIndicesLimit().GetFirst(),
                                  low.GetSecond()));
  }
             
  // High limit	     
  AliMpMotifType* motifType = motifPosition->GetMotif()->GetMotifType();  
  AliMpIntPair high 
    = motifPosition->GetLowIndicesLimit()
      + AliMpIntPair(motifType->GetNofPadsX()-1, motifType->GetNofPadsY()-1);  
                
  motifPosition->SetHighIndicesLimit(high);

  // Increment index only if last motif position is processed 
  if (i != GetNofMotifs()-1) 
    return indices.GetFirst();
    //return 0;
  else
    return indices.GetFirst() + MaxNofPadsInRow();  
    //return MaxNofPadsInRow();  
}
//______________________________________________________________________________
void AliMpRowSegmentLSpecial::SetGlobalIndices(AliMpRow* rowBefore)
{
/// Set indices limits
/// The limits are defined as the limits of the smallest rectangle which
/// includes all pads of this special row segment.

  // Low ix
  Int_t ixl = GetLowIndicesLimit().GetFirst() + AliMpConstants::StartPadIndex();
      // the pads offset was already defined by Reader

  // High ix
  Int_t ixh = ixl + MaxNofPadsInRow() - 1;

  // Low iy
  Int_t iyl = AliMpConstants::StartPadIndex();
  if (rowBefore) {
    //if (constPadSizeDirection == kY) {
      iyl = rowBefore->GetHighIndicesLimit().GetSecond()+1;
    //} 
    /*
    else {
      AliMpVRowSegment* seg = rowBefore->FindRowSegment(ixl);	
      AliMpMotifPosition* motPos =  rowBefore->FindMotifPosition(seg, ixl);
      if (!motPos) 
        Fatal("SetGlobalIndices", "Motif position in rowBefore not found.");
      iyl = motPos->GetHighIndicesLimit().GetSecond()+1;
    }
    */
  }  

  // High iy
  Int_t iyh = iyl + GetNofPadRows() - 1;
  
  SetLowIndicesLimit(AliMpIntPair(ixl, iyl));
  SetHighIndicesLimit(AliMpIntPair(ixh, iyh));
}  




