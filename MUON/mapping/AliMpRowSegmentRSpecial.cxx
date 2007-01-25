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
// $MpId: AliMpRowSegmentRSpecial.cxx,v 1.7 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpRowSegmentRSpecial
// -----------------------------
// Class describing a special outer row segment composed of the 
// pad rows.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpRowSegmentRSpecial.h"
#include "AliMpRow.h"
#include "AliMpPadRow.h"
#include "AliMpVPadRowSegment.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMpRowSegmentRSpecial)
/// \endcond

//______________________________________________________________________________
AliMpRowSegmentRSpecial::AliMpRowSegmentRSpecial(AliMpRow* row, Double_t offsetX)
  : AliMpVRowSegmentSpecial(row, offsetX)
{
/// Standard constructor 
}

//______________________________________________________________________________
AliMpRowSegmentRSpecial::AliMpRowSegmentRSpecial() 
  : AliMpVRowSegmentSpecial()
{
/// Default constructor 
}

//______________________________________________________________________________
AliMpRowSegmentRSpecial::~AliMpRowSegmentRSpecial() 
{
/// Destructor  
}

//
// private methods  
//

//______________________________________________________________________________
AliMpVPadRowSegment*  
AliMpRowSegmentRSpecial::FindMostLeftPadRowSegment(Int_t motifPositionId) const
{
/// Find the most left pad row segment with this motifPositionId.

  AliMpVPadRowSegment* found = 0;

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if ( padRowSegment->GetMotifPositionId() == motifPositionId &&
           (!found || padRowSegment->LeftBorderX() < found->LeftBorderX()))
	   
        found = padRowSegment;  
    }
  }

  return found;   	
}

//______________________________________________________________________________
void AliMpRowSegmentRSpecial::SetGlobalIndicesLow()
{
/// Set global low indices
   
  // Last normal row segment in the row
  // (preceding this special row segment)
  AliMpVRowSegment* rowSegment 
    =  GetRow()->GetRowSegment(GetRow()->GetNofRowSegments()-2);
    
  // Set low indices limit to continue indices of the
  // preceding row segment  
  Int_t ix = rowSegment->GetHighIndicesLimit().GetFirst() + 1;
  Int_t iy = rowSegment->GetLowIndicesLimit().GetSecond();
  
  SetLowIndicesLimit(AliMpIntPair(ix, iy));
}

//
// protected methods  
//

//______________________________________________________________________________
TVector2 AliMpRowSegmentRSpecial::MotifCenterSlow(Int_t motifPositionId) const
{
/// Return the coordinates of the motif specified with
/// the given position identifier.                                           \n
/// !! Applicable only for motifs that have their most down pad in
/// this row segment.

  // Find the first (left, down) pad row segment with this motifPositionId.
  AliMpVPadRowSegment* downPadRowSegment 
    = FindPadRowSegment(motifPositionId);
  AliMpVPadRowSegment* leftPadRowSegment 
    = FindMostLeftPadRowSegment(motifPositionId);
  
  // Check if the motifPositionId is present 
  if (!downPadRowSegment || !leftPadRowSegment) {
    AliErrorStream() << "Outside row segment region" << endl;
    return 0;
  }

  // Check if both pad row segments have the same motif 
  if (downPadRowSegment->GetMotif() != leftPadRowSegment->GetMotif()) {
    AliFatal("Outside row segment region");
    return 0;
  }

  // Get position of found row segment
  Double_t x = leftPadRowSegment->LeftBorderX();       
  Double_t y = GetRow()->LowBorderY()  ;   
  
  for (Int_t i=0; i<downPadRowSegment->GetPadRow()->GetID(); i++)
    y += GetPadRow(i)->HalfSizeY()*2.;
    
  // Add motifs dimensions
  x += downPadRowSegment->GetMotif()->Dimensions().X();
  y += downPadRowSegment->GetMotif()->Dimensions().Y();
  
  return TVector2(x, y);
}

//
// public methods  
//

//______________________________________________________________________________
Double_t  AliMpRowSegmentRSpecial::LeftBorderX() const
{
/// Return the x coordinate of the left row segment border
/// in the global coordinate system.

  // The right edge of the last normal segment
  Double_t sameBorder = GetOffsetX();

  // Consistence check  
  Double_t leftBorder = DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border = padRow->GetPadRowSegment(0)->LeftBorderX();
    if (border < leftBorder) leftBorder =  border;
  }  

  if (TMath::Abs(sameBorder - leftBorder) > 1.e-04)  {
    AliErrorStream() << "WrongBorder" << endl;;
    return sameBorder;
  }  
  
  return leftBorder;

}

//______________________________________________________________________________
Double_t  AliMpRowSegmentRSpecial::RightBorderX() const
{
/// Return the x coordinate of the right row segment border
/// in the global coordinate system.

  Double_t rightBorder = -DBL_MAX;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
    Double_t border 
      = padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)
        ->RightBorderX();
      
    if (border > rightBorder) rightBorder =  border;
  }  
  
  return rightBorder;
}

//______________________________________________________________________________
TVector2 AliMpRowSegmentRSpecial::Position() const
{
/// Return the position of the row segment centre.
/// The centre is defined as the centre of the rectangular
/// row segment envelope.

  // The right edge of the last normal segment
  Double_t x = GetOffsetX()  + Dimensions().X();
  Double_t y = GetRow()->Position().Y();  
    
  return TVector2(x, y);   
}

//______________________________________________________________________________
Int_t AliMpRowSegmentRSpecial::SetIndicesToMotifPosition(Int_t i, 
                                         const AliMpIntPair& indices)
{
/// Set global indices to i-th motif position and returns next index in x.

  // Update low indices limit for this row segment
  SetGlobalIndicesLow();

  // Check for consistence
  if (GetLowIndicesLimit().GetFirst() != indices.GetFirst()) 
    AliFatal("Inconsistent indices");

  // Get motif position
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(i));
    
  // Set limits only once
  if ( motifPosition->GetHighIndicesLimit().IsValid() ) 
    return indices.GetFirst();; 

  // Low limit
  //
  Int_t ixl = GetLowIndicesLimit().GetFirst();
  Int_t iyl = GetLowIndicesLimit().GetSecond();

  // Find the most down pad row segment with this motifPositionId.
  AliMpVPadRowSegment* padRowSegment = FindPadRowSegment(GetMotifPositionId(i));
  Int_t padRowID = padRowSegment->GetPadRow()->GetID();
  iyl += padRowID; 

  // Add pads offset of this motif position in the row segment
  for (Int_t im=0; im<i; im++) {
    AliMpVPadRowSegment* rs = GetPadRow(padRowID)->GetPadRowSegment(im);
    if ( rs->GetMotifPositionId() == GetMotifPositionId(i) ) break; 
    ixl += rs->GetNofPads();
  }  
  motifPosition->SetLowIndicesLimit(AliMpIntPair(ixl, iyl));

  // High limit	
  //     
  AliMpMotifType* motifType = motifPosition->GetMotif()->GetMotifType();  
  AliMpIntPair high 
    = motifPosition->GetLowIndicesLimit()
      + AliMpIntPair(motifType->GetNofPadsX()-1, motifType->GetNofPadsY()-1);            
  motifPosition->SetHighIndicesLimit(high);

  // No increment index needed (this is always the last element)
  return indices.GetFirst();
}

//______________________________________________________________________________
void AliMpRowSegmentRSpecial::SetGlobalIndices(AliMpRow* rowBefore)
{
/// Set indices limits.
/// The limits are defined as the limits of the smallest rectangle which
/// includes all pads of this special row segment.

  // Get first motif position
  AliMpMotifPosition* firstMotifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(GetMotifPositionId(0));
    
  // Low ix
  Int_t ixl = firstMotifPosition->GetLowIndicesLimit().GetFirst();
              // We have to take the motif position limit
	      // as it can overlap over more rows and the indices
	      // of the right border of the precedent normal segment
	      // differ from one row to another  

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


