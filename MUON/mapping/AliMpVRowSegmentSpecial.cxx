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
// $MpId: AliMpVRowSegmentSpecial.cxx,v 1.9 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpVRowSegmentSpecial
// ----------------------------
// Class describing a special row segment composed of the 
// pad rows.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpVRowSegmentSpecial.h"
#include "AliMpRow.h"
#include "AliMpPadRow.h"
#include "AliMpVPadRowSegment.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpConstants.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpVRowSegmentSpecial)
/// \endcond

const Int_t  AliMpVRowSegmentSpecial::fgkMaxNofMotifPositionIds = 20;

//______________________________________________________________________________
AliMpVRowSegmentSpecial::AliMpVRowSegmentSpecial(AliMpRow* row, Double_t offsetX)
  : AliMpVRowSegment(),
    fRow(row),
    fOffsetX(offsetX),
    fPadRows(),
    fMotifs(),
    fMotifPositionIds(),
    fNofMotifPositionIds(0)
{
/// Standard constructor  
}

//______________________________________________________________________________
AliMpVRowSegmentSpecial::AliMpVRowSegmentSpecial() 
  : AliMpVRowSegment(),
    fRow(0),
    fOffsetX(0.),
    fPadRows(),
    fMotifs(),
    fMotifPositionIds(),
    fNofMotifPositionIds(0)
{
/// Default constructor  

   fMotifPositionIds.Set(fgkMaxNofMotifPositionIds);
}

//______________________________________________________________________________
AliMpVRowSegmentSpecial::~AliMpVRowSegmentSpecial() 
{
/// Destructor
  
  for (Int_t i=0; i<GetNofPadRows(); i++)
    delete fPadRows[i];
}

//
// protected methods  
//

//______________________________________________________________________________
AliMpPadRow*  AliMpVRowSegmentSpecial::FindPadRow(Double_t y) const
{
/// Find the pad row in the given y coordinate.

  Double_t lowBorder =  fRow->LowBorderY();
  Double_t highBorder = fRow->LowBorderY();
  
  for (Int_t i=0; i<GetNofPadRows(); i++) {    

    AliMpPadRow* padRow = GetPadRow(i);
    highBorder += 2.*padRow->HalfSizeY();

    if ( y >= lowBorder &&  y <= highBorder)
      return padRow;

    lowBorder = highBorder;
  }
  
  return 0;   	
}

//______________________________________________________________________________
AliMpVPadRowSegment*  
AliMpVRowSegmentSpecial::FindPadRowSegment(Int_t motifPositionId) const
{
/// Find the most down pad row segment with this motifPositionId.

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);    

    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) { 
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);

      if (padRowSegment->GetMotifPositionId() == motifPositionId) 
        return padRowSegment;
    }
  }
  return 0;   	
}

//______________________________________________________________________________
MpPair_t 
AliMpVRowSegmentSpecial::FindRelativeLowIndicesOf(Int_t motifPositionId) const 
{ 
/// Return the lowest pad indices where the motif of the given position ID
/// exist in this segment.

  Int_t ans0Ix = 0;
  Int_t ans0Iy = 1000;
  Int_t ansIx = 0;
  Int_t ansIy = 1000;
  
  Int_t maxNofPadsX=0;
  
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);

    Int_t nofPadsX=0;
    for ( Int_t j=0; j<padRow->GetNofPadRowSegments(); j++ ) {
      AliMpVPadRowSegment* padRowSegment = padRow->GetPadRowSegment(j);
      nofPadsX += padRowSegment->GetNofPads();
      if (padRowSegment->GetMotifPositionId() == motifPositionId) {
         if ( ansIx < nofPadsX ) ansIx = nofPadsX;
         if ( ansIy >i ) ansIy = i;
                  // ansIx = max (nof pads of this pos ID)
                  // ansIy = min of pad row number
      }
    }  
    if ( nofPadsX > maxNofPadsX ) maxNofPadsX = nofPadsX;
  }    
  if ( ansIx == ans0Ix && ansIy == ans0Iy ) return -1;
  
  return AliMp::Pair(maxNofPadsX-ansIx, ansIy);
}
 
//______________________________________________________________________________
Int_t  AliMpVRowSegmentSpecial::MaxNofPadsInRow() const 
{ 
/// Return the maximum number of pads in this row segment along the X direction

  Int_t maxNofPads = 0;    

  for (Int_t i=0; i<GetNofPadRows(); i++){
    Int_t nofPads = GetPadRow(i)->GetNofPads(); 

    // Find maximum
    if (nofPads > maxNofPads) maxNofPads = nofPads;
  }
    
  return maxNofPads;
}

//______________________________________________________________________________
Bool_t AliMpVRowSegmentSpecial::HasMotif(const AliMpVMotif* motif) const
{
/// Return true if the specified motif is already in fMotifs vector,
/// returns false otherwise.

  for (Int_t i=0; i<fMotifs.GetEntriesFast(); i++)
    if (fMotifs[i] == (const TObject*)motif) return true;

  return false;	 
}

//______________________________________________________________________________
Int_t AliMpVRowSegmentSpecial::GetNofPadRows() const
{
/// Return number of pad rows.

  return fPadRows.GetEntriesFast();
}  

//______________________________________________________________________________
AliMpPadRow* AliMpVRowSegmentSpecial::GetPadRow(Int_t i) const
{
/// Return number of pad rows.

  return (AliMpPadRow*)fPadRows[i];
}  

//
// public methods  
//

//______________________________________________________________________________
void  AliMpVRowSegmentSpecial::AddPadRow(AliMpPadRow* padRow)
{
/// Add a pad row.

  padRow->SetOffsetX(fOffsetX);
  padRow->SetID(GetNofPadRows());

  fPadRows.Add(padRow);
}  

//______________________________________________________________________________
void AliMpVRowSegmentSpecial::UpdateMotifVector()
{
/// Add motifs associated with the pad row segments in the specified
/// pad row in the fMotifs vector.

  for (Int_t i=0; i<GetNofPadRows(); i++) {
    AliMpPadRow* padRow = GetPadRow(i);
 
    for (Int_t j=0; j<padRow->GetNofPadRowSegments(); j++) {
      AliMpVMotif* motif = padRow->GetPadRowSegment(j)->GetMotif();            

      if (!HasMotif(motif)) {
        fMotifs.Add(motif);
	
	// resize array if needed
	if (fNofMotifPositionIds<fgkMaxNofMotifPositionIds)
	  fMotifPositionIds.Set(fMotifPositionIds.GetSize()+
	                        fgkMaxNofMotifPositionIds);	 
        fMotifPositionIds.AddAt(
          padRow->GetPadRowSegment(j)->GetMotifPositionId(),
	  fNofMotifPositionIds++);
      }
    }  
  }
}

//______________________________________________________________________________
Double_t  AliMpVRowSegmentSpecial::HalfSizeY() const
{
/// Return the size in y of this row segment.

  Double_t halfSizeY = 0.;
  for (Int_t i=0; i<GetNofPadRows(); i++) {
    halfSizeY += GetPadRow(i)->HalfSizeY();
  }  
  
  return halfSizeY;
}

//______________________________________________________________________________
AliMpVMotif*  AliMpVRowSegmentSpecial::FindMotif(const TVector2& position) const
{
/// Return the motif of this row; 

  AliMpPadRow* padRow 
    = FindPadRow(position.Y());
  
  if (!padRow) return 0;

  AliMpVPadRowSegment* padRowSegment 
    = padRow->FindPadRowSegment(position.X());
    
  if (!padRowSegment) return 0;

  return padRowSegment->GetMotif();
}  

//______________________________________________________________________________
Int_t AliMpVRowSegmentSpecial::FindMotifPositionId(const TVector2& position) const
{
/// Return the motif position identified for the given
/// geometric position.

  AliMpPadRow* padRow 
    = FindPadRow(position.Y());
  
  if (!padRow) return 0;

  AliMpVPadRowSegment* padRowSegment 
    = padRow->FindPadRowSegment(position.X());
    
  if (!padRowSegment) return 0;

  return padRowSegment->GetMotifPositionId();
}

//______________________________________________________________________________
Bool_t AliMpVRowSegmentSpecial::HasMotifPosition(Int_t motifPositionId) const
{
/// Return true if the motif specified with the given position identifier
/// is in this segment.

  if (FindPadRowSegment(motifPositionId))
    return true;
  else  
    return false;   	
}

//______________________________________________________________________________
TVector2 AliMpVRowSegmentSpecial::MotifCenter(Int_t motifPositionId) const
{
/// Return the coordinates of the motif specified with
/// the given position identifier.

  // Try to get the motif position from the motif map first
  AliMpMotifPosition* motifPosition
    = GetRow()->GetMotifMap()->FindMotifPosition(motifPositionId);
  if (motifPosition) return motifPosition->Position();

  // Use slow method otherwise
  return MotifCenterSlow(motifPositionId);
}

//______________________________________________________________________________
TVector2 AliMpVRowSegmentSpecial::Dimensions() const
{
/// Return the halflengths in x, y of the row segment rectangular envelope.

  Double_t x = 0.;		    
  Double_t y = 0.;  
  for (Int_t i=0; i<GetNofPadRows(); i++) {    
    AliMpPadRow* padRow = GetPadRow(i); 
    
    // Add all pad rows y halfsizes   
    y += padRow->HalfSizeY();

    // Find the biggest pad rows x halfsize
    Double_t xx 
      = (padRow->GetPadRowSegment(0)->RightBorderX() -
         padRow->GetPadRowSegment(padRow->GetNofPadRowSegments()-1)->LeftBorderX())/2.;
    if (xx > x) x = xx;		   
  }                  
    
  return TVector2(x, y);   
}

//______________________________________________________________________________
AliMpRow*  AliMpVRowSegmentSpecial::GetRow() const
{
/// Return the row.which this row segment belongs to.

  return fRow;
}  

//______________________________________________________________________________
Int_t  AliMpVRowSegmentSpecial::GetNofMotifs() const 
{ 
/// Return the number of different motifs present in this row segment.

  return fMotifs.GetEntriesFast();
}  

//______________________________________________________________________________
AliMpVMotif* AliMpVRowSegmentSpecial::GetMotif(Int_t i) const  
{
/// Return the i-th motif present in this row segment.

   return (AliMpVMotif*)fMotifs[i]; 
}

//______________________________________________________________________________
Int_t  AliMpVRowSegmentSpecial::GetMotifPositionId(Int_t i) const 
{ 
/// Return the i-th motif position Id present in this row segment.

   return fMotifPositionIds[i]; 
} 

