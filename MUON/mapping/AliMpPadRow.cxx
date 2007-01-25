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
// $MpId: AliMpPadRow.cxx,v 1.8 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpPadRow
// ------------------
// Class describing a pad row composed of the pad row segments.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpPadRow.h"
#include "AliMpPadRowLSegment.h"
#include "AliMpPadRowRSegment.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpPadRow)
/// \endcond

//_____________________________________________________________________________
AliMpPadRow::AliMpPadRow(AliMp::XDirection direction) 
  : TObject(),
    fDirection(direction), 
    fID(0),
    fOffsetX(0),
    fSegments() 
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpPadRow::AliMpPadRow() 
  : TObject(),
    fDirection(AliMp::kLeft), 
    fID(0),
    fOffsetX(0),
    fSegments() 
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpPadRow::~AliMpPadRow() 
{
/// Destructor  

  for (Int_t i=0; i<GetNofPadRowSegments() ; i++)
    delete fSegments[i];
}

//
// private methods
//

//_____________________________________________________________________________
Double_t AliMpPadRow::CurrentBorderX() const
{
/// Return the left/right x border 
/// (depending on the direction which the row segments are filled in).

  if (GetNofPadRowSegments() == 0)
      return fOffsetX;
  else 
    if (fDirection == AliMp::kLeft)
      return GetPadRowSegment(GetNofPadRowSegments()-1)->LeftBorderX();
    else  
      return GetPadRowSegment(GetNofPadRowSegments()-1)->RightBorderX();
}

//
// public methods
//

//_____________________________________________________________________________
AliMpVPadRowSegment* 
AliMpPadRow::AddPadRowSegment(AliMpMotif* motif, Int_t motifPositionId,
                              Int_t nofPads)
{
/// Add a pad row segment.

  AliMpVPadRowSegment* padRowSegment = 0;

  if (fDirection == AliMp::kLeft) {
    padRowSegment 
      = new AliMpPadRowLSegment(this, motif, motifPositionId, nofPads);
  }    
  else  {
    padRowSegment 
      = new AliMpPadRowRSegment(this, motif, motifPositionId, nofPads);
  }     

  // Set pad row segment offset
  padRowSegment->SetOffsetX(CurrentBorderX());

  // Adds the pad row segment
#ifdef WITH_STL
  fSegments.push_back(padRowSegment);
#endif

#ifdef WITH_ROOT
  fSegments.Add(padRowSegment);
#endif
  
  return padRowSegment;
}  
  
//_____________________________________________________________________________
AliMpVPadRowSegment* AliMpPadRow::FindPadRowSegment(Double_t x) const
{
/// Find the row segment for the specified x position;
/// return 0 if no row segment is found.

  for (Int_t i=0; i<GetNofPadRowSegments(); i++) {
    AliMpVPadRowSegment* rs = GetPadRowSegment(i);
    if (x >= rs->LeftBorderX() && x <= rs->RightBorderX())
      return rs;
  }
  
  return 0;    
}    

//_____________________________________________________________________________
Double_t  AliMpPadRow::HalfSizeY() const
{
/// Return the half size in y

  return GetPadRowSegment(0)->HalfSizeY();
}

//_____________________________________________________________________________
void  AliMpPadRow::SetID(Int_t id)
{
/// Set the ID.

  fID = id;
}    

//_____________________________________________________________________________
void  AliMpPadRow::SetOffsetX(Double_t offsetX)
{
/// Set the x offset.

  fOffsetX = offsetX;
}    

//_____________________________________________________________________________
Int_t AliMpPadRow::GetID() const 
{
/// Return the pad row ID.

  return fID;
}  

//_____________________________________________________________________________
Int_t AliMpPadRow::GetNofPadRowSegments() const 
{
/// Return the number of pad row segments.

#ifdef WITH_STL
  return fSegments.size();
#endif

#ifdef WITH_ROOT
  return fSegments.GetEntriesFast();
#endif
}  

//_____________________________________________________________________________
AliMpVPadRowSegment* AliMpPadRow::GetPadRowSegment(Int_t i) const 
{
/// Return the pad row segment with the specified number.

  if (i<0 || i>=GetNofPadRowSegments()) {
    AliWarningStream() << "Index outside range" << endl;
    return 0;
  }
  
#ifdef WITH_STL
  return fSegments[i];  
#endif

#ifdef WITH_ROOT
  return (AliMpVPadRowSegment*)fSegments[i];  
#endif
}

//_____________________________________________________________________________
Int_t AliMpPadRow::GetNofPads() const 
{
/// Return the number of pads in this pad row.

  Int_t nofPads=0;
  for (Int_t i=0; i<GetNofPadRowSegments(); i++)
    nofPads += GetPadRowSegment(i)->GetNofPads();

  return nofPads;
}  

