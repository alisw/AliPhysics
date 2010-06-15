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
// $MpId: AliMpZone.cxx,v 1.7 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpZone
// ---------------
// Class describing a zone composed of the zone segments.
// The zone contains pads of the same dimensions.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpZone.h"
#include "AliMpSubZone.h"

/// \cond CLASSIMP
ClassImp(AliMpZone)
/// \endcond

//_____________________________________________________________________________
AliMpZone::AliMpZone(Int_t id) 
  : TObject(),
    fID(id),
    fSubZones(),
    fPadDimensionX(0.),
    fPadDimensionY(0.)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpZone::AliMpZone() 
  : TObject(),
    fID(0),
    fSubZones(),
    fPadDimensionX(0.),
    fPadDimensionY(0.)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpZone::~AliMpZone() 
{
/// Destructor

  for (Int_t i=0; i<GetNofSubZones(); i++)
    delete fSubZones[i];  
}

//
// public methods
//

//_____________________________________________________________________________
void AliMpZone::AddSubZone(AliMpSubZone* subZone)
{
/// Add row segment.

  fSubZones.Add(subZone);
}  
  
//_____________________________________________________________________________
AliMpSubZone* AliMpZone::FindSubZone(const AliMpVMotif* motif) const
{
/// Find a subzone with a specified motif;
/// return 0 if not found.

  for (Int_t i=0; i<GetNofSubZones(); i++) {
    AliMpSubZone* subZone = GetSubZone(i);
    if (subZone->GetMotif() == motif) return subZone;
  }
  
  return 0;  
}

//_____________________________________________________________________________
void AliMpZone::SetPadDimensions(Double_t dx, Double_t dy)
{ 
/// Set pad dimensions

  fPadDimensionX = dx; 
  fPadDimensionY = dy; 
}

//_____________________________________________________________________________
Int_t AliMpZone::GetNofSubZones() const 
{
/// Return number of row segments.

  return fSubZones.GetEntriesFast();
}  

//_____________________________________________________________________________
AliMpSubZone* AliMpZone::GetSubZone(Int_t i) const 
{
/// Return i-th sub zone.

  if (i<0 || i>=GetNofSubZones()) {
    Warning("GetSubZone", "Index outside range");
    return 0;
  }
  
  return (AliMpSubZone*)fSubZones[i];  
}
