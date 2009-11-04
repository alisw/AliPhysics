// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTracklet.h"

#include <AliVVertex.h>
#include <TEveTrackPropagator.h>

//______________________________________________________________________________
// AliEveTracklet is a representation of SPD tracklet.
// It inherits from AliEveTrack to allow for common functionality
// regarding track counting.
//
// TEveTrack::fV - stores primary vertex.
// TEveTrack::fP - stores vector in direction of the tracklet with
//                 transverse component equal to 1.

ClassImp(AliEveTracklet)

Float_t AliEveTracklet::fgDefaultRadius = 10;

//______________________________________________________________________________
Float_t AliEveTracklet::GetDefaultRadius()
{
  // Static - return defualt extrapolation radius.

  return fgDefaultRadius;
}

//______________________________________________________________________________
void AliEveTracklet::SetDefaultRadius(Float_t r)
{
  // Static - set defualt extrapolation radius.

  fgDefaultRadius = r;
}

//==============================================================================

//______________________________________________________________________________
AliEveTracklet::AliEveTracklet(Int_t index, const AliVVertex* pv, Float_t theta, Float_t phi,
                               TEveTrackPropagator* prop) :
  AliEveTrack()
{
  // Constructor.

  fIndex = index;
  fV.Set(pv->GetX(), pv->GetY(), pv->GetZ());
  fP.Set(TMath::Cos(phi), TMath::Sin(phi), 1.0/TMath::Tan(theta));

  if (prop) SetPropagator(prop);
}

//==============================================================================

//______________________________________________________________________________
void AliEveTracklet::MakeTrack(Bool_t recurse)
{
  // Make track -- just make a line to radius specified in propagator
  // or use the default if it is not set.

  Float_t r = fPropagator ? fPropagator->GetMaxR() : fgDefaultRadius;
  Reset(2);
  SetPoint(0, fV.fX, fV.fY, fV.fZ);
  SetPoint(1, fV.fX + r*fP.fX, fV.fY + r*fP.fY, fV.fZ + r*fP.fZ);

  if (recurse)
  {
    for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i)
    {
      TEveTrack* t = dynamic_cast<TEveTrack*>(*i);
      if (t) t->MakeTrack(recurse);
    }
  }
}

//______________________________________________________________________________
void AliEveTracklet::SecSelected(TEveTrack* track)
{
  // Emits "SecSelected(TEveTrack*)" signal.
  // Called from TEveTrackGL on secondary-selection.

  Emit("SecSelected(TEveTrack*)", (Long_t)track);
  SecSelectedTracklet((AliEveTracklet*) track);
}

//______________________________________________________________________________
void AliEveTracklet::SecSelectedTracklet(AliEveTracklet* track)
{
  // Emits "SecSelectedTracklet(AliEveTracklet*)" signal.

  Emit("SecSelectedTracklet(AliEveTracklet*)", (Long_t)track);
}
