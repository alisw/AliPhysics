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

/* $Id: AliTRDtrackletGTU.cxx 28397 2008-09-02 09:33:00Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  GTU tracklet                                                          //
//                                                                        //
//  Author: J. Klein (Jochen.Klein@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TClass.h"

#include "AliTRDtrackletGTU.h"
#include "AliTRDtrackletWord.h"
#include "AliTRDmcmTracklet.h"
#include "AliTRDtrackletMCM.h"
#include "AliLog.h"
#include "AliTRDgtuParam.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"

ClassImp(AliTRDtrackletGTU)

AliTRDtrackletBase* AliTRDtrackletGTU::fgkDummyTracklet = new AliTRDtrackletWord(0);

AliTRDtrackletGTU::AliTRDtrackletGTU() :
  AliTRDtrackletBase(),
  fGtuParam(AliTRDgtuParam::Instance()),
  fTracklet(0x0), //fgkDummyTracklet), 
  fSubChannel(0x0),
  fAssignedZ(kFALSE),
  fAlpha(0),
  fYProj(0),
  fYPrime(0),
  fIndex(0)
{
  // ctor for any tracklet deriving from AliTRDtrackletBase

  fSubChannel = new Int_t[fGtuParam->GetNZChannels()];
  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) 
    fSubChannel[zch] = 0;
}

AliTRDtrackletGTU::AliTRDtrackletGTU(AliTRDtrackletBase *tracklet) :
  AliTRDtrackletBase(*tracklet),
  fGtuParam(AliTRDgtuParam::Instance()),
  fTracklet(0x0), 
  fSubChannel(0x0),
  fAssignedZ(kFALSE),
  fAlpha(0),
  fYProj(0),
  fYPrime(0),
  fIndex(0)
{
  // ctor for any tracklet deriving from AliTRDtrackletBase

  fSubChannel = new Int_t[fGtuParam->GetNZChannels()];
  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) 
    fSubChannel[zch] = 0;
  fTracklet = tracklet;
  if ( fTracklet->IsA() == TClass::GetClass("AliTRDtrackletMCM")) {
      AliDebug(5,Form("label from mcm tracklet: %i", ((AliTRDtrackletMCM*) fTracklet)->GetLabel()));
  }
}

AliTRDtrackletGTU::AliTRDtrackletGTU(const AliTRDtrackletGTU& tracklet) :
  AliTRDtrackletBase(tracklet),
  fGtuParam(AliTRDgtuParam::Instance()),
  fTracklet(tracklet.fTracklet), 
  fSubChannel(0x0),
  fAssignedZ(tracklet.fAssignedZ),
  fAlpha(tracklet.fAlpha),
  fYProj(tracklet.fYProj),
  fYPrime(tracklet.fYPrime), 
  fIndex(tracklet.fIndex)
{
  // copy ctor

  fSubChannel = new Int_t[fGtuParam->GetNZChannels()];
  for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) 
    fSubChannel[zch] = tracklet.fSubChannel[zch];
}

AliTRDtrackletGTU& AliTRDtrackletGTU::operator=(const AliTRDtrackletGTU &rhs)
{
  if (&rhs != this) {
    fTracklet = rhs.fTracklet;
    for (Int_t zch = 0; zch < fGtuParam->GetNZChannels(); zch++) 
      fSubChannel[zch] = rhs.fSubChannel[zch];
    fIndex = rhs.fIndex;
    fYPrime = rhs.fYPrime;
    fYProj = rhs.fYProj;
    fAlpha = rhs.fAlpha;
    fAssignedZ = rhs.fAssignedZ;
  }
  
  return *this;
}

AliTRDtrackletGTU::~AliTRDtrackletGTU() 
{
  // dtor
  if (fSubChannel)
    delete [] fSubChannel; 
  fTracklet = 0x0;
}

Int_t AliTRDtrackletGTU::Compare(const TObject *o) const {
  // sorting w. r. t. Z, Y if no z-channel assigned (as needed in input unit)
  // otherwise w. r. t. Z, Y_proj
  // must be changed to Z-channel, Y_proj
  // will be changed

  if (!o) 
    return 0;

  if (!o->InheritsFrom("AliTRDtrackletGTU")) {
      AliError("Cannot compare to object not deriving from AliTRDtrackletGTU");
      return 0;
  }

  if (!fAssignedZ) {
    if ( GetZbin() < ((AliTRDtrackletGTU*) o)->GetZbin()) 
      return -1;
    else if (GetZbin() > ((AliTRDtrackletGTU*) o)->GetZbin()) 
      return 1;
    else 
      if (GetYbin() < ((AliTRDtrackletGTU*) o)->GetYbin())
	return -1;
      else if (GetYbin() > ((AliTRDtrackletGTU*) o)->GetYbin())
	return 1;
      else 
	return 0;
  }
  else {
    // sorting should be according to zsubindex, not to Z !!!
    // therefore this depends on the zch
    if (GetZbin() < ((AliTRDtrackletGTU*) o)->GetZbin()) 
      return -1;
    else if (GetZbin() > ((AliTRDtrackletGTU*) o)->GetZbin()) 
      return 1;
    else 
      if (GetYProj() < ((AliTRDtrackletGTU*) o)->GetYProj())
	return -1;
      else if (GetYProj() > ((AliTRDtrackletGTU*) o)->GetYProj())
	return 1;
      else 
	return 0;
  }
}

void AliTRDtrackletGTU::SetSubChannel(Int_t zch, Int_t subch) 
{
  // set the subchannel in the given z-channel
  fAssignedZ = kTRUE;
  fSubChannel[zch] = subch;
}

Int_t AliTRDtrackletGTU::GetSubChannel(Int_t zch) 
{
  // get the subchannel in the given z-channel
  return fSubChannel[zch];
}

Int_t AliTRDtrackletGTU::GetLabel() const
{
    if ( fTracklet->IsA() == TClass::GetClass("AliTRDtrackletMCM"))
	return ((AliTRDtrackletMCM*) fTracklet)->GetLabel();
    else
	return -1;
}

/*
Float_t AliTRDtrackletGTU::GetPhysX(Int_t layer) 
{
  // get the x-position (in the local system) assuming the tracklet is in the given layer
  return fGtuParam->GetGeo()->GetTime0(layer);
}

Float_t AliTRDtrackletGTU::GetPhysY() 
{
  // 
  return GetYbin() * 0.0160; 
}

Float_t AliTRDtrackletGTU::GetPhysAlpha() 
{
  return GetAlpha() * 0.01; // wrong factor!
}

Float_t AliTRDtrackletGTU::GetPhysZ(Int_t stack, Int_t layer) 
{
  return fGtuParam->GetGeo()->GetPadPlane(layer, stack)->GetRowPos(GetZbin()); // not the middle of a pad!
}
*/
