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

/* $Id$ */ 

// AliFlowTrack:
// A track class for use in AliFlowEvent for flow analysis
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include "AliVParticle.h"
#include "AliFlowTrack.h"
#include "AliFlowTrackSimple.h"

ClassImp(AliFlowTrack)

//-----------------------------------------------------------------------
AliFlowTrack::AliFlowTrack():
  AliFlowTrackSimple(),
  fTrackSourceBits() 
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrack::AliFlowTrack(const AliVParticle* p):
  AliFlowTrackSimple(),
  fTrackSourceBits()
{
  //constructor 
  Set(p);
}

//-----------------------------------------------------------------------
void AliFlowTrack::Set(const AliVParticle* p)
{
  //set from an AliVParticle
  SetPhi(p->Phi());
  SetEta(p->Eta());
  SetPt(p->Pt());
  SetCharge(p->Charge());
}

//-----------------------------------------------------------------------
AliFlowTrack::AliFlowTrack(const AliFlowTrack& aTrack):
  AliFlowTrackSimple(aTrack),
  fTrackSourceBits(aTrack.fTrackSourceBits)
{
  //copy constructor 
}

//-----------------------------------------------------------------------
AliFlowTrack* AliFlowTrack::Clone(const char* /*option*/) const
{
  //clone "constructor"
  return new AliFlowTrack(*this);
}

//-----------------------------------------------------------------------
AliFlowTrack& AliFlowTrack::operator=(const AliFlowTrack& aTrack)
{
  //assignment
  AliFlowTrackSimple::operator=(aTrack);
  fTrackSourceBits = aTrack.fTrackSourceBits;
  return *this;
}

////-----------------------------------------------------------------------
//AliFlowTrackSimple& AliFlowTrack::operator=(const AliFlowTrackSimple& aTrack)
//{
//  //polymorphic assignment
//  AliFlowTrackSimple::operator=(aTrack);
//  const AliFlowTrack* pft = dynamic_cast<const AliFlowTrack*>(&aTrack);
//  if (pft)
//  {
//    fTrackSourceBits = pft->fTrackSourceBits;
//  }
//  else
//  {
//    fTrackSourceBits.ResetAllBits();
//  }
//  return *this;
//}

//----------------------------------------------------------------------- 
AliFlowTrack::~AliFlowTrack()
{
  //destructor
}

