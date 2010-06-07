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
// A simple track class to the the AliFlowEventSimple for flow analysis
//
//
// author: N. van der Kolk (kolk@nikhef.nl)
// mods: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include "AliFlowTrack.h"

ClassImp(AliFlowTrack)

//-----------------------------------------------------------------------
AliFlowTrack::AliFlowTrack():
  AliFlowTrackSimple(0),
  fTrackSourceBits(),
  fFMDmultiplicity(0.)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrack::AliFlowTrack(const AliFlowTrack& aTrack):
  AliFlowTrackSimple(aTrack),
  fTrackSourceBits(aTrack.fTrackSourceBits),
  fFMDmultiplicity(aTrack.fFMDmultiplicity)
{
  //copy constructor 
}

//-----------------------------------------------------------------------
AliFlowTrack& AliFlowTrack::operator=(const AliFlowTrack& aTrack)
{
  AliFlowTrackSimple::operator=(aTrack);
  fTrackSourceBits = aTrack.fTrackSourceBits;
  fFMDmultiplicity = aTrack.fFMDmultiplicity;

  return *this;
}


//----------------------------------------------------------------------- 
AliFlowTrack::~AliFlowTrack()
{
  //destructor
}

