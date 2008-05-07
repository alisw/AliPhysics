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

#include "AliFlowTrackSimple.h"


// AliFlowTrackSimple:
// A simple track class to the the AliFlowEventSimple for flow analysis
//
//
// author: N. van der Kolk (kolk@nikhef.nl)


ClassImp(AliFlowTrackSimple)

//-----------------------------------------------------------------------

AliFlowTrackSimple::AliFlowTrackSimple(): 
  fEta(0),
  fPt(0),
  fPhi(0),
  fFlowBits(0)
{
  //constructor 
  
}

//-----------------------------------------------------------------------

AliFlowTrackSimple::AliFlowTrackSimple(const AliFlowTrackSimple& aTrack):
  TObject(),
  fEta(aTrack.fEta),
  fPt(aTrack.fPt),
  fPhi(aTrack.fPhi),
  fFlowBits(aTrack.fFlowBits)
{
  //copy constructor 
}
//-----------------------------------------------------------------------

AliFlowTrackSimple& AliFlowTrackSimple::operator=(const AliFlowTrackSimple& aTrack)
{
  fEta = aTrack.fEta;
  fPt = aTrack.fPt;
  fPhi = aTrack.fPhi;
  fFlowBits = aTrack.fFlowBits;

  return *this;

}


//----------------------------------------------------------------------- 

AliFlowTrackSimple::~AliFlowTrackSimple()
{
  //destructor
  
}

