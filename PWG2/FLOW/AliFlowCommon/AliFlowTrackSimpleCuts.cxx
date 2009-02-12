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
#include "TObject.h"
#include "AliFlowTrackSimpleCuts.h"


// AliFlowTrackSimpleCuts:
// A simple track cut class to the the AliFlowTrackSimple 
// for basic kinematic cuts
//
// author: N. van der Kolk (kolk@nikhef.nl)


ClassImp(AliFlowTrackSimpleCuts)

//-----------------------------------------------------------------------

AliFlowTrackSimpleCuts::AliFlowTrackSimpleCuts():
  fPtMax(0.),
  fPtMin(0.),
  fEtaMax(0.),
  fEtaMin(0.),
  fPhiMax(0.),
  fPhiMin(0.),
  fPID(0)
  
{
  //constructor 
  
}

//-----------------------------------------------------------------------

AliFlowTrackSimpleCuts::AliFlowTrackSimpleCuts(const AliFlowTrackSimpleCuts& someCuts):
  TObject(),
  fPtMax(someCuts.fPtMax),
  fPtMin(someCuts.fPtMin),
  fEtaMax(someCuts.fEtaMax),
  fEtaMin(someCuts.fEtaMin),
  fPhiMax(someCuts.fPhiMax),
  fPhiMin(someCuts.fPhiMin),
  fPID(someCuts.fPID)
{
  //copy constructor 
}

//-----------------------------------------------------------------------

AliFlowTrackSimpleCuts& AliFlowTrackSimpleCuts::operator=(const AliFlowTrackSimpleCuts& someCuts)
{
  fPtMax  = someCuts.fPtMax;
  fPtMin  = someCuts.fPtMin;
  fEtaMax = someCuts.fEtaMax;
  fEtaMin = someCuts.fEtaMin;
  fPhiMax = someCuts.fPhiMax;
  fPhiMin = someCuts.fPhiMin;
  fPID    = someCuts.fPID;

  return *this;

}


//----------------------------------------------------------------------- 

AliFlowTrackSimpleCuts::~AliFlowTrackSimpleCuts()
{
  //destructor
  
}
