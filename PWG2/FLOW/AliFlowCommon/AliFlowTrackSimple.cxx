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

// AliFlowTrackSimple:
// A simple track class to the the AliFlowEventSimple for flow analysis
//
//
// author: N. van der Kolk (kolk@nikhef.nl)
// mods: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include "TNamed.h"
#include "TParticle.h"
#include "AliFlowTrackSimple.h"
#include "TRandom.h"
#include "AliLog.h"

ClassImp(AliFlowTrackSimple)

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple():
  fEta(0),
  fPt(0),
  fPhi(0),
  fFlowBits(0),
  fSubEventBits(0)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(Double_t phi, Double_t eta, Double_t pt):
  fEta(eta),
  fPt(pt),
  fPhi(phi),
  fFlowBits(0),
  fSubEventBits(0)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(const TParticle* p):
  fEta(p->Eta()),
  fPt(p->Pt()),
  fPhi(p->Phi()),
  fFlowBits(0),
  fSubEventBits(0)
{
  //ctor
}

//-----------------------------------------------------------------------

AliFlowTrackSimple::AliFlowTrackSimple(const AliFlowTrackSimple& aTrack):
  TNamed(),
  fEta(aTrack.fEta),
  fPt(aTrack.fPt),
  fPhi(aTrack.fPhi),
  fFlowBits(aTrack.fFlowBits),
  fSubEventBits(aTrack.fSubEventBits)
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
  fSubEventBits = aTrack.fSubEventBits;

  return *this;
}


//----------------------------------------------------------------------- 
AliFlowTrackSimple::~AliFlowTrackSimple()
{
  //destructor
  
}

//----------------------------------------------------------------------- 
void AliFlowTrackSimple::ResolutionPt(Double_t res)
{
  //smear the pt by a gaussian with sigma=res
  fPt += gRandom->Gaus(0.,res);
}

//----------------------------------------------------------------------- 
void AliFlowTrackSimple::AddV2( Double_t v2, Double_t reactionPlaneAngle, Double_t precisionPhi, Int_t maxNumberOfIterations )
{
  //afterburner, adds v2, uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f,fp,v2sin,v2cos,phiprev;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    v2sin = v2*TMath::Sin(2.*(fPhi-reactionPlaneAngle));
    v2cos = v2*TMath::Cos(2.*(fPhi-reactionPlaneAngle));
    f = fPhi-phi0+v2sin;
    fp = 1.+2.*v2cos; //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}
