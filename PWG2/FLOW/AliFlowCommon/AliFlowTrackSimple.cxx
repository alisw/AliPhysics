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

#include "TObject.h"
#include "TParticle.h"
#include "AliFlowTrackSimple.h"
#include "TRandom.h"
#include "AliLog.h"

ClassImp(AliFlowTrackSimple)

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple():
  TObject(),
  fEta(0),
  fPt(0),
  fPhi(0),
  fTrackWeight(1.),
  fFlowBits(0),
  fSubEventBits(0)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(Double_t phi, Double_t eta, Double_t pt, Double_t weight):
  TObject(),
  fEta(eta),
  fPt(pt),
  fPhi(phi),
  fTrackWeight(weight),
  fFlowBits(0),
  fSubEventBits(0)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(const TParticle* p):
  TObject(),
  fEta(p->Eta()),
  fPt(p->Pt()),
  fPhi(p->Phi()),
  fTrackWeight(1.),
  fFlowBits(0),
  fSubEventBits(0)
{
  //ctor
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(const AliFlowTrackSimple& aTrack):
  TObject(aTrack),
  fEta(aTrack.fEta),
  fPt(aTrack.fPt),
  fPhi(aTrack.fPhi),
  fTrackWeight(aTrack.fTrackWeight),
  fFlowBits(aTrack.fFlowBits),
  fSubEventBits(aTrack.fSubEventBits)
{
  //copy constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple* AliFlowTrackSimple::Clone(const char* /*option*/) const
{
  //clone "constructor"
  return new AliFlowTrackSimple(*this);
}

//-----------------------------------------------------------------------
AliFlowTrackSimple& AliFlowTrackSimple::operator=(const AliFlowTrackSimple& aTrack)
{
  fEta = aTrack.fEta;
  fPt = aTrack.fPt;
  fPhi = aTrack.fPhi;
  fTrackWeight = aTrack.fTrackWeight;
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
void AliFlowTrackSimple::AddV1( Double_t v1,
                                Double_t reactionPlaneAngle,
                                Double_t precisionPhi,
                                Int_t maxNumberOfIterations )
{
  //afterburner, adds v1, uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f,fp,phiprev;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0+2.0*v1*TMath::Sin(fPhi-reactionPlaneAngle);
    fp = 1.0+2.0*v1*TMath::Cos(fPhi-reactionPlaneAngle); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//----------------------------------------------------------------------- 
void AliFlowTrackSimple::AddV2( Double_t v2,
                                Double_t reactionPlaneAngle,
                                Double_t precisionPhi,
                                Int_t maxNumberOfIterations )
{
  //afterburner, adds v2, uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f,fp,phiprev;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0+v2*TMath::Sin(2.*(fPhi-reactionPlaneAngle));
    fp = 1.0+2.0*v2*TMath::Cos(2.*(fPhi-reactionPlaneAngle)); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//----------------------------------------------------------------------- 
void AliFlowTrackSimple::AddV4( Double_t v4,
                                Double_t reactionPlaneAngle,
                                Double_t precisionPhi,
                                Int_t maxNumberOfIterations )
{
  //afterburner, adds v4, uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f,fp,phiprev;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0+0.5*v4*TMath::Sin(4.*(fPhi-reactionPlaneAngle));
    fp = 1.0+2.0*v4*TMath::Cos(4.*(fPhi-reactionPlaneAngle)); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//______________________________________________________________________________
void AliFlowTrackSimple::AddFlow( Double_t v1,
                                  Double_t v2,
                                  Double_t v4,
                                  Double_t reactionPlaneAngle,
                                  Double_t precisionPhi,
                                  Int_t maxNumberOfIterations )
{
  //afterburner, adds v1,v2,v4 uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f,fp,phiprev;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0
        +2.0*v1*TMath::Sin(fPhi-reactionPlaneAngle)
        +    v2*TMath::Sin(2.*(fPhi-reactionPlaneAngle))
        +0.5*v4*TMath::Sin(4.*(fPhi-reactionPlaneAngle))
        ;
    fp =  1.0
         +2.0*(
           +v1*TMath::Cos(fPhi-reactionPlaneAngle)
           +v2*TMath::Cos(2.*(fPhi-reactionPlaneAngle))
           +v4*TMath::Cos(4.*(fPhi-reactionPlaneAngle))
         ); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//______________________________________________________________________________
void AliFlowTrackSimple::Print( Option_t* /*option*/ ) const
{
  //print stuff
  printf("Phi: %.3f, Eta: %+.3f, Pt: %.3f, weight: %.3f",fPhi,fEta,fPt,fTrackWeight);
  if (InRPSelection()) printf(", RP");
  if (InPOISelection()) printf(", POI");
  for (Int_t i=0; i<2; i++)
  {
    if (InSubevent(i)) printf(", subevent %i",i);
  }
  printf("\n");
}
