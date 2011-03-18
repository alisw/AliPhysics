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
#include "TParticlePDG.h"
#include "AliFlowTrackSimple.h"
#include "TRandom.h"
#include "TMath.h"

ClassImp(AliFlowTrackSimple)

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple():
  TObject(),
  fEta(0),
  fPt(0),
  fPhi(0),
  fTrackWeight(1.),
  fCharge(0),
  fFlowBits(0),
  fSubEventBits(0),
  fID(-1)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(Double_t phi, Double_t eta, Double_t pt, Double_t weight, Int_t charge):
  TObject(),
  fEta(eta),
  fPt(pt),
  fPhi(phi),
  fTrackWeight(weight),
  fCharge(charge),
  fFlowBits(0),
  fSubEventBits(0),
  fID(-1)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple( TParticle* p ):
  TObject(),
  fEta(p->Eta()),
  fPt(p->Pt()),
  fPhi(p->Phi()),
  fTrackWeight(1.),
  fCharge(0),
  fFlowBits(0),
  fSubEventBits(0),
  fID(-1)
{
  //ctor
  TParticlePDG* ppdg = p->GetPDG();
  fCharge = TMath::Nint(ppdg->Charge()/3.0);
}

//-----------------------------------------------------------------------
AliFlowTrackSimple::AliFlowTrackSimple(const AliFlowTrackSimple& aTrack):
  TObject(aTrack),
  fEta(aTrack.fEta),
  fPt(aTrack.fPt),
  fPhi(aTrack.fPhi),
  fTrackWeight(aTrack.fTrackWeight),
  fCharge(aTrack.fCharge),
  fFlowBits(aTrack.fFlowBits),
  fSubEventBits(aTrack.fSubEventBits),
  fID(aTrack.fID)
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
  fCharge = aTrack.fCharge;
  fFlowBits = aTrack.fFlowBits;
  fSubEventBits = aTrack.fSubEventBits;
  fID = aTrack.fID;

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
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

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
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

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
void AliFlowTrackSimple::AddV3( Double_t v3,
                                Double_t reactionPlaneAngle,
                                Double_t precisionPhi,
                                Int_t maxNumberOfIterations )
{
  //afterburner, adds v3, uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0+2./3.*v3*TMath::Sin(3.*(fPhi-reactionPlaneAngle));
    fp = 1.0+2.0*v3*TMath::Cos(3.*(fPhi-reactionPlaneAngle)); //first derivative
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
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0+0.5*v4*TMath::Sin(4.*(fPhi-reactionPlaneAngle));
    fp = 1.0+2.0*v4*TMath::Cos(4.*(fPhi-reactionPlaneAngle)); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//----------------------------------------------------------------------- 
void AliFlowTrackSimple::AddV5( Double_t v5,
                                Double_t reactionPlaneAngle,
                                Double_t precisionPhi,
                                Int_t maxNumberOfIterations )
{
  //afterburner, adds v4, uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0+0.4*v5*TMath::Sin(5.*(fPhi-reactionPlaneAngle));
    fp = 1.0+2.0*v5*TMath::Cos(5.*(fPhi-reactionPlaneAngle)); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//______________________________________________________________________________
void AliFlowTrackSimple::AddFlow( Double_t v1,
                                  Double_t v2,
                                  Double_t v3,
                                  Double_t v4,
                                  Double_t v5,
                                  Double_t rp1,
                                  Double_t rp2,
                                  Double_t rp3,
                                  Double_t rp4,
                                  Double_t rp5,
                                  Double_t precisionPhi,
                                  Int_t maxNumberOfIterations )
{
  //afterburner, adds v1,v2,v4 uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0
        +2.0*  v1*TMath::Sin(    fPhi-rp1)
        +      v2*TMath::Sin(2.*(fPhi-rp2))
        +2./3.*v3*TMath::Sin(3.*(fPhi-rp3))
        +0.5*  v4*TMath::Sin(4.*(fPhi-rp4))
        +0.4*  v5*TMath::Sin(5.*(fPhi-rp5))
        ;
    fp =  1.0
         +2.0*(
           +v1*TMath::Cos(    fPhi-rp1)
           +v2*TMath::Cos(2.*(fPhi-rp2))
           +v3*TMath::Cos(3.*(fPhi-rp3))
           +v4*TMath::Cos(4.*(fPhi-rp4))
           +v5*TMath::Cos(5.*(fPhi-rp5))
         ); //first derivative
    fPhi -= f/fp;
    if (TMath::AreEqualAbs(phiprev,fPhi,precisionPhi)) break;
  }
}

//______________________________________________________________________________
void AliFlowTrackSimple::AddFlow( Double_t v1,
                                  Double_t v2,
                                  Double_t v3,
                                  Double_t v4,
                                  Double_t v5,
                                  Double_t rp,
                                  Double_t precisionPhi,
                                  Int_t maxNumberOfIterations )
{
  //afterburner, adds v1,v2,v4 uses Newton-Raphson iteration
  AddFlow(v1,v2,v3,v4,v5,rp,rp,rp,rp,rp,precisionPhi,maxNumberOfIterations);
}

//______________________________________________________________________________
void AliFlowTrackSimple::AddFlow( Double_t v1,
                                  Double_t v2,
                                  Double_t v3,
                                  Double_t v4,
                                  Double_t reactionPlaneAngle,
                                  Double_t precisionPhi,
                                  Int_t maxNumberOfIterations )
{
  //afterburner, adds v1,v2,v4 uses Newton-Raphson iteration
  Double_t phi0=fPhi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

  for (Int_t i=0; i<maxNumberOfIterations; i++)
  {
    phiprev=fPhi; //store last value for comparison
    f =  fPhi-phi0
        +2.0*  v1*TMath::Sin(    fPhi-reactionPlaneAngle)
        +      v2*TMath::Sin(2.*(fPhi-reactionPlaneAngle))
        +2./3.*v3*TMath::Sin(3.*(fPhi-reactionPlaneAngle))
        +0.5*  v4*TMath::Sin(4.*(fPhi-reactionPlaneAngle))
        ;
    fp =  1.0
         +2.0*(
           +v1*TMath::Cos(    fPhi-reactionPlaneAngle)
           +v2*TMath::Cos(2.*(fPhi-reactionPlaneAngle))
           +v3*TMath::Cos(3.*(fPhi-reactionPlaneAngle))
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
