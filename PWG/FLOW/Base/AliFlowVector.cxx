/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "TMath.h"

//********************************************************************
// AliFlowVector:                                                    *
// Class to hold the flow vector and multiplicity for flow analysis. *
// Author: A. Bilandzic (anteb@nikhef.nl)                            *
// extended: M.Krzewicki (mikolaj.krzewicki@cern.ch)                 *
//********************************************************************

ClassImp(AliFlowVector)

//________________________________________________________________________

AliFlowVector::AliFlowVector():
  TVector2(0.0,0.0),
  fMult(0.0),
  fHarmonic(2),
  fPOItype(0),
  fSubeventNumber(-1)
{
  // default constructor
}

//________________________________________________________________________

AliFlowVector::AliFlowVector(const AliFlowVector& aVector):
  TVector2(aVector),
  fMult(aVector.fMult),
  fHarmonic(aVector.fHarmonic),
  fPOItype(aVector.fPOItype),
  fSubeventNumber(aVector.fSubeventNumber)
{
  // copy constructor
}

//________________________________________________________________________

AliFlowVector::AliFlowVector(Double_t *y, Double_t m, Int_t h, Int_t t, Int_t s):
  TVector2(y),
  fMult(m),
  fHarmonic(h),
  fPOItype(t),
  fSubeventNumber(s)
{
  // Analogue of TVector2 constructor. Sets (x,y) and multiplicity 1.
}

 //________________________________________________________________________

AliFlowVector::AliFlowVector(const TVector2 &v, Double_t m, Int_t h, Int_t t, Int_t s):
  TVector2(v),
  fMult(m),
  fHarmonic(h),
  fPOItype(t),
  fSubeventNumber(s)
{
  // custom constructor, Sets vector and multiplicity
}

 //________________________________________________________________________

AliFlowVector::AliFlowVector(Double_t x, Double_t y, Double_t m, Int_t h, Int_t t, Int_t s):
  TVector2(x,y),
  fMult(m),
  fHarmonic(h),
  fPOItype(t),
  fSubeventNumber(s)
{
  // custom constructor analogue of TVector2 constructor
}

//________________________________________________________________________ 

AliFlowVector::~AliFlowVector()
{
  // default destructor 
}

void AliFlowVector::SetMagPhi(Double_t size, Double_t angle, Double_t mult)
{
   // Analogue to SetMagPhi for a TVector2 but here including a sum of weights
   TVector2::SetMagPhi(size,angle);
   SetMult(mult);
}

//________________________________________________________________________

AliFlowVector& AliFlowVector::operator=(const AliFlowVector& aVector)
{
  // assignement operator
  if (this==&aVector) return *this;
  fX = aVector.X();
  fY = aVector.Y();
  fMult = aVector.GetMult();
  return *this;
}

//________________________________________________________________________

AliFlowVector& AliFlowVector::operator+=(const AliFlowVector& aVector)
{
  // addition operator
  fX += aVector.X(); 
  fY += aVector.Y(); 
  fMult += aVector.GetMult(); 
  return *this;
}

AliFlowVector& AliFlowVector::operator-=(const AliFlowVector& aVector)
{
  // subtraction operator
  fX -= aVector.X(); 
  fY -= aVector.Y(); 
  fMult -= aVector.GetMult();
  return *this;
}

AliFlowVector& AliFlowVector::operator*=(Double_t w)
{
   // multiply by a weight operator
   fX*=w;
   fY*=w;
   fMult*=w;
   return *this;
}

//________________________________________________________________________
void AliFlowVector::Clear(Option_t* /*option*/)
{
  //clear
  fX=0.;
  fY=0.;
  fMult=0;
  fHarmonic=2;
  fPOItype=AliFlowTrackSimple::kRP;
  fSubeventNumber=-1;
}

//________________________________________________________________________
Int_t AliFlowVector::SubtractTrackWithDaughters( const AliFlowTrackSimple* track, 
                                                Double_t extraWeight 
                                              )
{
  //subtract a track and all its daughters, only if tagged with flowTag and in specified
  //subevent (-1 for no subevent selection)
  //to only subtract if it was actually used in the construction of the vector)
  //TODO: maybe make recursive if it ever becomes needed
  //for complicated decay topologies
  Bool_t inSubEvent=kTRUE;
  if (fSubeventNumber>=0) 
  {
    inSubEvent = track->InSubevent(fSubeventNumber);
  }
  if (track->IsPOItype(fPOItype) && inSubEvent )
  {
    fX -= extraWeight * track->Weight() * TMath::Cos(fHarmonic*track->Phi());
    fY -= extraWeight * track->Weight() * TMath::Sin(fHarmonic*track->Phi());
  }
  
  Int_t numberOfsubtractedDaughters=0;
  for (Int_t i=0; i<track->GetNDaughters(); i++)
  {
    AliFlowTrackSimple* daughter = track->GetDaughter(i);
    if (!daughter) continue;
    inSubEvent=kTRUE;
    if (fSubeventNumber>=0) 
    {
      inSubEvent = daughter->InSubevent(fSubeventNumber);
    }
    if (daughter->IsPOItype(fPOItype) && inSubEvent )
    {
      fX -= extraWeight * daughter->Weight() * TMath::Cos(fHarmonic*daughter->Phi());
      fY -= extraWeight * daughter->Weight() * TMath::Sin(fHarmonic*daughter->Phi());
      numberOfsubtractedDaughters++;
    }
  }
  return numberOfsubtractedDaughters;
}
