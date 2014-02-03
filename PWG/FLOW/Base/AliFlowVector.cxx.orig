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

//********************************************************************
// AliFlowVector:                                                    *
// Class to hold the flow vector and multiplicity for flow analysis. *
// Author: A. Bilandzic (anteb@nikhef.nl)                            *
//********************************************************************

ClassImp(AliFlowVector)

//________________________________________________________________________

AliFlowVector::AliFlowVector():
 TVector2(0.0,0.0),fMult(0.0)
{
  // default constructor
}

//________________________________________________________________________

AliFlowVector::AliFlowVector(const AliFlowVector& aVector):
  TVector2(aVector),
  fMult(aVector.fMult)
{
  // copy constructor
}

//________________________________________________________________________

AliFlowVector::AliFlowVector(Double_t *y, Double_t m):
  TVector2(y),
  fMult(m)
{
  // Analogue of TVector2 constructor. Sets (x,y) and multiplicity 1.
}

 //________________________________________________________________________

AliFlowVector::AliFlowVector(const TVector2 &v, Double_t m):
  TVector2(v),
  fMult(m)
{
  // custom constructor, Sets vector and multiplicity
}

 //________________________________________________________________________

AliFlowVector::AliFlowVector(Double_t x, Double_t y, Double_t m):
  TVector2(x,y),
  fMult(m)
{
  // custom constructor analogue of TVector2 constructor
}

//________________________________________________________________________ 

AliFlowVector::~AliFlowVector()
{
  // default destructor 
}

void AliFlowVector::SetMagPhi(double size, double angle, double mult)
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

AliFlowVector& AliFlowVector::operator*=(double w)
{
   // multiply by a weight operator
   fX*=w;
   fY*=w;
   fMult*=w;
   return *this;
}
