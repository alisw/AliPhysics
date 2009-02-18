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

#define AliFlowVector_cxx

#include "AliFlowVector.h"

//********************************************************************
// AliFlowVector:                                                    *
// Class to hold the flow vector and multiplicity for flow analysis. *
// Author: A. Bilandzic (anteb@nikhef.nl)                            *
//********************************************************************

ClassImp(AliFlowVector)

//________________________________________________________________________

AliFlowVector::AliFlowVector():
 fMult(0),
 fSumOfWeightsToPower2(0),
 fSumOfWeightsToPower3(0),
 fSumOfWeightsToPower4(0),
 fSumOfWeightsToPower5(0),
 fSumOfWeightsToPower6(0),
 fSumOfWeightsToPower7(0),
 fSumOfWeightsToPower8(0) 
 {
  // default contructor
 }

AliFlowVector::AliFlowVector(const AliFlowVector& aVector):
 TVector2(aVector),
 fMult(aVector.fMult),
 fSumOfWeightsToPower2(aVector.fSumOfWeightsToPower2),
 fSumOfWeightsToPower3(aVector.fSumOfWeightsToPower3),
 fSumOfWeightsToPower4(aVector.fSumOfWeightsToPower4),
 fSumOfWeightsToPower5(aVector.fSumOfWeightsToPower5),
 fSumOfWeightsToPower6(aVector.fSumOfWeightsToPower6),
 fSumOfWeightsToPower7(aVector.fSumOfWeightsToPower7),
 fSumOfWeightsToPower8(aVector.fSumOfWeightsToPower8)
 {
  // copy constructor
 }
 
AliFlowVector::AliFlowVector(const TVector2 &v, const Double_t m, const Double_t sumPow2w, const Double_t sumPow3w, const Double_t sumPow4w, const Double_t sumPow5w, const Double_t sumPow6w, const Double_t sumPow7w, const Double_t sumPow8w):
 TVector2(v),
 fMult(m),
 fSumOfWeightsToPower2(sumPow2w),
 fSumOfWeightsToPower3(sumPow3w),
 fSumOfWeightsToPower4(sumPow4w),
 fSumOfWeightsToPower5(sumPow5w),
 fSumOfWeightsToPower6(sumPow6w),
 fSumOfWeightsToPower7(sumPow7w),
 fSumOfWeightsToPower8(sumPow8w) 
 {
  // custom constructor
 }
 
AliFlowVector::~AliFlowVector()
{
 // default constructor 
}
AliFlowVector& AliFlowVector::operator=(const AliFlowVector& aVector)
{
 // assignement operator
 fX = aVector.X();
 fY = aVector.Y();
 fMult = aVector.GetMult();
 fSumOfWeightsToPower2 = aVector.GetSumOfWeightsToPower2();
 fSumOfWeightsToPower3 = aVector.GetSumOfWeightsToPower3();
 fSumOfWeightsToPower4 = aVector.GetSumOfWeightsToPower4();
 fSumOfWeightsToPower5 = aVector.GetSumOfWeightsToPower5();
 fSumOfWeightsToPower6 = aVector.GetSumOfWeightsToPower6();
 fSumOfWeightsToPower7 = aVector.GetSumOfWeightsToPower7();
 fSumOfWeightsToPower8 = aVector.GetSumOfWeightsToPower8();
 
 return *this;
}


