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

AliFlowVector::AliFlowVector():fMult(0){}

AliFlowVector::AliFlowVector(const AliFlowVector& aVector):TVector2(aVector),
fMult(aVector.fMult){}
AliFlowVector::AliFlowVector(const TVector2 &v, Double_t m):TVector2(v),fMult(m){}
AliFlowVector::~AliFlowVector(){}
AliFlowVector& AliFlowVector::operator=(const AliFlowVector& aVector)
{
   //assignement operator
   fX=aVector.X();
   fY=aVector.Y();
   fMult=aVector.GetMult();
   return *this;
}


