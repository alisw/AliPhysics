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

//-------------------------------------------------------------------------
//                Implementation of the HLT ITS track class
//
//          Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliL3ITStrack.h"

ClassImp(AliL3ITStrack)

//____________________________________________________________________________
AliL3ITStrack::AliL3ITStrack()
              :AliITStrackV2()
{
  //------------------------------------------------------------------
  //Constructor
  //------------------------------------------------------------------
}

//____________________________________________________________________________
AliL3ITStrack::AliL3ITStrack(AliESDtrack& t)
              :AliITStrackV2(t)
{
  //------------------------------------------------------------------
  //Constructor
  //------------------------------------------------------------------
}

//____________________________________________________________________________
AliL3ITStrack::AliL3ITStrack(const AliL3ITStrack& t) 
              : AliITStrackV2(t)
{
  //------------------------------------------------------------------
  //Copy constructor
  //------------------------------------------------------------------
}

//_____________________________________________________________________________
Int_t AliL3ITStrack::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliL3ITStrack *t=(AliL3ITStrack*)o;
  Double_t co=TMath::Abs(t->Get1Pt());
  Double_t c =TMath::Abs(Get1Pt());
  //  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  //  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

Bool_t AliL3ITStrack::GetPxPyPzAt(Double_t x,Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // at the position "x" using the helix track approximation
  //---------------------------------------------------------------------
  p[0]=fP4;
  p[1]=fP2+(x-fX)*fP4/AliKalmanTrack::GetConvConst(); 
  p[2]=fP3;

  if (TMath::Abs(p[0])<=0)        return kFALSE;
  if (TMath::Abs(p[1])> 0.999999) return kFALSE;

  Double_t pt=1./TMath::Abs(p[0]);
  Double_t cs=TMath::Cos(fAlpha), sn=TMath::Sin(fAlpha);
  Double_t r=TMath::Sqrt(1 - p[1]*p[1]);
  p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];

  return kTRUE;
}
