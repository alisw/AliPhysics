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
//                Implementation of the AliKalmanTrack class
//
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliKalmanTrack.h"
#include "AliCluster.h"
#include <TMath.h>
#include <iostream.h>

ClassImp(AliKalmanTrack)

//_____________________________________________________________________________
AliKalmanTrack::AliKalmanTrack(const AliKalmanTrack& t) {
  //-----------------------------------------------------------------
  // This is a copy constructor.
  //-----------------------------------------------------------------
  fLab=t.fLab;

  fP0=t.fP0; fP1=t.fP1; fP2=t.fP2; fP3=t.fP3; fP4=t.fP4;

  fC00=t.fC00;
  fC10=t.fC10;  fC11=t.fC11;
  fC20=t.fC20;  fC21=t.fC21;  fC22=t.fC22;
  fC30=t.fC30;  fC31=t.fC31;  fC32=t.fC32;  fC33=t.fC33;
  fC40=t.fC40;  fC41=t.fC41;  fC42=t.fC42;  fC43=t.fC43;  fC44=t.fC44;

  fChi2=t.fChi2;
  fN=t.fN;
}

//_____________________________________________________________________________
Int_t AliKalmanTrack::Compare(TObject *o) {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliKalmanTrack *t=(AliKalmanTrack*)o;
  Double_t co=TMath::Abs(t->GetPt());
  Double_t c =TMath::Abs(GetPt());
  if (c<co) return 1;
  else if (c>co) return -1;
  return 0;
}

//_____________________________________________________________________________
Double_t AliKalmanTrack::GetPredictedChi2(const AliCluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    if (fN>4) cerr<<fN<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}

//_____________________________________________________________________________
void AliKalmanTrack::GetCovariance(Double_t cc[15]) const {
  // return covariance maxtrix
  cc[0 ]=fC00;
  cc[1 ]=fC10;  cc[2 ]=fC11;
  cc[3 ]=fC20;  cc[4 ]=fC21;  cc[5 ]=fC22;
  cc[6 ]=fC30;  cc[7 ]=fC31;  cc[8 ]=fC32;  cc[9 ]=fC33;
  cc[10]=fC40;  cc[11]=fC41;  cc[12]=fC42;  cc[13]=fC43;  cc[14]=fC44;
}



