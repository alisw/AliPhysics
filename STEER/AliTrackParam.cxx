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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// base class for track parameters                                           //
//                                                                           //
// The idea of having a separate class for the track parametrisation and     //
// not including it in the track classes itself is to not duplicate code.    //
// Track classes which use the same parametrisation can include a track      //
// parameters object for their common parametrisation. The code, e.g. for    //
// the propagation, does not need to be duplicated.                          //
//                                                                           //
// The AliTrackParam and its derived classes                                 //
// - know the current track parameters and their covariance matrix as well   //
//   as the local x coordinate and azimuthal angle alpha of the current      //
//   parametrisation                                                         //
// - can rotate their local coordinate system                                //
// - can create a parametrisation in external format from itself             //
// - can propagate through material or vacuum                                //
// - can calculate a chi^2 for a cluster                                     //
// - can update the parameters using the position and covariance of a        //
//   cluster                                                                 //
//                                                                           //
// In addition some methods to get quantities useful for analysis, like      //
// the momentum, are implemented.                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliTrackParam.h"
#include "AliExternalTrackParam.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliMagF.h"


Double_t FastMath::fgFastAsin[20000];

FastMath::FastMath()
{
  for (Int_t i=0;i<10000;i++){
    fgFastAsin[2*i] = TMath::ASin(i/10000.);
    fgFastAsin[2*i+1] = (TMath::ASin((i+1)/10000.)-fgFastAsin[2*i]);
  }
}

Double_t FastMath::FastAsin(Double_t x)
{
  if (x>0){
    Int_t index = int(x*10000);
    return fgFastAsin[2*index]+(x*10000.-index)*fgFastAsin[2*index+1];
  }
  x*=-1;
  Int_t index = int(x*10000);
  return -(fgFastAsin[2*index]+(x*10000.-index)*fgFastAsin[2*index+1]);
}

FastMath gFastMath;



ClassImp(AliTrackParam)




//_____________________________________________________________________________
Bool_t AliTrackParam::RotateAndPropagateTo(Double_t alpha, Double_t x, 
					   Double_t* length)
{
// Rotate the reference axis for the parametrisation to the given angle and
// propagate the track parameters to the given x coordinate assuming vacuum.
// If length is not NULL, the change of track length is added to it.

  if (!RotateTo(alpha)) return kFALSE;
  if (!PropagateTo(x, length)) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Double_t AliTrackParam::GetDsdx() const
{
// get the change of track length s per step in x: ds/dx

  TVector3 x(TMath::Cos(Alpha()), TMath::Sin(Alpha()), 0);
  TVector3 p = Momentum();
  Double_t xp = x*p;
  if (xp == 0) return 1.E6; 
  return p.Mag() / xp;
}


//_____________________________________________________________________________
Double_t AliTrackParam::Phi() const
{
// get the azimuthal angre

  return Momentum().Phi();
}

//_____________________________________________________________________________
Double_t AliTrackParam::SigmaPhi() const
{
// get the error of the azimuthal angle

  AliExternalTrackParam* param = CreateExternalParam();
  Double_t result = param->SigmaPhi();
  delete param;
  return result;
}

//_____________________________________________________________________________
Double_t AliTrackParam::Theta() const
{
// the the polar angle

  return Momentum().Theta();
}

//_____________________________________________________________________________
Double_t AliTrackParam::SigmaTheta() const
{
// get the error of the polar angle

  AliExternalTrackParam* param = CreateExternalParam();
  Double_t result = param->SigmaTheta();
  delete param;
  return result;
}

//_____________________________________________________________________________
Double_t AliTrackParam::Eta() const
{
// get the pseudorapidity

  return Momentum().Eta();
}

//_____________________________________________________________________________
Double_t AliTrackParam::Px() const
{
// get the x component of the momentum

  return Momentum().Px();
}

//_____________________________________________________________________________
Double_t AliTrackParam::Py() const
{
// get the y component of the momentum

  return Momentum().Py();
}

//_____________________________________________________________________________
Double_t AliTrackParam::Pz() const
{
// get the z component of the momentum

  return Momentum().Pz();
}

//_____________________________________________________________________________
Double_t AliTrackParam::Pt() const
{
// get the transversal component of the momentum

  return Momentum().Pt();
}

//_____________________________________________________________________________
Double_t AliTrackParam::SigmaPt() const
{
// get the error of the transversal component of the momentum

  AliExternalTrackParam* param = CreateExternalParam();
  Double_t result = param->SigmaPt();
  delete param;
  return result;
}

//_____________________________________________________________________________
Double_t AliTrackParam::P() const
{
// get the absolute momentum

  return Momentum().Mag();
}

//_____________________________________________________________________________
TVector3 AliTrackParam::Momentum() const
{
// get the momentum vector

  AliExternalTrackParam* param = CreateExternalParam();
  TVector3 result = param->Momentum();
  delete param;
  return result;
}

//_____________________________________________________________________________
TVector3 AliTrackParam::Position() const
{
// get the current spatial position in global coordinates

  Double_t sinAlpha = TMath::Sin(Alpha());
  Double_t cosAlpha = TMath::Cos(Alpha());
  return TVector3(X()*cosAlpha - Y()*sinAlpha, 
		  X()*sinAlpha + Y()*cosAlpha, 
		  Z());
}

//_____________________________________________________________________________
TVector3 AliTrackParam::PositionAt(Double_t x) const
{
// get the spatial position at x in global coordinates

  Double_t y;
  Double_t z;
  if (!GetProlongationAt(x, y, z)) return TVector3(0,0,0);
  Double_t sinAlpha = TMath::Sin(Alpha());
  Double_t cosAlpha = TMath::Cos(Alpha());
  return TVector3(x*cosAlpha - y*sinAlpha, x*sinAlpha + y*cosAlpha, z);
}

