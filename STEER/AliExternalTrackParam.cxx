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
// track parameters in "external" format                                     //
//                                                                           //
// The track parameters are:                                                //
// - local y coordinate                                                      //
// - local z coordinate                                                      //
// - sin of azimutal angle                                                   //
// - tan of dip angle                                                        //
// - charge/pt                                                               //
// The parametrisation is given at the local x coordinate fX and the         //
// azimuthal angle fAlpha.                                                   //
//                                                                           //
// The external parametrisation can be used to exchange track parameters     //
// between different detectors.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliExternalTrackParam.h"
#include "AliKalmanTrack.h"

ClassImp(AliExternalTrackParam)


//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam() :
  AliTrackParam(),
  fX(0),
  fAlpha(0)
{
  //
  // default constructor
  //
  for (Int_t i = 0; i < 5; i++) fParam[i] = 0;
  for (Int_t i = 0; i < 15; i++) fCovar[i] = 0;
}

//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam(Double_t x, Double_t alpha, 
					     const Double_t param[5], 
					     const Double_t covar[15]) :
  AliTrackParam(),
  fX(x),
  fAlpha(alpha)
{
  //
  // create external track parameters from given arguments
  //
  for (Int_t i = 0; i < 5; i++) fParam[i] = param[i];
  for (Int_t i = 0; i < 15; i++) fCovar[i] = covar[i];
}

//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam(const AliKalmanTrack& track) :
  fX(0),
  fAlpha(track.GetAlpha())
{
  //
  //
  track.GetExternalParameters(fX,fParam);
  track.GetExternalCovariance(fCovar);
}

//_____________________________________________________________________________
const Double_t* AliExternalTrackParam::GetParameter() const
{
// get a pointer to the array of track parameters

  return fParam;
}

//_____________________________________________________________________________
const Double_t* AliExternalTrackParam::GetCovariance() const
{
// get a pointer to the array of the track parameter covariance matrix

  return fCovar;
}

//_____________________________________________________________________________
AliExternalTrackParam* AliExternalTrackParam::CreateExternalParam() const
{
// copy this instance

  return new AliExternalTrackParam(fX, fAlpha, fParam, fCovar);
}

//_____________________________________________________________________________
void AliExternalTrackParam::ResetCovariance(Double_t factor,
					    Bool_t clearOffDiagonal)
{
// reset the covariance matrix ("forget" track history)

  Int_t k = 0;
  for (Int_t i = 0; i < 5; i++) {
    for (Int_t j = 0; j < i; j++) {  // off diagonal elements
      if (clearOffDiagonal) {
	fCovar[k++] = 0;
      } else {
	fCovar[k++] *= factor;
      }
    }
    fCovar[k++] *= factor;     // diagonal elements
  }
}


//_____________________________________________________________________________
Bool_t AliExternalTrackParam::PropagateTo(Double_t /*x*/, Double_t* /*length*/)
{
// Propagate the track parameters to the given x coordinate assuming vacuum.
// If length is not NULL, the change of track length is added to it.
//
// NOT IMPLEMENTED for this class

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::RotateTo(Double_t /*alpha*/)
{
// Rotate the reference axis for the parametrisation to the given angle.
//
// NOT IMPLEMENTED for this class

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::CorrectForMaterial(Double_t /*dAngle*/, 
						 Double_t /*dPrel*/)
{
// Take into account material effects assuming:
//   dAngle2: mean multiple scattering angle in rad squared
//   dPrel  : mean relative momentum gain (if > 0) or loss (if < 0)
//
// NOT IMPLEMENTED for this class

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::GetProlongationAt(Double_t /*x*/, 
						Double_t& /*y*/, 
						Double_t& /*z*/) const
{
// Get the local y and z coordinates at the given x value
//
// NOT IMPLEMENTED for this class

  return kFALSE;
}

//_____________________________________________________________________________
Double_t AliExternalTrackParam::GetXAtVertex(Double_t /*x*/, 
					     Double_t /*y*/) const
{
// Get the x coordinate at the given vertex (x,y)
//
// NOT IMPLEMENTED for this class

  return 0;
}


//_____________________________________________________________________________
Double_t AliExternalTrackParam::GetPredictedChi2(const AliCluster* /*cluster*/)
{
// calculate the chi2 contribution of the given cluster
//
// NOT IMPLEMENTED for this class

  return -1;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::Update(const AliCluster* /*cluster*/)
{
// update the track parameters using the position and error 
// of the given cluster
//
// NOT IMPLEMENTED for this class

  return kFALSE;
}


//_____________________________________________________________________________
Double_t AliExternalTrackParam::SigmaPhi() const
{
// get the error of the azimuthal angle

  return TMath::Sqrt(TMath::Abs(fCovar[5] / (1. - fParam[2]*fParam[2])));
}

//_____________________________________________________________________________
Double_t AliExternalTrackParam::SigmaTheta() const
{
// get the error of the polar angle

  return TMath::Sqrt(TMath::Abs(fCovar[9])) / (1. + fParam[3]*fParam[3]);
}

//_____________________________________________________________________________
Double_t AliExternalTrackParam::SigmaPt() const
{
// get the error of the transversal component of the momentum

  return TMath::Sqrt(fCovar[14]) / TMath::Abs(fParam[4]);
}

//_____________________________________________________________________________
TVector3 AliExternalTrackParam::Momentum() const
{
// get the momentum vector

  Double_t phi = TMath::ASin(fParam[2]) + fAlpha;
  Double_t pt = 1. / TMath::Abs(fParam[4]);
  return TVector3(pt * TMath::Cos(phi), 
		  pt * TMath::Sin(phi), 
		  pt * fParam[3]);
}

//_____________________________________________________________________________
TVector3 AliExternalTrackParam::Position() const
{
// get the current spatial position in global coordinates

  return TVector3(fX * TMath::Cos(fAlpha) - fParam[0] * TMath::Sin(fAlpha),
		  fX * TMath::Sin(fAlpha) + fParam[0] * TMath::Cos(fAlpha),
		  fParam[1]);
}


//_____________________________________________________________________________
void AliExternalTrackParam::Print(Option_t* /*option*/) const
{
// print the parameters and the covariance matrix

  printf("AliExternalTrackParam: x = %-12g  alpha = %-12g\n", fX, fAlpha);
  printf("  parameters: %12g %12g %12g %12g %12g\n",
	 fParam[0], fParam[1], fParam[2], fParam[3], fParam[4]);
  printf("  covariance: %12g\n", fCovar[0]);
  printf("              %12g %12g\n", fCovar[1], fCovar[2]);
  printf("              %12g %12g %12g\n", fCovar[3], fCovar[4], fCovar[5]);
  printf("              %12g %12g %12g %12g\n", 
	 fCovar[6], fCovar[7], fCovar[8], fCovar[9]);
  printf("              %12g %12g %12g %12g %12g\n", 
	 fCovar[10], fCovar[11], fCovar[12], fCovar[13], fCovar[14]);
}
