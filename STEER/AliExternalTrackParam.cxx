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

#include <TMath.h>
#include <TVector3.h>

#include "AliExternalTrackParam.h"
#include "AliKalmanTrack.h"
#include "AliTrackReference.h"

ClassImp(AliExternalTrackParam)

const AliMagF *AliExternalTrackParam::fgkFieldMap=0;
Double_t AliExternalTrackParam::fgConvConst=0.;


//_____________________________________________________________________________
AliExternalTrackParam::AliExternalTrackParam() :
  fX(0),
  fAlpha(0),
  fLocalConvConst(0)
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
  fX(x),
  fAlpha(alpha),
  fLocalConvConst(0)
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


AliExternalTrackParam * AliExternalTrackParam::MakeTrack(const AliTrackReference *ref, Double_t mass)
{
  //
  // Make dummy track from the track reference 
  // negative mass means opposite charge 
  //
  Double_t xx[5];
  Double_t cc[15];
  for (Int_t i=0;i<15;i++) cc[i]=0;
  Double_t x = ref->X(), y = ref->Y(), z = ref->Z();
  Double_t alpha = TMath::ATan2(y,x);
  Double_t xr = TMath::Sqrt(x*x+y*y);
  xx[0] = 0;
  xx[1] = z;
  xx[3] = ref->Pz()/ref->Pt();
  Float_t b[3];
  Float_t xyz[3]={x,y,z};
  Float_t convConst = 0;
  (AliKalmanTrack::GetFieldMap())->Field(xyz,b);
  convConst=1000/0.299792458/(1e-13 - b[2]);
  xx[4] = 1./(convConst*ref->Pt()); // curvature rpresentation
  if (mass<0) xx[4]*=-1.;  // negative mass - negative direction
  Double_t alphap = TMath::ATan2(ref->Py(),ref->Px())-alpha;
  if (alphap> TMath::Pi()) alphap-=TMath::Pi();
  if (alphap<-TMath::Pi()) alphap+=TMath::Pi();
  xx[2] = TMath::Sin(alphap);
  xx[4]*=convConst;   // 1/pt representation 
  //  AliExternalTrackParam * track = new  AliExternalTrackParam(xx,cc,xr,alpha);
  AliExternalTrackParam * track = new  AliExternalTrackParam(xr,alpha,xx,cc);
  track->fMass = TMath::Abs(mass);
  //track->StartTimeIntegral();  
  track->SaveLocalConvConst(); 
  return track;
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
Bool_t AliExternalTrackParam::PropagateTo(Double_t xk, Double_t x0, Double_t rho)
{
  //
  // Propagate the track parameters to the given x coordinate assuming vacuum.
  // If length is not NULL, the change of track length is added to it.
  //
  
  Double_t lcc=GetLocalConvConst();  
  Double_t cur = fParam[4]/lcc;
  Double_t x1=fX, x2=xk, dx=x2-x1;
  Double_t f1=fParam[2], f2=f1 + cur*dx;
  if (TMath::Abs(f2) >= 0.98) {
    // MI change  - don't propagate highly inclined tracks
    //              covariance matrix distorted
    return kFALSE;
  }

  // old position [SR, GSI, 17.02.2003]
  Double_t oldX = fX, oldY = fParam[0], oldZ = fParam[1];
  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);  
  fParam[0] += dx*(f1+f2)/(r1+r2);
  fParam[1] += dx*(f1+f2)/(f1*r2 + f2*r1)*fParam[3];
  fParam[2] += dx*cur;
  // transform error matrix to the curvature
  fCovar[10]/=lcc;
  fCovar[11]/=lcc;
  fCovar[12]/=lcc;
  fCovar[13]/=lcc;
  fCovar[14]/=lcc*lcc;

  //f = F - 1
  
  Double_t f02=    dx/(r1*r1*r1);
  Double_t f04=0.5*dx*dx/(r1*r1*r1);
  Double_t f12=    dx*fParam[3]*f1/(r1*r1*r1);
  Double_t f14=0.5*dx*dx*fParam[3]*f1/(r1*r1*r1);
  Double_t f13=    dx/r1;
  Double_t f24=    dx; 
  
  //b = C*ft
  Double_t b00=f02*fCovar[3] + f04*fCovar[10], b01=f12*fCovar[3] + f14*fCovar[10] + f13*fCovar[6];
  Double_t b02=f24*fCovar[10];
  Double_t b10=f02*fCovar[4] + f04*fCovar[11], b11=f12*fCovar[4] + f14*fCovar[11] + f13*fCovar[7];
  Double_t b12=f24*fCovar[11];
  Double_t b20=f02*fCovar[5] + f04*fCovar[12], b21=f12*fCovar[5] + f14*fCovar[12] + f13*fCovar[8];
  Double_t b22=f24*fCovar[12];
  Double_t b40=f02*fCovar[12] + f04*fCovar[14], b41=f12*fCovar[12] + f14*fCovar[14] + f13*fCovar[13];
  Double_t b42=f24*fCovar[14];
  Double_t b30=f02*fCovar[8] + f04*fCovar[13], b31=f12*fCovar[8] + f14*fCovar[13] + f13*fCovar[9];
  Double_t b32=f24*fCovar[13];
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f24*b42;

  //F*C*Ft = C + (b + bt + a)
  fCovar[0] += b00 + b00 + a00;
  fCovar[1] += b10 + b01 + a01; 
  fCovar[2] += b11 + b11 + a11;
  fCovar[3] += b20 + b02 + a02;
  fCovar[4] += b21 + b12 + a12;
  fCovar[5] += b22 + b22 + a22;
  fCovar[6] += b30;
  fCovar[7] += b31; 
  fCovar[8] += b32;
  fCovar[10] += b40;
  fCovar[11] += b41;
  fCovar[12] += b42;

  fX=x2;

  //Change of the magnetic field *************
  SaveLocalConvConst();
  // transform back error matrix from curvature to the 1/pt
  fCovar[10]*=lcc;
  fCovar[11]*=lcc;
  fCovar[12]*=lcc;
  fCovar[13]*=lcc;
  fCovar[14]*=lcc*lcc;

  Double_t dist = TMath::Sqrt((fX-oldX)*(fX-oldX)+(fParam[0]-oldY)*(fParam[0]-oldY)+
			      (fParam[1]-oldZ)*(fParam[1]-oldZ));
  if (!CorrectForMaterial(dist,x0,rho)) return 0;

  // Integrated Time [SR, GSI, 17.02.2003]
 //  if (IsStartedTimeIntegral() && fX>oldX) {
//     Double_t l2 = (fX-oldX)*(fX-oldX)+(fParam[0]-oldY)*(fParam[0]-oldY)+
//                   (fParam[1]-oldZ)*(fParam[1]-oldZ);
//     AddTimeStep(TMath::Sqrt(l2));
//   }
  //

  return kTRUE;
}

Bool_t     AliExternalTrackParam::PropagateToDCA(Double_t xd, Double_t yd,  Double_t x0, Double_t rho){
  //
  // Propagate the track parameters to the nearest point of given xv yv coordinate
  //
  Double_t a=fAlpha;
  Double_t cs=TMath::Cos(a),sn=TMath::Sin(a);

  Double_t xv= xd*cs + yd*sn;
  Double_t yv=-xd*sn + yd*cs;   // vertex position in local frame
  //  
  Double_t c=fParam[4]/GetLocalConvConst(), snp=fParam[2];
  //
  Double_t x=fX, y=fParam[1];
  Double_t tgfv=-(c*(x-xv)-snp)/(c*(y-yv) + TMath::Sqrt(1.-snp*snp));
  Double_t fv=TMath::ATan(tgfv);
  cs=TMath::Cos(fv); sn=TMath::Sin(fv);
  x = xv*cs + yv*sn;
  yv=-xv*sn + yv*cs; xv=x;
  RotateTo(fv+a);
  return PropagateTo(xv,x0,rho);  
}


//_____________________________________________________________________________
Bool_t AliExternalTrackParam::RotateTo(Double_t alp)
{
  // Rotate the reference axis for the parametrisation to the given angle.
  //
  Double_t  x=fX;
  Double_t p0=fParam[0];
  //
  if      (alp < -TMath::Pi()) alp += 2*TMath::Pi();
  else if (alp >= TMath::Pi()) alp -= 2*TMath::Pi();
  Double_t ca=TMath::Cos(alp-fAlpha), sa=TMath::Sin(alp-fAlpha);
  Double_t sf=fParam[2], cf=TMath::Sqrt(1.- fParam[2]*fParam[2]);
  // **** rotation **********************
  
  fAlpha = alp;
  fX =  x*ca + p0*sa;
  fParam[0]= -x*sa + p0*ca;
  fParam[2]=  sf*ca - cf*sa;
  Double_t rr=(ca+sf/cf*sa);  
  //
  fCovar[0] *= (ca*ca);
  fCovar[1] *= ca; 
  fCovar[3] *= ca*rr;
  fCovar[6] *= ca;
  fCovar[10] *= ca;
  fCovar[4] *= rr;
  fCovar[5] *= rr*rr;
  fCovar[7] *= rr;
  fCovar[11] *= rr;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::CorrectForMaterial(Double_t d, Double_t x0, Double_t rho)
{
  //
  // Take into account material effects assuming:
  // x0  - mean rad length
  // rho - mean density

  //
  // multiple scattering
  //
  Double_t p2=(1.+ fParam[3]*fParam[3])/(fParam[4]*fParam[4]);
  Double_t beta2=p2/(p2 + fMass*fMass);
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho;
  //
  fCovar[5] += theta2*(1.- fParam[2]*fParam[2])*(1. + fParam[3]*fParam[3]);
  fCovar[9] += theta2*(1. + fParam[3]*fParam[3])*(1. + fParam[3]*fParam[3]);
  fCovar[13] += theta2*fParam[3]*fParam[4]*(1. + fParam[3]*fParam[3]);
  fCovar[14] += theta2*fParam[3]*fParam[4]*fParam[3]*fParam[4];
  
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2+1e-10)) - beta2)*d*rho;  
  fParam[4] *=(1.- TMath::Sqrt(p2+fMass*fMass)/p2*dE);
  //
  Double_t sigmade = 0.02*TMath::Sqrt(TMath::Abs(dE));   // energy loss fluctuation 
  Double_t sigmac2 = sigmade*sigmade*fParam[4]*fParam[4]*(p2+fMass*fMass)/(p2*p2);
  fCovar[14] += sigmac2;
  //
  //
  


  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliExternalTrackParam::GetProlongationAt(Double_t xk, 
						Double_t& y, 
						Double_t& z) const
{
  //
  // Get the local y and z coordinates at the given x value
  //
  Double_t lcc=GetLocalConvConst();  
  Double_t cur = fParam[4]/lcc;
  Double_t x1=fX, x2=xk, dx=x2-x1; 
  Double_t f1=fParam[2], f2=f1 + cur*dx;
  //  
  if (TMath::Abs(f2) >= 0.98) {
    // MI change  - don't propagate highly inclined tracks
    //              covariance matrix distorted
    return kFALSE;
  }
  Double_t r1=sqrt(1.- f1*f1), r2=sqrt(1.- f2*f2);
  y = fParam[0] + dx*(f1+f2)/(r1+r2);
  z = fParam[1] + dx*(f1+f2)/(f1*r2 + f2*r1)*fParam[3];
  return kTRUE;
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


// //_____________________________________________________________________________
// Double_t AliExternalTrackParam::GetPredictedChi2(const AliCluster* /*cluster*/)
// {
// // calculate the chi2 contribution of the given cluster
// //
// // NOT IMPLEMENTED for this class

//   return -1;
// }

// //_____________________________________________________________________________
// Bool_t AliExternalTrackParam::Update(const AliCluster* /*cluster*/)
// {
// // update the track parameters using the position and error 
// // of the given cluster
// //
// // NOT IMPLEMENTED for this class

//   return kFALSE;
// }


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
