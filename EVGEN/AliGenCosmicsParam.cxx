/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


// Generator for muons according to kinematic parametrizations at ALICE
// (not at the surface).
// Origin: andrea.dainese@lnl.infn.it

#include <TParticle.h>
#include <TF1.h>

#include "AliRun.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliGenCosmicsParam.h"

ClassImp(AliGenCosmicsParam)

//-----------------------------------------------------------------------------
AliGenCosmicsParam::AliGenCosmicsParam():
AliGenerator(),
fParamMI(kFALSE),
fParamACORDE(kFALSE),
fYOrigin(600.),
fMaxAngleWRTVertical(-99.),
fBkG(0.),
fTPC(kFALSE),
fITS(kFALSE),
fSPDouter(kFALSE),
fSPDinner(kFALSE)
{
  //
  // Default constructor
  //
}
//-----------------------------------------------------------------------------
void AliGenCosmicsParam::Generate()
{
  //
  // Generate one muon
  //
  
  //
  Float_t origin[3];
  Float_t p[3];
  Int_t nt;
  Double_t ptot=0,pt=0,angleWRTVertical=0;
  Bool_t okMom=kFALSE,okAngle=kFALSE;
  //

  // mu+ or mu-
  Int_t ipart=13;
  if(gRandom->Rndm()<0.5) ipart *= -1;

  if(fParamACORDE) { // extracted from AliGenACORDE events
    // sample total momentum only once (to speed up)
    TF1 *dNdpACORDE = new TF1("dNdpACORDE","x/(1.+(x/12.8)*(x/12.8))^1.96",fPMin,fPMax);
    ptot = (Double_t)dNdpACORDE->GetRandom();
    delete dNdpACORDE;
    dNdpACORDE = 0;
  }

  Int_t trials=0;

  while(1) {
    trials++;
    // origin
    origin[0]  = fYOrigin*TMath::Tan(fMaxAngleWRTVertical)*(-1.+2.*gRandom->Rndm());
    origin[1]  = fYOrigin;
    origin[2]  = fYOrigin*TMath::Tan(fMaxAngleWRTVertical)*(-1.+2.*gRandom->Rndm());

    // momentum
    while(1) {
      okMom=kFALSE; okAngle=kFALSE;

      if(fParamMI) { // parametrization by M.Ivanov of LEP cosmics data
	Float_t	pref  = 1. + gRandom->Exp(30.);
	p[1] = -pref; 
	p[0] = gRandom->Gaus(0.0,0.2)*pref;
	p[2] = gRandom->Gaus(0.0,0.2)*pref;
	if(gRandom->Rndm()>0.9) {
	  p[0] = gRandom->Gaus(0.0,0.4)*pref;
	  p[2] = gRandom->Gaus(0.0,0.4)*pref;
	}
	ptot=TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	pt=TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
      } else if(fParamACORDE) { // extracted from AliGenACORDE events
	Float_t theta,phi;
	while(1) {
	  theta = gRandom->Gaus(0.5*TMath::Pi(),0.42);
	  if(TMath::Abs(theta-0.5*TMath::Pi())<fMaxAngleWRTVertical) break;
	}
	while(1) {
	  phi = gRandom->Gaus(-0.5*TMath::Pi(),0.42);
	  if(TMath::Abs(phi+0.5*TMath::Pi())<fMaxAngleWRTVertical) break;
	}
	pt = ptot*TMath::Sin(theta);
	p[0] = pt*TMath::Cos(phi); 
	p[1] = pt*TMath::Sin(phi); 
	p[2] = ptot*TMath::Cos(theta);
      } else {
	AliFatal("Parametrization not set: use SetParamMI or SetParamACORDE");
      }

      
      // check kinematic cuts
      if(TestBit(kMomentumRange)) {
	if(ptot>fPMin && ptot<fPMax) okMom=kTRUE;
      }

      angleWRTVertical=TMath::ACos(TMath::Abs(p[1])/ptot); // acos(|py|/ptot)
      if(angleWRTVertical<fMaxAngleWRTVertical) okAngle=kTRUE;
      
      if(okAngle&&okMom) break;
    }

    // acceptance
    if(!fTPC && !fITS && !fSPDinner && !fSPDouter) break;
    if(fTPC) if(IntersectCylinder(250.,250.,ipart,origin,p)) break;
    if(fITS) if(IntersectCylinder(50.,50.,ipart,origin,p)) break;
    if(fSPDouter) if(IntersectCylinder(6.5,14.,ipart,origin,p)) break;
    if(fSPDinner) if(IntersectCylinder(3.5,14.0,ipart,origin,p)) break;
  }

  Float_t polarization[3]= {0,0,0};
  PushTrack(fTrackIt,-1,ipart,p,origin,polarization,0,kPPrimary,nt);

  //printf("TRIALS %d\n",trials);
  

  return;
}
//-----------------------------------------------------------------------------
void AliGenCosmicsParam::Init()
{
  // 
  // Initialisation, check consistency of selected ranges
  //
  if(TestBit(kPtRange)) 
    AliFatal("You cannot set the pt range for this generator! Only momentum range");
  if(fPMin<8.) { 
    fPMin=8.; 
    if(TestBit(kMomentumRange)) 
      AliWarning("Minimum momentum cannot be < 8 GeV/c"); 
  }
  if(fMaxAngleWRTVertical<0.) 
    AliFatal("You must use SetMaxAngleWRTVertical() instead of SetThetaRange(), SetPhiRange()");

  printf("************ AliGenCosmicsParam ****************\n");
  printf("***** Muons generated at Y = %f cm\n",fYOrigin);
  printf("************************************************\n");

  return;
}
//-----------------------------------------------------------------------------
Bool_t AliGenCosmicsParam::IntersectCylinder(Float_t r,Float_t z,Int_t pdg,
					     Float_t o[3],Float_t p[3]) const
{
  //
  // Intersection between muon and cylinder [-z,+z] with radius r
  //

  Float_t en = TMath::Sqrt(0.105*0.105+p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  TParticle part(pdg,0,0,0,0,0,p[0],p[1],p[2],en,o[0],o[1],o[2],0);
  AliESDtrack track(&part);
  Double_t pos[3]={0.,0.,0.},sigma[3]={0.,0.,0.};
  AliESDVertex origin(pos,sigma);

  track.RelateToVertex(&origin,fBkG,10000.);

  Float_t d0z0[2],covd0z0[3];
  track.GetImpactParameters(d0z0,covd0z0);

  // check rphi 
  if(TMath::Abs(d0z0[0])>r) return kFALSE;
  // check z
  if(TMath::Abs(d0z0[1])>z) return kFALSE;

  /*
    if(TMath::Abs(fB)<0.01) {  // NO FIELD
    Float_t drphi = TMath::Abs(o[1]-p[1]/p[0]*o[0])/
    TMath::Sqrt(p[1]*p[1]/p[0]/p[0]+1.);
    if(drphi>r) return kFALSE;
    Float_t dz = o[2]-p[2]/p[0]*o[0]+p[2]/p[0]*
    (p[1]*p[1]/p[0]/p[0]*o[0]-p[1]/p[0]*o[1])/(1.+p[1]*p[1]/p[0]/p[0]);
    if(TMath::Abs(dz)>z) return kFALSE;
    }
  */

  return kTRUE;
}
//-----------------------------------------------------------------------------
