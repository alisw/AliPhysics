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

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliTOFtrack class                                                       //
//                                                                         //
// Authors: Bologna-CERN-ITEP-Salerno Group                                //
//                                                                         //
// Description: class for handling ESD extracted tracks for TOF matching.  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "AliESDtrack.h" 
#include "AliTracker.h" 

#include "AliTOFGeometryV5.h"
#include "AliTOFGeometry.h"
#include "AliTOFtrack.h" 

ClassImp(AliTOFtrack)

//_____________________________________________________________________________
AliTOFtrack::AliTOFtrack() : 
  AliKalmanTrack(),
  fSeedInd(-1),
  fSeedLab(-1),
  fTOFgeometry(0)
{
  //
  // Default constructor.
  //
}                                

//_____________________________________________________________________________
AliTOFtrack::AliTOFtrack(const AliTOFtrack& t) : 
  AliKalmanTrack(t),
  fSeedInd(t.fSeedInd),
  fSeedLab(t.fSeedLab),
  fTOFgeometry(new AliTOFGeometryV5()) 
{
  //
  // Copy constructor.
  //
}                                

//_____________________________________________________________________________
AliTOFtrack::AliTOFtrack(const AliESDtrack& t) :
  AliKalmanTrack(), 
  fSeedInd(-1),
  fSeedLab(-1),
  fTOFgeometry(new AliTOFGeometryV5()) 
{
  //
  // Constructor from AliESDtrack
  //
  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetMass(t.GetMass());

  Set(t.GetX(),t.GetAlpha(),t.GetParameter(),t.GetCovariance());

  if ((t.GetStatus()&AliESDtrack::kTIME) == 0) return;
  StartTimeIntegral();
  Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());

}              

//____________________________________________________________________________
AliTOFtrack& AliTOFtrack::operator=(const AliTOFtrack &source)
{
  // ass. op.

  this->fTOFgeometry=source.fTOFgeometry;
  return *this;

}

//_____________________________________________________________________________
Bool_t AliTOFtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho)
{
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0

  if (xk == GetX()) return kTRUE;
  
  Double_t oldX=GetX(), oldY=GetY(), oldZ=GetZ();

  Double_t bz=GetBz();
  if (!AliExternalTrackParam::PropagateTo(xk,bz)) return kFALSE;

  Double_t d = TMath::Sqrt((GetX()-oldX)*(GetX()-oldX) + 
                           (GetY()-oldY)*(GetY()-oldY) + 
                           (GetZ()-oldZ)*(GetZ()-oldZ));
  if (IsStartedTimeIntegral() && GetX()>oldX) AddTimeStep(d);

  if (!AliExternalTrackParam::CorrectForMaterial(d*rho/x0,x0,GetMass())) 
     return kFALSE;

  /*
  //Energy losses************************
  if((5940*beta2/(1-beta2+1e-10) - beta2) < 0){return 0;}

  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2+1e-10)) - beta2)*d*rho;
  //
  // suspicious part - think about it ?
  Double_t kinE =  TMath::Sqrt(p2);
  if (dE>0.8*kinE) dE = 0.8*kinE;  //      
  if (dE<0)        dE = 0.0;       // not valid region for Bethe bloch 
  */

  return kTRUE;            
}     

//_____________________________________________________________________________
Bool_t AliTOFtrack::PropagateToInnerTOF(Bool_t holes)
{
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0


  Double_t ymax=fTOFgeometry->RinTOF()*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
  Bool_t skip = kFALSE;
  Double_t y=GetYat(fTOFgeometry->RinTOF(),skip);
  if (skip) {
    return kFALSE;
  }
  if (y > ymax) {
    if (!Rotate(AliTOFGeometry::GetAlpha())) {
      return kFALSE;
    }
  } else if (y <-ymax) {
    if (!Rotate(-AliTOFGeometry::GetAlpha())) {
      return kFALSE;
    }
  }
  
  
  Double_t x = GetX();
  Int_t nsteps=Int_t((370.-x)/0.5); // 0.5 cm Steps
  for (Int_t istep=0;istep<nsteps;istep++){
    Float_t xp = x+istep*0.5; 
    Double_t param[2];  
    GetPropagationParameters(holes,param);  
    PropagateTo(xp,param[0],param[1]);
    
  }
  
  if(!PropagateTo(fTOFgeometry->RinTOF()))return 0;
  
  return kTRUE;
  
}     

//_________________________________________________________________________
Double_t AliTOFtrack::GetYat(Double_t xk, Bool_t & skip) const {     
//-----------------------------------------------------------------
// This function calculates the Y-coordinate of a track at the plane x=xk.
// Needed for matching with the TOF (I.Belikov)
//-----------------------------------------------------------------
     Double_t y=0.;
     skip=(!GetYAt(xk,GetBz(),y));
     return y;
}

//_____________________________________________________________________________
Int_t AliTOFtrack::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliTOFtrack *t=(AliTOFtrack*)o;
  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

Double_t AliTOFtrack::GetBz() const {
  //
  // returns Bz component of the magnetic field (kG)
  //
  if (AliTracker::UniformField()) return AliTracker::GetBz();
  Double_t r[3]; GetXYZ(r);
  return AliTracker::GetBz(r);
}

//_____________________________________________________________________________
void AliTOFtrack::GetPropagationParameters(Bool_t holes, Double_t *param) {

 //Get average medium density, x0 while propagating the track

  //For TRD holes description

  Double_t thetamin = (90.-31.1) * TMath::Pi()/180.;
  Double_t thetamax = (90.+31.1) * TMath::Pi()/180.;

  Double_t zmin = -55.;
  Double_t zmax =  55.;

  // Detector inner/outer radii
  Double_t rTPC    = 261.53;
  Double_t rTPCTRD = 294.5;
  Double_t rTRD    = 369.1;

  // Medium parameters
  Double_t x0TPC = 40.;
  Double_t rhoTPC =0.06124;

  Double_t x0Air = 36.66;
  Double_t rhoAir =1.2931e-3;

  Double_t x0TRD = 171.7;
  Double_t rhoTRD =0.33;

  Int_t isec = GetSector();
  Double_t r[3]; GetXYZ(r);
  Float_t thetatr = TMath::ATan2(TMath::Sqrt(r[0]*r[0]+r[1]*r[1]),r[2]);

  if(holes){
    if (isec == 0 || isec == 1 || isec == 2 ) {
      if( thetatr>=thetamin && thetatr<=thetamax){ 
	x0TRD= x0Air;
	rhoTRD = rhoAir;
      }
    }
    if (isec == 11 || isec == 12 || isec == 13 || isec == 14 || isec == 15 ) {
      if( r[2]>=zmin && r[2]<=zmax){ 
	x0TRD= x0Air;
	rhoTRD = rhoAir;
      }
    }
  }

  if(GetX() <= rTPC)
    {param[0]=x0TPC;param[1]=rhoTPC;}
  else if(GetX() > rTPC &&  GetX() < rTPCTRD)
    {param[0]=x0Air;param[1]=rhoAir;}
  else if(GetX() >= rTPCTRD &&  GetX() < rTRD)
    {param[0]=x0TRD;param[1]=rhoTRD;}
  else if(GetX() >= rTRD )
    {param[0]=x0Air;param[1]=rhoAir;}
}
