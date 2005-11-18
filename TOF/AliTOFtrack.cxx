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
////////////////////////////////////////////////////////////////////////
//
// AliTOFtrack class
//
// Authors: Bologna-CERN-ITEP-Salerno Group
//
// Description: class for handling ESD extracted tracks for TOF matching.
/* $Id$ */

#include <Riostream.h>
#include <TObject.h>   
#include "AliLog.h" 
#include "AliTOFtrack.h" 
#include "AliESDtrack.h" 

ClassImp(AliTOFtrack)

//_____________________________________________________________________________
AliTOFtrack::AliTOFtrack(const AliTOFtrack& t) : AliKalmanTrack(t) {
  //
  // Copy constructor.
  //
  
  SetSeedIndex(t.GetSeedIndex());
  SetLabel(t.GetLabel());
  fSeedLab=t.GetSeedLabel();
  SetChi2(t.GetChi2());

  fAlpha=t.fAlpha;
  fX=t.fX;

  fY=t.fY; fZ=t.fZ; fE=t.fE; fT=t.fT; fC=t.fC;

  fCyy=t.fCyy;
  fCzy=t.fCzy;  fCzz=t.fCzz;
  fCey=t.fCey;  fCez=t.fCez;  fCee=t.fCee;
  fCty=t.fCty;  fCtz=t.fCtz;  fCte=t.fCte;  fCtt=t.fCtt;
  fCcy=t.fCcy;  fCcz=t.fCcz;  fCce=t.fCce;  fCct=t.fCct;  fCcc=t.fCcc;  


}                                

//_____________________________________________________________________________
AliTOFtrack::AliTOFtrack(const AliESDtrack& t) 
           :AliKalmanTrack() {
  //
  // Constructor from AliESDtrack
  //

  SetSeedIndex(-1);
  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetMass(t.GetMass());

  fAlpha = t.GetAlpha();
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();
  Double_t x, p[5]; t.GetExternalParameters(x,p);

  fX=x;

  fY=p[0];
  fZ=p[1]; SaveLocalConvConst();
  fT=p[3]; x=GetLocalConvConst();
  fC=p[4]/x;
  fE=fC*fX - p[2];   

  //Conversion of the covariance matrix
  Double_t c[15]; t.GetExternalCovariance(c);

  c[10]/=x; c[11]/=x; c[12]/=x; c[13]/=x; c[14]/=x*x;

  Double_t c22=fX*fX*c[14] - 2*fX*c[12] + c[5];
  Double_t c32=fX*c[13] - c[8];
  Double_t c20=fX*c[10] - c[3], c21=fX*c[11] - c[4], c42=fX*c[14] - c[12];

  fCyy=c[0 ];
  fCzy=c[1 ];   fCzz=c[2 ];
  fCey=c20;     fCez=c21;     fCee=c22;
  fCty=c[6 ];   fCtz=c[7 ];   fCte=c32;   fCtt=c[9 ];
  fCcy=c[10];   fCcz=c[11];   fCce=c42;   fCct=c[13]; fCcc=c[14];  

  if ((t.GetStatus()&AliESDtrack::kTIME) == 0) return;
  StartTimeIntegral();
  Double_t times[10]; t.GetIntegratedTimes(times); SetIntegratedTimes(times);
  SetIntegratedLength(t.GetIntegratedLength());


}              
//____________________________________________________________________________
void AliTOFtrack::GetExternalParameters(Double_t& xr, Double_t x[5]) const {
  //
  // This function returns external TOF track representation
  //
     xr=fX;
     x[0]=GetY();
     x[1]=GetZ();
     x[2]=GetSnp();
     x[3]=GetTgl();
     x[4]=Get1Pt();
}           

//_____________________________________________________________________________
void AliTOFtrack::GetExternalCovariance(Double_t cc[15]) const {
  //
  // This function returns external representation of the covriance matrix.
  //
  Double_t a=GetLocalConvConst();
  Double_t c22=fX*fX*fCcc-2*fX*fCce+fCee;
  Double_t c32=fX*fCct-fCte;
  Double_t c20=fX*fCcy-fCey, c21=fX*fCcz-fCez, c42=fX*fCcc-fCce;

  cc[0 ]=fCyy;
  cc[1 ]=fCzy;   cc[2 ]=fCzz;
  cc[3 ]=c20;    cc[4 ]=c21;    cc[5 ]=c22;
  cc[6 ]=fCty;   cc[7 ]=fCtz;   cc[8 ]=c32;   cc[9 ]=fCtt;
  cc[10]=fCcy*a; cc[11]=fCcz*a; cc[12]=c42*a; cc[13]=fCct*a; cc[14]=fCcc*a*a; 
  
}               
                       

//_____________________________________________________________________________
void AliTOFtrack::GetCovariance(Double_t cc[15]) const {
  //
  // Returns the covariance matrix.
  //

  cc[0]=fCyy;
  cc[1]=fCzy;  cc[2]=fCzz;
  cc[3]=fCey;  cc[4]=fCez;  cc[5]=fCee;
  cc[6]=fCcy;  cc[7]=fCcz;  cc[8]=fCce;  cc[9]=fCcc;
  cc[10]=fCty; cc[11]=fCtz; cc[12]=fCte; cc[13]=fCct; cc[14]=fCtt;
  
}    


//_____________________________________________________________________________
Int_t AliTOFtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho)
{
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0

  if (xk == fX) return 1;
  
  if (TMath::Abs(fC*xk - fE) >= 0.90000) {
    return 0;
  }
  Double_t lcc=GetLocalConvConst();

  // track Length measurement [SR, GSI, 17.02.2003]

  Double_t oldX = fX, oldY = fY, oldZ = fZ;  

  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1, y1=fY, z1=fZ;
  Double_t c1=fC*x1 - fE;
  if((c1*c1) > 1){
    return 0;}
  Double_t r1=sqrt(1.- c1*c1);
  Double_t c2=fC*x2 - fE; 
  if((c2*c2) > 1) {
    return 0;
  }
  Double_t r2=sqrt(1.- c2*c2);

  fY += dx*(c1+c2)/(r1+r2);
  fZ += dx*(c1+c2)/(c1*r2 + c2*r1)*fT;

  //f = F - 1
  Double_t rr=r1+r2, cc=c1+c2, xx=x1+x2;
  Double_t f02=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  Double_t f04= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  Double_t cr=c1*r2+c2*r1;
  Double_t f12=-dx*fT*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  Double_t f13= dx*cc/cr;
  Double_t f14=dx*fT*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);

  //b = C*ft
  Double_t b00=f02*fCey + f04*fCcy, b01=f12*fCey + f14*fCcy + f13*fCty;
  Double_t b10=f02*fCez + f04*fCcz, b11=f12*fCez + f14*fCcz + f13*fCtz;
  Double_t b20=f02*fCee + f04*fCce, b21=f12*fCee + f14*fCce + f13*fCte;
  Double_t b30=f02*fCte + f04*fCct, b31=f12*fCte + f14*fCct + f13*fCtt;
  Double_t b40=f02*fCce + f04*fCcc, b41=f12*fCce + f14*fCcc + f13*fCct;

  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a11=f12*b21+f14*b41+f13*b31;

  //F*C*Ft = C + (a + b + bt)
  fCyy += a00 + 2*b00;
  fCzy += a01 + b01 + b10;
  fCey += b20;
  fCty += b30;
  fCcy += b40;
  fCzz += a11 + 2*b11;
  fCez += b21;
  fCtz += b31;
  fCcz += b41;

  fX=x2;                                                     

  //Change of the magnetic field *************
  SaveLocalConvConst();
  cc=fC;
  fC*=lcc/GetLocalConvConst();
  fE+=fX*(fC-cc);

  //Multiple scattering  ******************
  Double_t d=sqrt((x1-fX)*(x1-fX)+(y1-fY)*(y1-fY)+(z1-fZ)*(z1-fZ));
  Double_t p2=(1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt());
  Double_t beta2=p2/(p2 + GetMass()*GetMass());
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho;

  Double_t ey=fC*fX - fE, ez=fT;
  Double_t xz=fC*ez, zz1=ez*ez+1, xy=fE+ey;
  
  fCee += (2*ey*ez*ez*fE+1-ey*ey+ez*ez+fE*fE*ez*ez)*theta2;
  fCte += ez*zz1*xy*theta2;
  fCtt += zz1*zz1*theta2;
  fCce += xz*ez*xy*theta2;
  fCct += xz*zz1*theta2;
  fCcc += xz*xz*theta2;
  /*
  Double_t dc22 = (1-ey*ey+xz*xz*fX*fX)*theta2;
  Double_t dc32 = (xz*fX*zz1)*theta2;
  Double_t dc33 = (zz1*zz1)*theta2;
  Double_t dc42 = (xz*fX*xz)*theta2;
  Double_t dc43 = (zz1*xz)*theta2;
  Double_t dc44 = (xz*xz)*theta2; 
  fCee += dc22;
  fCte += dc32;
  fCtt += dc33;
  fCce += dc42;
  fCct += dc43;
  fCcc += dc44;
  */
  //Energy losses************************
  if((5940*beta2/(1-beta2+1e-10) - beta2) < 0){return 0;}

  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2+1e-10)) - beta2)*d*rho;
  if (x1 < x2) dE=-dE;
  cc=fC;
  fC*=(1.- sqrt(p2+GetMass()*GetMass())/p2*dE);
  fE+=fX*(fC-cc);    

  // track time measurement [SR, GSI 17.02.2002]
  if (x1 < x2)
  if (IsStartedTimeIntegral()) {
    Double_t l2 = (fX-oldX)*(fX-oldX) + (fY-oldY)*(fY-oldY) + (fZ-oldZ)*(fZ-oldZ);
    AddTimeStep(TMath::Sqrt(l2));
  }

  return 1;            
}     

//_____________________________________________________________________________
Int_t AliTOFtrack::PropagateToInnerTOF( Bool_t holes)
{
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0


  Double_t ymax=AliTOFGeometry::RinTOF()*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
  Bool_t skip = kFALSE;
  Double_t y=GetYat(AliTOFGeometry::RinTOF(),skip);
  if(skip){
    return 0;
  }
  if (y > ymax) {
    if (!Rotate(AliTOFGeometry::GetAlpha())) {
      return 0;
    }
  } else if (y <-ymax) {
    if (!Rotate(-AliTOFGeometry::GetAlpha())) {
      return 0;
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
  
  if(!PropagateTo(AliTOFGeometry::RinTOF()))return 0;
  
  return 1;
  
}     

//_____________________________________________________________________________
Int_t AliTOFtrack::Rotate(Double_t alpha)
{
  // Rotates track parameters in R*phi plane
  

  fAlpha += alpha;
  if (fAlpha<-TMath::Pi()) fAlpha += 2*TMath::Pi();
  if (fAlpha>=TMath::Pi()) fAlpha -= 2*TMath::Pi();

  Double_t x1=fX, y1=fY;
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=fC*fX - fE;

  fX = x1*ca + y1*sa;
  fY =-x1*sa + y1*ca;
  if((r1*r1) > 1) return 0;
  fE=fE*ca + (fC*y1 + sqrt(1.- r1*r1))*sa;

  Double_t r2=fC*fX - fE;
  if (TMath::Abs(r2) >= 0.90000) {
    AliWarning("Rotation failed !");
    return 0;
  }

  if((r2*r2) > 1) return 0;
  Double_t y0=fY + sqrt(1.- r2*r2)/fC;
  if ((fY-y0)*fC >= 0.) {
    AliWarning("Rotation failed !!!");
    return 0;
  }

  //f = F - 1
  Double_t f00=ca-1,    f24=(y1 - r1*x1/sqrt(1.- r1*r1))*sa,
           f20=fC*sa,  f22=(ca + sa*r1/sqrt(1.- r1*r1))-1;

  //b = C*ft
  Double_t b00=fCyy*f00, b02=fCyy*f20+fCcy*f24+fCey*f22;
  Double_t b10=fCzy*f00, b12=fCzy*f20+fCcz*f24+fCez*f22;
  Double_t b20=fCey*f00, b22=fCey*f20+fCce*f24+fCee*f22;
  Double_t b30=fCty*f00, b32=fCty*f20+fCct*f24+fCte*f22;
  Double_t b40=fCcy*f00, b42=fCcy*f20+fCcc*f24+fCce*f22;

  //a = f*b = f*C*ft
  Double_t a00=f00*b00, a02=f00*b02, a22=f20*b02+f24*b42+f22*b22;

  //F*C*Ft = C + (a + b + bt)
  fCyy += a00 + 2*b00;
  fCzy += b10;
  fCey += a02+b20+b02;
  fCty += b30;
  fCcy += b40;
  fCez += b12;
  fCte += b32;
  fCee += a22 + 2*b22;
  fCce += b42;

  return 1;                            
}                         

//_________________________________________________________________________
Double_t AliTOFtrack::GetYat(Double_t xk, Bool_t & skip) const {     
//-----------------------------------------------------------------
// This function calculates the Y-coordinate of a track at the plane x=xk.
// Needed for matching with the TOF (I.Belikov)
//-----------------------------------------------------------------
     skip=kFALSE;
     Double_t c1=fC*fX - fE, r1=TMath::Sqrt(TMath::Abs(1.- c1*c1));
     Double_t c2=fC*xk - fE, r2=TMath::Sqrt(TMath::Abs(1.- c2*c2));
      if( ((1.- c2*c2)<0) || ((1.- c1*c1)<0) ) skip=kTRUE;
      return fY + (xk-fX)*(c1+c2)/(r1+r2);
}
//_________________________________________________________________________
void AliTOFtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const
{
  // Returns reconstructed track momentum in the global system.

  Double_t pt=TMath::Abs(GetPt()); // GeV/c
  Double_t r=fC*fX-fE;

  Double_t y0; 
  if(r > 1) { py = pt; px = 0; }
  else if(r < -1) { py = -pt; px = 0; }
  else {
    y0=fY + sqrt(1.- r*r)/fC;  
    px=-pt*(fY-y0)*fC;    //cos(phi);
    py=-pt*(fE-fX*fC);    //sin(phi);
  }
  pz=pt*fT;
  Double_t tmp=px*TMath::Cos(fAlpha) - py*TMath::Sin(fAlpha);
  py=px*TMath::Sin(fAlpha) + py*TMath::Cos(fAlpha);
  px=tmp;            

}                                

//_________________________________________________________________________
void AliTOFtrack::GetGlobalXYZ(Double_t& x, Double_t& y, Double_t& z) const
{
  // Returns reconstructed track coordinates in the global system.

  x = fX; y = fY; z = fZ; 
  Double_t tmp=x*TMath::Cos(fAlpha) - y*TMath::Sin(fAlpha);
  y=x*TMath::Sin(fAlpha) + y*TMath::Cos(fAlpha);
  x=tmp;            

}                                

//_________________________________________________________________________
void AliTOFtrack::ResetCovariance() {
  //
  // Resets covariance matrix
  //

  fCyy*=10.;
  fCzy=0.;  fCzz*=10.;
  fCey=0.;  fCez=0.;  fCee*=10.;
  fCty=0.;  fCtz=0.;  fCte=0.;  fCtt*=10.;
  fCcy=0.;  fCcz=0.;  fCce=0.;  fCct=0.;  fCcc*=10.;  
}                                                         


//_________________________________________________________________________
void AliTOFtrack::ResetCovariance(Float_t mult) {
  //
  // Resets covariance matrix
  //

  fCyy*=mult;
  fCzy*=0.;  fCzz*=mult;
  fCey*=0.;  fCez*=0.;  fCee*=mult;
  fCty*=0.;  fCtz*=0.;  fCte*=0.;  fCtt*=mult;
  fCcy*=0.;  fCcz*=0.;  fCce*=0.;  fCct*=0.;  fCcc*=mult;  
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
  Double_t xtr,ytr,ztr;
  GetGlobalXYZ(xtr,ytr,ztr);
  Float_t thetatr = TMath::ATan2(TMath::Sqrt(xtr*xtr+ytr*ytr),ztr);

  if(holes){
    if (isec == 0 || isec == 1 || isec == 2 ) {
      if( thetatr>=thetamin && thetatr<=thetamax){ 
	x0TRD= x0Air;
	rhoTRD = rhoAir;
      }
    }
    if (isec == 11 || isec == 12 || isec == 13 || isec == 14 || isec == 15 ) {
      if( ztr>=zmin && ztr<=zmax){ 
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
