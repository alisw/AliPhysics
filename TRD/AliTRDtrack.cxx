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

/*
$Log$
Revision 1.13  2002/10/22 15:53:08  alibrary
Introducing Riostream.h

Revision 1.12  2002/10/14 14:57:44  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.8.10.2  2002/07/24 10:09:31  alibrary
Updating VirtualMC

RRevision 1.11  2002/06/13 12:09:58  hristov
Minor corrections

Revision 1.10  2002/06/12 09:54:35  cblume
Update of tracking code provided by Sergei

Revision 1.8  2001/05/30 12:17:47  hristov
Loop variables declared once

Revision 1.7  2001/05/28 17:07:58  hristov
Last minute changes; ExB correction in AliTRDclusterizerV1; taking into account of material in G10 TEC frames and material between TEC planes (C.Blume,S.Sedykh)

Revision 1.4  2000/12/08 16:07:02  cblume
Update of the tracking by Sergei

Revision 1.3  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.2  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.2.1  2000/09/22 14:47:52  cblume
Add the tracking code

*/                                                        

#include <Riostream.h>
#include <TObject.h>   

#include "AliTRDgeometry.h" 
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "../TPC/AliTPCtrack.h" 

ClassImp(AliTRDtrack)


//_____________________________________________________________________________

AliTRDtrack::AliTRDtrack(const AliTRDcluster *c, UInt_t index, 
                         const Double_t xx[5], const Double_t cc[15], 
                         Double_t xref, Double_t alpha) : AliKalmanTrack() {
  //-----------------------------------------------------------------
  // This is the main track constructor.
  //-----------------------------------------------------------------

  fSeedLab = -1;

  fAlpha=alpha;
  if (fAlpha<-TMath::Pi()) fAlpha += 2*TMath::Pi();
  if (fAlpha>=TMath::Pi()) fAlpha -= 2*TMath::Pi();   

  fX=xref;

  fY=xx[0]; fZ=xx[1]; fC=xx[2]; fE=xx[3]; fT=xx[4];

  fCyy=cc[0];
  fCzy=cc[1];  fCzz=cc[2];
  fCcy=cc[3];  fCcz=cc[4];  fCcc=cc[5];
  fCey=cc[6];  fCez=cc[7];  fCec=cc[8];  fCee=cc[9];
  fCty=cc[10]; fCtz=cc[11]; fCtc=cc[12]; fCte=cc[13]; fCtt=cc[14];

  fIndex[0]=index;
  SetNumberOfClusters(1);

  fdEdx=0.;

  Double_t q = TMath::Abs(c->GetQ());
  Double_t s = fX*fC - fE, t=fT;
  if(s*s < 1) q *= TMath::Sqrt((1-s*s)/(1+t*t));

  fdQdl[0] = q;
}                              
           
//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliTRDtrack& t) : AliKalmanTrack(t) {
  //
  // Copy constructor.
  //

  SetLabel(t.GetLabel());
  fSeedLab=t.GetSeedLabel();

  SetChi2(t.GetChi2());
  fdEdx=t.fdEdx;

  fAlpha=t.fAlpha;
  fX=t.fX;

  fY=t.fY; fZ=t.fZ; fC=t.fC; fE=t.fE; fT=t.fT;

  fCyy=t.fCyy;
  fCzy=t.fCzy;  fCzz=t.fCzz;
  fCcy=t.fCcy;  fCcz=t.fCcz;  fCcc=t.fCcc;
  fCey=t.fCey;  fCez=t.fCez;  fCec=t.fCec;  fCee=t.fCee;
  fCty=t.fCty;  fCtz=t.fCtz;  fCtc=t.fCtc;  fCte=t.fCte;  fCtt=t.fCtt;

  Int_t n=t.GetNumberOfClusters(); 
  SetNumberOfClusters(n);
  for (Int_t i=0; i<n; i++) {
    fIndex[i]=t.fIndex[i];
    fdQdl[i]=t.fdQdl[i];
  }
}                                

//_____________________________________________________________________________
AliTRDtrack::AliTRDtrack(const AliKalmanTrack& t, Double_t alpha) {
  //
  // Constructor from AliTPCtrack or AliITStrack .
  //

  SetLabel(t.GetLabel());
  SetChi2(0.);
  SetNumberOfClusters(0);

  fdEdx=0;

  fAlpha = alpha;
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();

  Double_t x, p[5]; t.GetExternalParameters(x,p);

  fX=x;

  x = GetConvConst();  

  fY=p[0]; fZ=p[1]; fC=p[4]/x;
  fE=fX*fC-p[2]; fT=p[3];

  //Conversion of the covariance matrix
  Double_t c[15]; t.GetExternalCovariance(c);

  c[10]/=x; c[11]/=x; c[12]/=x; c[13]/=x; c[14]/=x*x;

  Double_t c22=fX*fX*c[14] - 2*fX*c[12] + c[5];
  Double_t c32=fX*c[13] - c[8];
  Double_t c20=fX*c[10] - c[3], c21=fX*c[11] - c[4], c42=fX*c[14] - c[12];

  fCyy=c[0 ];
  fCzy=c[1 ];   fCzz=c[2 ];
  fCcy=c[10];   fCcz=c[11];   fCcc=c[14];
  fCey=c20;     fCez=c21;     fCec=c42;     fCee=c22;
  fCty=c[6 ];   fCtz=c[7 ];   fCtc=c[13];   fCte=c32;   fCtt=c[9 ];

}              

//____________________________________________________________________________
void AliTRDtrack::GetExternalParameters(Double_t& xr, Double_t x[5]) const {
  //
  // This function returns external TRD track representation
  //
     xr=fX;
     x[0]=GetY();
     x[1]=GetZ();
     x[2]=GetSnp();
     x[3]=GetTgl();
     x[4]=fC*GetConvConst();
}           

//_____________________________________________________________________________
void AliTRDtrack::GetExternalCovariance(Double_t cc[15]) const {
  //
  // This function returns external representation of the covriance matrix.
  //
  Double_t a=GetConvConst();

  Double_t c22=fX*fX*fCcc-2*fX*fCec+fCee;
  Double_t c32=fX*fCtc-fCte;
  Double_t c20=fX*fCcy-fCey, c21=fX*fCcz-fCez, c42=fX*fCcc-fCec;

  cc[0 ]=fCyy;
  cc[1 ]=fCzy;   cc[2 ]=fCzz;
  cc[3 ]=c20;    cc[4 ]=c21;    cc[5 ]=c22;
  cc[6 ]=fCty;   cc[7 ]=fCtz;   cc[8 ]=c32;   cc[9 ]=fCtt;
  cc[10]=fCcy*a; cc[11]=fCcz*a; cc[12]=c42*a; cc[13]=fCtc*a; cc[14]=fCcc*a*a;
}               
                       

//_____________________________________________________________________________
void AliTRDtrack::GetCovariance(Double_t cc[15]) const {
  cc[0]=fCyy;
  cc[1]=fCzy;  cc[2]=fCzz;
  cc[3]=fCcy;  cc[4]=fCcz;  cc[5]=fCcc;
  cc[6]=fCey;  cc[7]=fCez;  cc[8]=fCec;  cc[9]=fCee;
  cc[10]=fCty; cc[11]=fCtz; cc[12]=fCtc; cc[13]=fCte; cc[14]=fCtt;
}    

//_____________________________________________________________________________
Int_t AliTRDtrack::Compare(const TObject *o) const {

// Compares tracks according to their Y2 or curvature

  AliTRDtrack *t=(AliTRDtrack*)o;
  //  Double_t co=t->GetSigmaY2();
  //  Double_t c =GetSigmaY2();

  Double_t co=TMath::Abs(t->GetC());
  Double_t c =TMath::Abs(GetC());  

  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}                

//_____________________________________________________________________________
void AliTRDtrack::CookdEdx(Double_t low, Double_t up) {
  //-----------------------------------------------------------------
  // Calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Int_t i;
  Int_t nc=GetNumberOfClusters(); 

  Float_t sorted[kMAX_CLUSTERS_PER_TRACK];
  for (i=0; i < nc; i++) {
    sorted[i]=fdQdl[i];
  }

  Int_t swap; 

  do {
    swap=0;
    for (i=0; i<nc-1; i++) {
      if (sorted[i]<=sorted[i+1]) continue;
      Float_t tmp=sorted[i];
      sorted[i]=sorted[i+1]; sorted[i+1]=tmp;
      swap++;
    }
  } while (swap);

  Int_t nl=Int_t(low*nc), nu=Int_t(up*nc);
  Float_t dedx=0;
  for (i=nl; i<=nu; i++) dedx += sorted[i];
  dedx /= (nu-nl+1);

  SetdEdx(dedx);
}                     




//_____________________________________________________________________________
Int_t AliTRDtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho,Double_t pm)
{
  // Propagates a track of particle with mass=pm to a reference plane 
  // defined by x=xk through media of density=rho and radiationLength=x0

  
  if (TMath::Abs(fC*xk - fE) >= 0.99999) {
    Int_t n=GetNumberOfClusters(); 
    if (n>4) cerr<<n<<" AliTRDtrack warning: Propagation failed !\n";
    return 0;
  }

  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1, y1=fY, z1=fZ;
  Double_t c1=fC*x1 - fE;
  if((c1*c1) > 1) return 0;
  Double_t r1=sqrt(1.- c1*c1);
  Double_t c2=fC*x2 - fE;
  if((c2*c2) > 1) return 0;
  Double_t r2=sqrt(1.- c2*c2);

  fY += dx*(c1+c2)/(r1+r2);
  fZ += dx*(c1+c2)/(c1*r2 + c2*r1)*fT;

  //f = F - 1
  Double_t rr=r1+r2, cc=c1+c2, xx=x1+x2;
  Double_t f02= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  Double_t f03=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  Double_t cr=c1*r2+c2*r1;
  Double_t f12= dx*fT*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);
  Double_t f13=-dx*fT*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  Double_t f14= dx*cc/cr;

  //b = C*ft
  Double_t b00=f02*fCcy + f03*fCey, b01=f12*fCcy + f13*fCey + f14*fCty;
  Double_t b10=f02*fCcz + f03*fCez, b11=f12*fCcz + f13*fCez + f14*fCtz;
  Double_t b20=f02*fCcc + f03*fCec, b21=f12*fCcc + f13*fCec + f14*fCtc;
  Double_t b30=f02*fCec + f03*fCee, b31=f12*fCec + f13*fCee + f14*fCte;
  Double_t b40=f02*fCtc + f03*fCte, b41=f12*fCtc + f13*fCte + f14*fCtt;

  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f03*b30,a01=f02*b21+f03*b31,a11=f12*b21+f13*b31+f14*b41;

  //F*C*Ft = C + (a + b + bt)
  fCyy += a00 + 2*b00;
  fCzy += a01 + b01 + b10;
  fCcy += b20;
  fCey += b30;
  fCty += b40;
  fCzz += a11 + 2*b11;
  fCcz += b21;
  fCez += b31;
  fCtz += b41;                  

  fX=x2;


  //Multiple scattering  ******************

  Double_t d=sqrt((x1-fX)*(x1-fX)+(y1-fY)*(y1-fY)+(z1-fZ)*(z1-fZ));
  Double_t p2=GetPt()*GetPt()*(1.+fT*fT);
  p2 = TMath::Min(p2,1e+08);  // to avoid division by (1-1) for stiff tracks
  Double_t beta2=p2/(p2 + pm*pm);

  Double_t ey=fC*fX - fE, ez=fT;
  Double_t xz=fC*ez, zz1=ez*ez+1, xy=fE+ey;

  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho;
  fCcc += xz*xz*theta2;
  fCec += xz*ez*xy*theta2;
  fCtc += xz*zz1*theta2;
  fCee += (2*ey*ez*ez*fE+1-ey*ey+ez*ez+fE*fE*ez*ez)*theta2;
  fCte += ez*zz1*xy*theta2;
  fCtt += zz1*zz1*theta2;


  //Energy losses************************
  if (x1 < x2) d=-d;
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d*rho;
  //PH  SetLength(GetLength()+d);

  cc = fC;
  fC*=(1.- sqrt(p2+pm*pm)/p2*dE);
  fE+=fX*(fC-cc);

  return 1;        

}     


//_____________________________________________________________________________
Int_t AliTRDtrack::Update(const AliTRDcluster *c, Double_t chisq, UInt_t index,                          Double_t h01)
{
  // Assignes found cluster to the track and updates track information


  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2()*100.0;

  r00+=(fCyy+2.0*h01*fCzy+h01*h01*fCzz);
  r01+=(fCzy+h01*fCzz);
  r11+=fCzz;

  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

//  Double_t k00=fCyy*r00+fCzy*r01, k01=fCyy*r01+fCzy*r11;
//  Double_t k10=fCzy*r00+fCzz*r01, k11=fCzy*r01+fCzz*r11;
//  Double_t k20=fCcy*r00+fCcz*r01, k21=fCcy*r01+fCcz*r11;
//  Double_t k30=fCey*r00+fCez*r01, k31=fCey*r01+fCez*r11;
//  Double_t k40=fCty*r00+fCtz*r01, k41=fCty*r01+fCtz*r11;

  Double_t k00=fCyy*r00+fCzy*(r01+h01*r00),k01=fCyy*r01+fCzy*(r11+h01*r01);
  Double_t k10=fCzy*r00+fCzz*(r01+h01*r00),k11=fCzy*r01+fCzz*(r11+h01*r01);
  Double_t k20=fCcy*r00+fCcz*(r01+h01*r00),k21=fCcy*r01+fCcz*(r11+h01*r01);
  Double_t k30=fCey*r00+fCez*(r01+h01*r00),k31=fCey*r01+fCez*(r11+h01*r01);
  Double_t k40=fCty*r00+fCtz*(r01+h01*r00),k41=fCty*r01+fCtz*(r11+h01*r01);

  Double_t dy=c->GetY() - fY, dz=c->GetZ() - fZ;

  dy=dy+h01*dz;

  Double_t cur=fC + k20*dy + k21*dz, eta=fE + k30*dy + k31*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.99999) {
  //  if (fN>4) cerr<<fN<<" AliTRDtrack warning: Filtering failed !\n";
    return 0;
  }

  fY += k00*dy + k01*dz;
  fZ += k10*dy + k11*dz;
  fC  = cur;
  fE  = eta;
  fT += k40*dy + k41*dz;


  k01+=h01*k00;
  k11+=h01*k10;
  k21+=h01*k20;
  k31+=h01*k30;
  k41+=h01*k40;

  Double_t c01=fCzy, c02=fCcy, c03=fCey, c04=fCty;
  Double_t c12=fCcz, c13=fCez, c14=fCtz;

  fCyy-=k00*fCyy+k01*fCzy; fCzy-=k00*c01+k01*fCzz;
  fCcy-=k00*c02+k01*c12; fCey-=k00*c03+k01*c13;
  fCty-=k00*c04+k01*c14;

  fCzz-=k10*c01+k11*fCzz;
  fCcz-=k10*c02+k11*c12; fCez-=k10*c03+k11*c13;
  fCtz-=k10*c04+k11*c14;

  fCcc-=k20*c02+k21*c12; fCec-=k20*c03+k21*c13;
  fCtc-=k20*c04+k21*c14;

  fCee-=k30*c03+k31*c13;
  fCte-=k30*c04+k31*c14;        

  fCtt-=k40*c04+k41*c14;

  Int_t n=GetNumberOfClusters();  
  fIndex[n]=index;
  SetNumberOfClusters(n+1);  

  SetChi2(GetChi2()+chisq); 

  //  cerr<<"in update: fIndex["<<fN<<"] = "<<index<<endl;

  return 1;

}                     


//_____________________________________________________________________________
Int_t AliTRDtrack::Rotate(Double_t alpha)
{
  // Rotates track parameters in R*phi plane

  fAlpha += alpha;

  Double_t x1=fX, y1=fY;
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=fC*fX - fE;

  if (TMath::Abs(r1) >= 0.99999) {
    Int_t n=GetNumberOfClusters(); 
    if (n>4) cerr<<n<<" AliTRDtrack warning: Rotation failed !\n";
    return 0;
  }

  fX = x1*ca + y1*sa;
  fY=-x1*sa + y1*ca;
  fE=fE*ca + (fC*y1 + sqrt(1.- r1*r1))*sa;

  Double_t r2=fC*fX - fE;
  if (TMath::Abs(r2) >= 0.99999) {
    Int_t n=GetNumberOfClusters(); 
    if (n>4) cerr<<n<<" AliTRDtrack warning: Rotation failed !\n";
    return 0;
  }

  Double_t y0=fY + sqrt(1.- r2*r2)/fC;
  if ((fY-y0)*fC >= 0.) {
    Int_t n=GetNumberOfClusters(); 
    if (n>4) cerr<<n<<" AliTRDtrack warning: Rotation failed !!!\n";
    return 0;
  }

  //f = F - 1
  Double_t f00=ca-1,    f32=(y1 - r1*x1/sqrt(1.- r1*r1))*sa,
           f30=fC*sa, f33=(ca + sa*r1/sqrt(1.- r1*r1))-1;

  //b = C*ft
  Double_t b00=fCyy*f00, b03=fCyy*f30+fCcy*f32+fCey*f33;
  Double_t b10=fCzy*f00, b13=fCzy*f30+fCcz*f32+fCez*f33;
  Double_t b20=fCcy*f00, b23=fCcy*f30+fCcc*f32+fCec*f33;
  Double_t b30=fCey*f00, b33=fCey*f30+fCec*f32+fCee*f33;
  Double_t b40=fCty*f00, b43=fCty*f30+fCtc*f32+fCte*f33;

  //a = f*b = f*C*ft
  Double_t a00=f00*b00, a03=f00*b03, a33=f30*b03+f32*b23+f33*b33;

  // *** Double_t dy2=fCyy;  
          
  //F*C*Ft = C + (a + b + bt)
  fCyy += a00 + 2*b00;
  fCzy += b10;
  fCcy += b20;
  fCey += a03+b30+b03;
  fCty += b40;
  fCez += b13;
  fCec += b23;
  fCee += a33 + 2*b33;
  fCte += b43;

  // *** fCyy+=dy2*sa*sa*r1*r1/(1.- r1*r1);
  // *** fCzz+=d2y*sa*sa*fT*fT/(1.- r1*r1);   

  return 1;
}                         


//_____________________________________________________________________________
Double_t AliTRDtrack::GetPredictedChi2(const AliTRDcluster *c, Double_t h01) const
{
  
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();

  r00+=fCyy; r01+=fCzy; r11+=fCzz;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    Int_t n=GetNumberOfClusters(); 
    if (n>4) cerr<<n<<" AliTRDtrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;

  Double_t dy=c->GetY() - fY, dz=c->GetZ() - fZ;

  dy=dy+h01*dz;

  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det; 

}      


//_________________________________________________________________________
void AliTRDtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const
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
    py=-pt*(fE-fX*fC);   //sin(phi);
  }
  pz=pt*fT;
  Double_t tmp=px*TMath::Cos(fAlpha) - py*TMath::Sin(fAlpha);
  py=px*TMath::Sin(fAlpha) + py*TMath::Cos(fAlpha);
  px=tmp;            

}                                

//_________________________________________________________________________
void AliTRDtrack::GetGlobalXYZ(Double_t& x, Double_t& y, Double_t& z) const
{
  // Returns reconstructed track coordinates in the global system.

  x = fX; y = fY; z = fZ; 
  Double_t tmp=x*TMath::Cos(fAlpha) - y*TMath::Sin(fAlpha);
  y=x*TMath::Sin(fAlpha) + y*TMath::Cos(fAlpha);
  x=tmp;            

}                                

//_________________________________________________________________________
void AliTRDtrack::ResetCovariance() {
  //
  // Resets covariance matrix
  //

  fCyy*=10.;
  fCzy=0.;   fCzz*=10.;
  fCcy=0.;   fCcz=0.;   fCcc*=10.;
  fCey=0.;   fCez=0.;   fCec=0.;   fCee*=10.;
  fCty=0.;   fCtz=0.;   fCtc=0.;   fCte=0.;   fCtt*=10.;

}                                                         
