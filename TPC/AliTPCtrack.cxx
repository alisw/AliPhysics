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

//-----------------------------------------------------------------
//           Implementation of the TPC track class
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include <iostream.h>

#include "AliCluster.h"
#include "AliTPCtrack.h"
#include "AliTPCcluster.h"
#include "AliTPCClustersRow.h"
#include "AliTPCClustersArray.h"

ClassImp(AliTPCtrack)
//_________________________________________________________________________
AliTPCtrack::AliTPCtrack(UInt_t index, const Double_t xx[5],
const Double_t cc[15], Double_t xref, Double_t alpha) : AliKalmanTrack() {
  //-----------------------------------------------------------------
  // This is the main track constructor.
  //-----------------------------------------------------------------
  fdEdx=0.;

  fAlpha=alpha;
  fX=xref;

  fP0=xx[0]; fP1=xx[1]; fP2=xx[2]; fP3=xx[3]; fP4=xx[4];

  fC00=cc[0];
  fC10=cc[1];  fC11=cc[2];
  fC20=cc[3];  fC21=cc[4];  fC22=cc[5];
  fC30=cc[6];  fC31=cc[7];  fC32=cc[8];  fC33=cc[9];
  fC40=cc[10]; fC41=cc[11]; fC42=cc[12]; fC43=cc[13]; fC44=cc[14];

  fIndex[fN++]=index;
}

//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(const AliTPCtrack& t) : AliKalmanTrack(t) {
  //-----------------------------------------------------------------
  // This is a track copy constructor.
  //-----------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) fIndex[i]=t.fIndex[i];

  fdEdx=t.fdEdx;

  fAlpha=t.fAlpha;
  fX=t.fX;
}

//_____________________________________________________________________________
Int_t AliTPCtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho,Double_t pm)
{
  //-----------------------------------------------------------------
  // This function propagates a track to a reference plane x=xk.
  //-----------------------------------------------------------------
  if (TMath::Abs(fP3*xk - fP2) >= 0.99999) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Propagation failed !\n";
    return 0;
  }

  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1, y1=fP0, z1=fP1;
  Double_t c1=fP3*x1 - fP2, r1=sqrt(1.- c1*c1);
  Double_t c2=fP3*x2 - fP2, r2=sqrt(1.- c2*c2);
  
  fP0 += dx*(c1+c2)/(r1+r2);
  fP1 += dx*(c1+c2)/(c1*r2 + c2*r1)*fP4;

  //f = F - 1
  Double_t rr=r1+r2, cc=c1+c2, xx=x1+x2;
  Double_t f02=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  Double_t f03= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  Double_t cr=c1*r2+c2*r1;
  Double_t f12=-dx*fP4*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  Double_t f13=dx*fP4*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);
  Double_t f14= dx*cc/cr; 

  //b = C*ft
  Double_t b00=f02*fC20 + f03*fC30, b01=f12*fC20 + f13*fC30 + f14*fC40;
  Double_t b10=f02*fC21 + f03*fC31, b11=f12*fC21 + f13*fC31 + f14*fC41;
  Double_t b20=f02*fC22 + f03*fC32, b21=f12*fC22 + f13*fC32 + f14*fC42;
  Double_t b30=f02*fC32 + f03*fC33, b31=f12*fC32 + f13*fC33 + f14*fC43;
  Double_t b40=f02*fC42 + f03*fC43, b41=f12*fC42 + f13*fC43 + f14*fC44;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f03*b30,a01=f02*b21+f03*b31,a11=f12*b21+f13*b31+f14*b41;

  //F*C*Ft = C + (a + b + bt)
  fC00 += a00 + 2*b00;
  fC10 += a01 + b01 + b10; 
  fC20 += b20;
  fC30 += b30;
  fC40 += b40;
  fC11 += a11 + 2*b11;
  fC21 += b21; 
  fC31 += b31; 
  fC41 += b41; 

  fX=x2;

  //Multiple scattering******************
  Double_t d=sqrt((x1-fX)*(x1-fX)+(y1-fP0)*(y1-fP0)+(z1-fP1)*(z1-fP1));
  Double_t p2=GetP()*GetP();
  Double_t beta2=p2/(p2 + pm*pm);
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho/2.;
  //Double_t theta2=1.0259e-6*10*10/20/(beta2*p2)*d*rho;

  Double_t ey=fP3*fX - fP2, ez=fP4;
  Double_t xz=fP3*ez, zz1=ez*ez+1, xy=fP2+ey;

  fC33 += xz*xz*theta2;
  fC32 += xz*ez*xy*theta2;
  fC43 += xz*zz1*theta2;
  fC22 += (2*ey*ez*ez*fP2+1-ey*ey+ez*ez+fP2*fP2*ez*ez)*theta2;
  fC42 += ez*zz1*xy*theta2;
  fC44 += zz1*zz1*theta2;

  //Energy losses************************
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d*rho;
  if (x1 < x2) dE=-dE;
  cc=fP3;
  fP3*=(1.- sqrt(p2+pm*pm)/p2*dE);
  fP2+=fX*(fP3-cc);

  return 1;
}

//_____________________________________________________________________________
Int_t AliTPCtrack::PropagateToVertex(Double_t x0,Double_t rho,Double_t pm) 
{
  //-----------------------------------------------------------------
  // This function propagates tracks to the "vertex".
  //-----------------------------------------------------------------
  Double_t c=fP3*fX - fP2;
  Double_t tgf=-fP2/(fP3*fP0 + sqrt(1-c*c));
  Double_t snf=tgf/sqrt(1.+ tgf*tgf);
  Double_t xv=(fP2+snf)/fP3;
  return PropagateTo(xv,x0,rho,pm);
}

//_____________________________________________________________________________
void AliTPCtrack::Update(const AliCluster *c, Double_t chisq, UInt_t index)
{
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  Double_t cur=fP3 + k30*dy + k31*dz, eta=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.99999) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Filtering failed !\n";
    return;
  }

  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = eta;
  fP3  = cur;
  fP4 += k40*dy + k41*dz;

  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k30*c04+k31*c14; 

  fC44-=k40*c04+k41*c14; 

  fIndex[fN++]=index;
  fChi2 += chisq;
}

//_____________________________________________________________________________
Int_t AliTPCtrack::Rotate(Double_t alpha)
{
  //-----------------------------------------------------------------
  // This function rotates this track.
  //-----------------------------------------------------------------
  fAlpha += alpha;
  
  Double_t x1=fX, y1=fP0;
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=fP3*fX - fP2;
  
  fX = x1*ca + y1*sa;
  fP0=-x1*sa + y1*ca;
  fP2=fP2*ca + (fP3*y1 + sqrt(1.- r1*r1))*sa;
  
  Double_t r2=fP3*fX - fP2;
  if (TMath::Abs(r2) >= 0.99999) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Rotation failed !\n";
    return 0;
  }
  
  Double_t y0=fP0 + sqrt(1.- r2*r2)/fP3;
  if ((fP0-y0)*fP3 >= 0.) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Rotation failed !!!\n";
    return 0;
  }

  //f = F - 1
  Double_t f00=ca-1,    f23=(y1 - r1*x1/sqrt(1.- r1*r1))*sa, 
           f20=fP3*sa,  f22=(ca + sa*r1/sqrt(1.- r1*r1))-1;

  //b = C*ft
  Double_t b00=fC00*f00, b02=fC00*f20+fC30*f23+fC20*f22;
  Double_t b10=fC10*f00, b12=fC10*f20+fC31*f23+fC21*f22;
  Double_t b20=fC20*f00, b22=fC20*f20+fC32*f23+fC22*f22;
  Double_t b30=fC30*f00, b32=fC30*f20+fC33*f23+fC32*f22;
  Double_t b40=fC40*f00, b42=fC40*f20+fC43*f23+fC42*f22;

  //a = f*b = f*C*ft
  Double_t a00=f00*b00, a02=f00*b02, a22=f20*b02+f23*b32+f22*b22;

  // *** Double_t dy2=fCyy;

  //F*C*Ft = C + (a + b + bt)
  fC00 += a00 + 2*b00;
  fC10 += b10;
  fC20 += a02+b20+b02;
  fC30 += b30;
  fC40 += b40;
  fC21 += b12;
  fC32 += b32;
  fC22 += a22 + 2*b22;
  fC42 += b42; 

  // *** fCyy+=dy2*sa*sa*r1*r1/(1.- r1*r1);
  // *** fCzz+=d2y*sa*sa*fT*fT/(1.- r1*r1);

  return 1;
}

//_____________________________________________________________________________
void AliTPCtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const 
{
  //-----------------------------------------------------------------
  // This function returns reconstructed track momentum in the global system.
  //-----------------------------------------------------------------
  Double_t pt=TMath::Abs(GetPt()); // GeV/c
  Double_t r=fP3*fX-fP2;
  Double_t y0=fP0 + sqrt(1.- r*r)/fP3;
  px=-pt*(fP0-y0)*fP3;    //cos(phi);
  py=-pt*(fP2-fX *fP3);   //sin(phi);
  pz=pt*fP4;
  Double_t tmp=px*TMath::Cos(fAlpha) - py*TMath::Sin(fAlpha);
  py=px*TMath::Sin(fAlpha) + py*TMath::Cos(fAlpha);
  px=tmp;  
}

//_____________________________________________________________________________
void AliTPCtrack::CookLabel(AliTPCClustersArray *ca) {
  //-----------------------------------------------------------------
  // This function cooks the track label. If label<0, this track is fake.
  //-----------------------------------------------------------------
  Int_t *lb=new Int_t[fN];
  Int_t *mx=new Int_t[fN];
  AliTPCcluster **clusters=new AliTPCcluster*[fN];

  Int_t i;
  Int_t sec,row,ncl;
  for (i=0; i<fN; i++) {
     lb[i]=mx[i]=0;
     GetCluster(i,sec,row,ncl);
     AliTPCClustersRow *clrow=ca->GetRow(sec,row);
     clusters[i]=(AliTPCcluster*)(*clrow)[ncl];      
  }
  
  Int_t lab=123456789;
  for (i=0; i<fN; i++) {
    AliTPCcluster *c=clusters[i];
    lab=TMath::Abs(c->GetLabel(0));
    Int_t j;
    for (j=0; j<fN; j++)
      if (lb[j]==lab || mx[j]==0) break;
    lb[j]=lab;
    (mx[j])++;
  }
  
  Int_t max=0;
  for (i=0; i<fN; i++) 
    if (mx[i]>max) {max=mx[i]; lab=lb[i];}
    
  for (i=0; i<fN; i++) {
    AliTPCcluster *c=clusters[i];
    if (TMath::Abs(c->GetLabel(1)) == lab ||
        TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }
  
  SetLabel(-lab);
  if (1.-Float_t(max)/fN <= 0.10) {
    //Int_t tail=Int_t(0.08*fN);
     Int_t tail=14;
     max=0;
     for (i=1; i<=tail; i++) {
       AliTPCcluster *c=clusters[fN-i];
       if (lab == TMath::Abs(c->GetLabel(0)) ||
           lab == TMath::Abs(c->GetLabel(1)) ||
           lab == TMath::Abs(c->GetLabel(2))) max++;
     }
     if (max >= Int_t(0.5*tail)) SetLabel(lab);
  }

  delete[] lb;
  delete[] mx;
  delete[] clusters;
}
