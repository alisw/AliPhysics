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
Revision 1.2  2000/06/30 12:07:50  kowal2
Updated from the TPC-PreRelease branch

Revision 1.1.2.2  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking

*/

//-----------------------------------------------------------------
//           Implementation of the TPC track class
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliTPCtrack.h"
#include "AliTPCcluster.h"
#include "AliTPCClustersRow.h"
#include "AliTPCClustersArray.h"

ClassImp(AliTPCtrack)
//_________________________________________________________________________
AliTPCtrack::AliTPCtrack(UInt_t index, const Double_t xx[5],
const Double_t cc[15], Double_t xref, Double_t alpha) {
  //-----------------------------------------------------------------
  // This is the main track constructor.
  //-----------------------------------------------------------------
  fLab=-1;
  fChi2=0.;
  fdEdx=0.;

  fAlpha=alpha;
  fX=xref;

  fY=xx[0]; fZ=xx[1]; fC=xx[2]; fE=xx[3]; fT=xx[4];

  fCyy=cc[0];
  fCzy=cc[1];  fCzz=cc[2];
  fCcy=cc[3];  fCcz=cc[4];  fCcc=cc[5];
  fCey=cc[6];  fCez=cc[7];  fCec=cc[8];  fCee=cc[9];
  fCty=cc[10]; fCtz=cc[11]; fCtc=cc[12]; fCte=cc[13]; fCtt=cc[14];

  fN=0;
  fIndex[fN++]=index;
}

//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(const AliTPCtrack& t) {
  //-----------------------------------------------------------------
  // This is a track copy constructor.
  //-----------------------------------------------------------------
  fLab=t.fLab;
  fChi2=t.fChi2;
  fdEdx=t.fdEdx;

  fAlpha=t.fAlpha;
  fX=t.fX;

  fY=t.fY; fZ=t.fZ; fC=t.fC; fE=t.fE; fT=t.fT;

  fCyy=t.fCyy;
  fCzy=t.fCzy;  fCzz=t.fCzz;
  fCcy=t.fCcy;  fCcz=t.fCcz;  fCcc=t.fCcc;
  fCey=t.fCey;  fCez=t.fCez;  fCec=t.fCec;  fCee=t.fCee;
  fCty=t.fCty;  fCtz=t.fCtz;  fCtc=t.fCtc;  fCte=t.fCte;  fCtt=t.fCtt;

  fN=t.fN;
  for (Int_t i=0; i<fN; i++) fIndex[i]=t.fIndex[i];
}

//_____________________________________________________________________________
void AliTPCtrack::GetCovariance(Double_t cc[15]) const {
  //just to calm down our rule checker
  cc[0]=fCyy;
  cc[1]=fCzy;  cc[2]=fCzz;
  cc[3]=fCcy;  cc[4]=fCcz;  cc[5]=fCcc;
  cc[6]=fCey;  cc[7]=fCez;  cc[8]=fCec;  cc[9]=fCee;
  cc[10]=fCty; cc[11]=fCtz; cc[12]=fCtc; cc[13]=fCte; cc[14]=fCtt;
}

//_____________________________________________________________________________
Int_t AliTPCtrack::Compare(TObject *o) {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliTPCtrack *t=(AliTPCtrack*)o;
  //Double_t co=t->GetSigmaY2();
  //Double_t c =GetSigmaY2();
  Double_t co=TMath::Abs(t->GetC());
  Double_t c =TMath::Abs(GetC());
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}

//_____________________________________________________________________________
Int_t AliTPCtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho,Double_t pm)
{
  //-----------------------------------------------------------------
  // This function propagates a track to a reference plane x=xk.
  //-----------------------------------------------------------------
  if (TMath::Abs(fC*xk - fE) >= 0.99999) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Propagation failed !\n";
    return 0;
  }

  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1, y1=fY, z1=fZ;
  Double_t c1=fC*x1 - fE, r1=sqrt(1.- c1*c1);
  Double_t c2=fC*x2 - fE, r2=sqrt(1.- c2*c2);
  
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

  //Multiple scattering******************
  Double_t d=sqrt((x1-fX)*(x1-fX)+(y1-fY)*(y1-fY)+(z1-fZ)*(z1-fZ));
  Double_t p2=GetPt()*GetPt()*(1.+fT*fT);
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
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d*rho;
  if (x1 < x2) dE=-dE;
  cc=fC;
  fC*=(1.- sqrt(p2+pm*pm)/p2*dE);
  fE+=fX*(fC-cc);

  return 1;
}

//_____________________________________________________________________________
void AliTPCtrack::PropagateToVertex(Double_t x0,Double_t rho,Double_t pm) 
{
  //-----------------------------------------------------------------
  // This function propagates tracks to the "vertex".
  //-----------------------------------------------------------------
  Double_t c=fC*fX - fE;
  Double_t tgf=-fE/(fC*fY + sqrt(1-c*c));
  Double_t snf=tgf/sqrt(1.+ tgf*tgf);
  Double_t xv=(fE+snf)/fC;
  PropagateTo(xv,x0,rho,pm);
}

//_____________________________________________________________________________
void AliTPCtrack::Update(const AliTPCcluster *c, Double_t chisq, UInt_t index)
{
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  r00+=fCyy; r01+=fCzy; r11+=fCzz;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

  Double_t k00=fCyy*r00+fCzy*r01, k01=fCyy*r01+fCzy*r11;
  Double_t k10=fCzy*r00+fCzz*r01, k11=fCzy*r01+fCzz*r11;
  Double_t k20=fCcy*r00+fCcz*r01, k21=fCcy*r01+fCcz*r11;
  Double_t k30=fCey*r00+fCez*r01, k31=fCey*r01+fCez*r11;
  Double_t k40=fCty*r00+fCtz*r01, k41=fCty*r01+fCtz*r11;

  Double_t dy=c->GetY() - fY, dz=c->GetZ() - fZ;
  Double_t cur=fC + k20*dy + k21*dz, eta=fE + k30*dy + k31*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.99999) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Filtering failed !\n";
    return;
  }

  fY += k00*dy + k01*dz;
  fZ += k10*dy + k11*dz;
  fC  = cur;
  fE  = eta;
  fT += k40*dy + k41*dz;

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
  
  Double_t x1=fX, y1=fY;
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=fC*fX - fE;
  
  fX = x1*ca + y1*sa;
  fY=-x1*sa + y1*ca;
  fE=fE*ca + (fC*y1 + sqrt(1.- r1*r1))*sa;
  
  Double_t r2=fC*fX - fE;
  if (TMath::Abs(r2) >= 0.99999) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Rotation failed !\n";
    return 0;
  }
  
  Double_t y0=fY + sqrt(1.- r2*r2)/fC;
  if ((fY-y0)*fC >= 0.) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Rotation failed !!!\n";
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
Double_t AliTPCtrack::GetPredictedChi2(const AliTPCcluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  r00+=fCyy; r01+=fCzy; r11+=fCzz;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    if (fN>4) cerr<<fN<<" AliTPCtrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=c->GetY() - fY, dz=c->GetZ() - fZ;
  
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}

//_____________________________________________________________________________
void AliTPCtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const 
{
  //-----------------------------------------------------------------
  // This function returns reconstructed track momentum in the global system.
  //-----------------------------------------------------------------
  Double_t pt=TMath::Abs(GetPt()); // GeV/c
  Double_t r=fC*fX-fE;
  Double_t y0=fY + sqrt(1.- r*r)/fC;
  px=-pt*(fY-y0)*fC;    //cos(phi);
  py=-pt*(fE-fX*fC);   //sin(phi);
  pz=pt*fT;
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

//____________________________________________________________________________
void AliTPCtrack::Streamer(TBuffer &R__b)
{
  //-----------------------------------------------------
  // This is AliTPCtrack streamer.
  //-----------------------------------------------------
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fLab;
      R__b >> fChi2;
      R__b >> fdEdx;
      R__b >> fAlpha;
      R__b >> fX;
      R__b >> fY;
      R__b >> fZ;
      R__b >> fC;
      R__b >> fE;
      R__b >> fT;
      R__b >> fCyy;
      R__b >> fCzy;
      R__b >> fCzz;
      R__b >> fCcy;
      R__b >> fCcz;
      R__b >> fCcc;
      R__b >> fCey;
      R__b >> fCez;
      R__b >> fCec;
      R__b >> fCee;
      R__b >> fCty;
      R__b >> fCtz;
      R__b >> fCtc;
      R__b >> fCte;
      R__b >> fCtt;
      R__b >> fN;
      for (Int_t i=0; i<fN; i++) R__b >> fIndex[i];
   } else {
      R__b.WriteVersion(AliTPCtrack::IsA());
      TObject::Streamer(R__b);
      R__b << fLab;
      R__b << fChi2;
      R__b << fdEdx;
      R__b << fAlpha;
      R__b << fX;
      R__b << fY;
      R__b << fZ;
      R__b << fC;
      R__b << fE;
      R__b << fT;
      R__b << fCyy;
      R__b << fCzy;
      R__b << fCzz;
      R__b << fCcy;
      R__b << fCcz;
      R__b << fCcc;
      R__b << fCey;
      R__b << fCez;
      R__b << fCec;
      R__b << fCee;
      R__b << fCty;
      R__b << fCtz;
      R__b << fCtc;
      R__b << fCte;
      R__b << fCtt;
      R__b << fN;
      for (Int_t i=0; i<fN; i++) R__b << fIndex[i];
   }
}


