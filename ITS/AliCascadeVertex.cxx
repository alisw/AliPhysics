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
//               Implementation of the cascade vertex class
//
//    Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr
//-------------------------------------------------------------------------
#include <Riostream.h>
#include <TMath.h>
#include <TPDGCode.h>

#include "AliCascadeVertex.h"
#include "AliITStrackV2.h"
#include "AliV0vertex.h"

ClassImp(AliCascadeVertex)

AliCascadeVertex::AliCascadeVertex() : TObject() {
  //--------------------------------------------------------------------
  // Default constructor  (Xi-)
  //--------------------------------------------------------------------
  fPdgCode=kXiMinus;
  fEffMass=1.32131;
  fChi2=1.e+33;
  fPos[0]=fPos[1]=fPos[2]=0.;
  fPosCov[0]=fPosCov[1]=fPosCov[2]=fPosCov[3]=fPosCov[4]=fPosCov[5]=0.;
}



inline Double_t det(Double_t a00, Double_t a01, Double_t a10, Double_t a11){
  // determinant 2x2
  return a00*a11 - a01*a10;
}

inline Double_t det (Double_t a00,Double_t a01,Double_t a02,
                     Double_t a10,Double_t a11,Double_t a12,
                     Double_t a20,Double_t a21,Double_t a22) {
  // determinant 3x3
  return 
  a00*det(a11,a12,a21,a22)-a01*det(a10,a12,a20,a22)+a02*det(a10,a11,a20,a21);
}



AliCascadeVertex::AliCascadeVertex(const AliV0vertex &v,const AliITStrackV2 &t) {
  //--------------------------------------------------------------------
  // Main constructor
  //--------------------------------------------------------------------
  fPdgCode=kXiMinus;

  fV0lab[0]=v.GetNlabel(); fV0lab[1]=v.GetPlabel();
  fBachLab=t.GetLabel(); 

  //Trivial estimation of the vertex parameters
  Double_t pt, phi, x, par[5];
  Double_t alpha, cs, sn;

  t.GetExternalParameters(x,par); alpha=t.GetAlpha();
  pt=1./TMath::Abs(par[4]);
  phi=TMath::ASin(par[2]) + alpha;  

  // momentum of the bachelor track

  Double_t px1=pt*TMath::Cos(phi), py1=pt*TMath::Sin(phi), pz1=pt*par[3];

  cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);

  Double_t x1=x*cs - par[0]*sn; // position of the bachelor at dca (bachelor,V0)
  Double_t y1=x*sn + par[0]*cs;
  Double_t z1=par[1];

  Double_t x2,y2,z2;          // position of the V0 
  v.GetXYZ(x2,y2,z2);
    
  Double_t px2,py2,pz2;       // momentum of V0
  v.GetPxPyPz(px2,py2,pz2);

  Double_t a2=((x1-x2)*px2+(y1-y2)*py2+(z1-z2)*pz2)/(px2*px2+py2*py2+pz2*pz2);

  Double_t xm=x2+a2*px2;
  Double_t ym=y2+a2*py2;
  Double_t zm=z2+a2*pz2;

  // position of the cascade decay
  
  fPos[0]=0.5*(x1+xm); fPos[1]=0.5*(y1+ym); fPos[2]=0.5*(z1+zm);
    

  // invariant mass of the cascade (default is Ximinus)
  
  Double_t e1=TMath::Sqrt(0.13957*0.13957 + px1*px1 + py1*py1 + pz1*pz1);
  Double_t e2=TMath::Sqrt(1.11568*1.11568 + px2*px2 + py2*py2 + pz2*pz2);
  
  fEffMass=TMath::Sqrt((e1+e2)*(e1+e2)-
    (px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));


  // momenta of the bachelor and the V0
  
  fBachMom[0]=px1; fBachMom[1]=py1; fBachMom[2]=pz1; 
  v.GetNPxPyPz(px2,py2,pz2);
  fV0mom[0][0]=px2; fV0mom[0][1]=py2; fV0mom[0][2]=pz2;
  v.GetPPxPyPz(px2,py2,pz2);
  fV0mom[1][0]=px2; fV0mom[1][1]=py2; fV0mom[1][2]=pz2;


  fChi2=7.;   

}

/*
Double_t AliCascadeVertex::ChangeMassHypothesis(Double_t &v0q, Int_t code) {
  //--------------------------------------------------------------------
  // This function changes the mass hypothesis for this cascade
  // and returns the "kinematical quality" of this hypothesis
  // together with the "quality" of associated V0 (argument v0q) 
  //--------------------------------------------------------------------
  Double_t nmass=0.13957, pmass=0.93827, des0=0.9437-0.1723; 
  Double_t bmass=0.13957, mass =1.3213,  des =1.1243-0.1970;

  fPdgCode=code;

  switch (code) {
  case 213: 
       bmass=0.93827; 
       break;
  case kXiMinus:
       break;
  case kXiPlusBar:
       nmass=0.93827; pmass=0.13957; des0=-des0; 
       des=-des;
       break;
  case kOmegaMinus: 
       bmass=0.49368; mass=1.67245; des=1.1355-0.5369;
       break;
  case kOmegaPlusBar: 
       nmass=0.93827; pmass=0.13957; des0=-des0; 
       bmass=0.49368; mass=1.67245; des=0.5369-1.1355;
       break;
  default:
       cerr<<"AliCascadeVertex::ChangeMassHypothesis: ";
       cerr<<"Invalide PDG code !  Assuming XiMinus's...\n";
       fPdgCode=kXiMinus;
    break;
  }

  Double_t pxn=fV0mom[0][0], pyn=fV0mom[0][1], pzn=fV0mom[0][2];
  Double_t pxp=fV0mom[1][0], pyp=fV0mom[1][1], pzp=fV0mom[1][2];
  Double_t en=TMath::Sqrt(nmass*nmass + pxn*pxn + pyn*pyn + pzn*pzn);
  Double_t ep=TMath::Sqrt(pmass*pmass + pxp*pxp + pyp*pyp + pzp*pzp);
  Double_t px0=pxn+pxp, py0=pyn+pyp, pz0=pzn+pzp;
  Double_t p0=TMath::Sqrt(px0*px0 + py0*py0 + pz0*pz0);

  Double_t gamma0=(en+ep)/1.11568, betagamma0=p0/1.11568;
  Double_t pln=(pxn*px0 + pyn*py0 + pzn*pz0)/p0;
  Double_t plp=(pxp*px0 + pyp*py0 + pzp*pz0)/p0;
  Double_t plps=gamma0*plp - betagamma0*ep;

  Double_t diff0=2*gamma0*plps + betagamma0*des0;


  v0q=plp-pln-diff0;


  Double_t pxb=fBachMom[0], pyb=fBachMom[1], pzb=fBachMom[2]; 

  Double_t e0=TMath::Sqrt(1.11568*1.11568 + p0*p0);
  Double_t eb=TMath::Sqrt(bmass*bmass + pxb*pxb + pyb*pyb + pzb*pzb);
  Double_t pxl=px0+pxb, pyl=py0+pyb, pzl=pz0+pzb;
  Double_t pl=TMath::Sqrt(pxl*pxl + pyl*pyl + pzl*pzl);
  
  fEffMass=TMath::Sqrt((e0+eb)*(e0+eb) - pl*pl);

  Double_t gamma=(e0+eb)/mass, betagamma=pl/mass;
  Double_t pl0=(px0*pxl + py0*pyl + pz0*pzl)/pl;
  Double_t plb=(pxb*pxl + pyb*pyl + pzb*pzl)/pl;
  Double_t pl0s=gamma*pl0 - betagamma*e0;

  Double_t diff=2*gamma*pl0s + betagamma*des;

  return (pl0-plb-diff);
}
*/

Double_t AliCascadeVertex::ChangeMassHypothesis(Double_t &v0q, Int_t code) {
  //--------------------------------------------------------------------
  // This function changes the mass hypothesis for this cascade
  // and returns the "kinematical quality" of this hypothesis
  // together with the "quality" of associated V0 (argument v0q) 
  //--------------------------------------------------------------------
  Double_t nmass=0.13957, pmass=0.93827, ps0=0.101; 
  Double_t bmass=0.13957, mass =1.3213,  ps =0.139;

  fPdgCode=code;

  switch (code) {
  case 213: 
       bmass=0.93827; 
       break;
  case kXiMinus:
       break;
  case kXiPlusBar:
       nmass=0.93827; pmass=0.13957; 
       break;
  case kOmegaMinus: 
       bmass=0.49368; mass=1.67245; ps=0.211;
       break;
  case kOmegaPlusBar: 
       nmass=0.93827; pmass=0.13957; 
       bmass=0.49368; mass=1.67245; ps=0.211;
       break;
  default:
       cerr<<"AliCascadeVertex::ChangeMassHypothesis: ";
       cerr<<"Invalide PDG code !  Assuming XiMinus's...\n";
       fPdgCode=kXiMinus;
    break;
  }

  Double_t pxn=fV0mom[0][0], pyn=fV0mom[0][1], pzn=fV0mom[0][2];
  Double_t pxp=fV0mom[1][0], pyp=fV0mom[1][1], pzp=fV0mom[1][2];
  Double_t px0=pxn+pxp, py0=pyn+pyp, pz0=pzn+pzp;
  Double_t p0=TMath::Sqrt(px0*px0 + py0*py0 + pz0*pz0);

  Double_t e0=TMath::Sqrt(1.11568*1.11568 + p0*p0);
  Double_t beta0=p0/e0;
  Double_t pln=(pxn*px0 + pyn*py0 + pzn*pz0)/p0;
  Double_t plp=(pxp*px0 + pyp*py0 + pzp*pz0)/p0;
  Double_t pt2=pxp*pxp + pyp*pyp + pzp*pzp - plp*plp;

  Double_t a=(plp-pln)/(plp+pln);
  a -= (pmass*pmass-nmass*nmass)/(1.11568*1.11568);
  a = 0.25*beta0*beta0*1.11568*1.11568*a*a + pt2;


  v0q=a - ps0*ps0;


  Double_t pxb=fBachMom[0], pyb=fBachMom[1], pzb=fBachMom[2]; 

  Double_t eb=TMath::Sqrt(bmass*bmass + pxb*pxb + pyb*pyb + pzb*pzb);
  Double_t pxl=px0+pxb, pyl=py0+pyb, pzl=pz0+pzb;
  Double_t pl=TMath::Sqrt(pxl*pxl + pyl*pyl + pzl*pzl);
  
  fEffMass=TMath::Sqrt((e0+eb)*(e0+eb) - pl*pl);

  Double_t beta=pl/(e0+eb);
  Double_t pl0=(px0*pxl + py0*pyl + pz0*pzl)/pl;
  Double_t plb=(pxb*pxl + pyb*pyl + pzb*pzl)/pl;
  pt2=p0*p0 - pl0*pl0;

  a=(pl0-plb)/(pl0+plb);
  a -= (1.11568*1.11568-bmass*bmass)/(mass*mass);
  a = 0.25*beta*beta*mass*mass*a*a + pt2;

  return (a - ps*ps);
}

void 
AliCascadeVertex::GetPxPyPz(Double_t &px, Double_t &py, Double_t &pz) const {
  //--------------------------------------------------------------------
  // This function returns the cascade momentum (global)
  //--------------------------------------------------------------------
  px=fV0mom[0][0]+fV0mom[1][0]+fBachMom[0]; 
  py=fV0mom[0][1]+fV0mom[1][1]+fBachMom[1]; 
  pz=fV0mom[0][2]+fV0mom[1][2]+fBachMom[2]; 
}

void AliCascadeVertex::GetXYZ(Double_t &x, Double_t &y, Double_t &z) const {
  //--------------------------------------------------------------------
  // This function returns cascade position (global)
  //--------------------------------------------------------------------
  x=fPos[0]; 
  y=fPos[1]; 
  z=fPos[2]; 
}

Double_t AliCascadeVertex::GetD(Double_t x0, Double_t y0, Double_t z0) const {
  //--------------------------------------------------------------------
  // This function returns the cascade impact parameter
  //--------------------------------------------------------------------

  Double_t x=fPos[0],y=fPos[1],z=fPos[2];
  Double_t px=fV0mom[0][0]+fV0mom[1][0]+fBachMom[0];
  Double_t py=fV0mom[0][1]+fV0mom[1][1]+fBachMom[1];
  Double_t pz=fV0mom[0][2]+fV0mom[1][2]+fBachMom[2];

  Double_t dx=(y0-y)*pz - (z0-z)*py; 
  Double_t dy=(x0-x)*pz - (z0-z)*px;
  Double_t dz=(x0-x)*py - (y0-y)*px;
  Double_t d=TMath::Sqrt((dx*dx+dy*dy+dz*dz)/(px*px+py*py+pz*pz));

  return d;
}

